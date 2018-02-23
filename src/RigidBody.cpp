#define TETLIBRARY
#include "RigidBody.h"
#include <iostream>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Particle.h"
#include "MatrixStack.h"
#include "Program.h"
#include "GLSL.h"
#include <unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h>

using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> ETriplet;
double inf = numeric_limits<double>::infinity();

VectorXd RBState::computeForces(double h) {
	return (M * Phi + h * (PhiT * M * Phi + B));
	
}

Matrix3d RBState::vec2crossmatrix(Vector3d a) {		// repackage a vector into a cross-product matrix
	Matrix3d A;
	A << 0, -a(2), a(1),
		a(2), 0, -a(0),
		-a(1), a(0), 0;
	return A;
}

Vector3d RBState::local2world(MatrixXd E, Vector3d x) {
	// homo coor
	VectorXd xh;
	xh.resize(4);
	xh.segment<3>(0) = x;
	xh(3) = 1.0;

	Vector3d xw = (E * xh).segment<3>(0);
	return xw;
}	

void RBState::updatePosition() {
	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->xo = nodes[i]->x;
		nodes[i]->x = local2world(E, nodes[i]->x0);
	}
}

void RBState::updateTransformationMatrix(double h) {
	E = E * (h * EE).exp();
	R = E.block(0, 0, 3, 3);
	p = E.block(0, 3, 3, 1);
	setBodyForce();
}


RigidBody::RigidBody() {
	this->numRB = 2;
	this->numJoints = 0;
	this->numVars = numRB * 6 + numJoints * 6; // add Lagrange multipliers
	ynormal << 0.0, 1.0, 0.0;
	g << 0.0, -9.8, 0.0;
	I.setIdentity();

	xl.resize(numVars);
	xu.resize(numVars);
	for (int i = 0; i < numVars; i++) {
		xl(i) = -inf;
		xu(i) = inf;
	}

	in.load_ply("rectcube");
	tetrahedralize("pqz", &in, &out);
	nVerts = out.numberofpoints;
	nTriFaces = out.numberoftrifaces;
	nEdges = out.numberofedges;
	mass = nVerts * 1.0;

	// Create rigid body structs
	for (int i = 0; i < numRB; i++) {
		bodies.push_back(make_shared<RBState>());
	}

	// Initialize rigid bodies: R,p->E, V, Omega, 
	VectorXd init_v, init_w, init_p;
	init_v.resize(3 * numRB);
	init_w.resize(3 * numRB);
	init_p.resize(3 * numRB);

	init_v.setZero();
	init_w.setZero();
	init_p << 0, 0, 0, 0, 6, 0;

	MatrixXd init_R, init_E;
	init_R.resize(3 * numRB, 3);

	for (int i = 0; i < numRB; i++) {
		bodies[i]->setMass(mass);
		bodies[i]->setLinearVelocity(init_v.segment<3>(3 * i));
		bodies[i]->setAngularVelocity(init_w.segment<3>(3 * i));
		bodies[i]->setRotational(I); // use init_R for different setting
		bodies[i]->setLocalFrameOrigin(init_p.segment<3>(3 * i));
		bodies[i]->setSpatialInertiaMatrix();
		

		for (int j = 0; j < nVerts; j++) {
			auto p = make_shared<Particle>();
			bodies[i]->nodes.push_back(p);
			p->x0 << out.pointlist[3 * j + 0],
					 out.pointlist[3 * j + 1],
					 out.pointlist[3 * j + 2];
			p->v.setZero();
		}
	}

	// Initialize Joints
	Eij.resize(4, 4);
	Eij.setZero();

	Eij.block<3, 3>(0, 0).setIdentity();
	Eij(3, 3) = 1.0;
	Eij.block<3, 1>(0, 3) << 0.0, 3.0, 0.0;


	RHS.resize(numVars);
	RHS.setZero();
	sol.resize(numVars);

	// Build buffers
	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();
	
	posBuf.resize(nTriFaces * 9 * numRB);
	norBuf.resize(nTriFaces * 9 * numRB);
	eleBuf.resize(nTriFaces * 3 * numRB);

	updatePosNor();
	for (int i = 0; i < numRB * nTriFaces; i++) {
		for (int j = 0; j < 3; j++) {
			eleBuf[3 * i + j] = 3 * i + j;
		}
	}
	
}

MatrixXd RigidBody::computeAdjoint(MatrixXd E) {
	Vector3d p = E.block<3, 1>(0, 3);
	Matrix3d R = E.block<3, 3>(0, 0);
	MatrixXd Ad;
	Ad.resize(6, 6);
	Ad.setZero();
	Ad.block<3, 3>(0, 0) = R;
	Ad.block<3, 3>(3, 3) = R;
	Ad.block<3, 3>(3, 0) = vec2crossmatrix(p)*R;
	return Ad;
}


void RigidBody::step(double h) {

	shared_ptr<QuadProgMosek> program_ = make_shared <QuadProgMosek>();
	program = program_;
	program_->setParamInt(MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT);
	program_->setParamInt(MSK_IPAR_LOG, 10);
	program_->setParamInt(MSK_IPAR_LOG_FILE, 1);
	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_DFEAS, 1e-8);
	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_INFEAS, 1e-10);
	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_MU_RED, 1e-8);
	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_NEAR_REL, 1e3);
	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_PFEAS, 1e-8);
	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_REL_GAP, 1e-8);

	program->setNumberOfVariables(numVars);
	program->setLowerVariableBound(xl);
	program->setUpperVariableBound(xu);

	A_.clear();
	// Initialize M matrix
	for (int i = 0; i < numRB; i++) {
		for (int j = 0; j < 6; j++) {
			A_.push_back(ETriplet(6 * i + j, 6 * i + j, bodies[i]->M(j, j)));
		}
	}
	// Initialize G and G' matrix
	//MatrixXd Ejk = Eij.inverse() * bodies[0]->E.inverse() * bodies[1]->E;
	//MatrixXd G;
	//G.resize(6, 12);
	//G.block<6, 6>(0, 0) = computeAdjoint(Eij.inverse());
	//G.block<6, 6>(0, 6) = computeAdjoint(Ejk);

	//MatrixXd GT;
	//GT.resize(12, 6);
	//GT = G.transpose();
	////cout << G << endl;

	//for (int i = 0; i < 6; i++) {
	//	for (int j = 0; j < 12; j++) {
	//		A_.push_back(ETriplet(6*numRB+i, j, G(i,j))); // Fill in G
	//		//A_.push_back(ETriplet(j, 6*numRB+i, G(i,j))); // Fill in G'
	//		
	//	}
	//}

	//for (int i = 0; i < 12; i++) {
	//	for (int j = 0; j < 6; j++) {
	//		A_.push_back(ETriplet(i, 6*numRB + j, GT(i, j))); // Fill in G
	//	}
	//}


	A.resize(numVars, numVars);
	A.setFromTriplets(A_.begin(), A_.end());


	program->setObjectiveMatrix(A);

	for (int i = 0; i < numRB; i++) {
		RHS.segment<6>(6 * i) = bodies[i]->computeForces(h);
	}
	//RHS.segment<6>(0).setZero();

	cout << RHS << endl;
	program->setObjectiveVector(-RHS);

	bool success = program->solve();
	sol = program->getPrimalSolution();

	cout << "sol" << endl;

	for (int i = 0; i < numRB; i++) {
		bodies[i]->setAngularVelocity(sol.segment<3>(6 * i + 0));
		bodies[i]->setLinearVelocity(sol.segment<3>(6 * i + 3));
		bodies[i]->updateTransformationMatrix(h);
	}

	updatePosNor();
}

void RigidBody::init() {
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &norBufID);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size() * sizeof(unsigned int), &eleBuf[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	assert(glGetError() == GL_NO_ERROR);
}

void RigidBody::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p)const {

	glUniform3fv(p->getUniform("kdFront"), 1, Vector3f(1.0, 0.0, 0.0).data());
	glUniform3fv(p->getUniform("kdBack"), 1, Vector3f(1.0, 1.0, 0.0).data());
	MV->pushMatrix();
	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	int h_pos = p->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);

	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);

	int h_nor = p->getAttribute("aNor");
	glEnableVertexAttribArray(h_nor);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glDrawElements(GL_TRIANGLES, 3 * nTriFaces * numRB, GL_UNSIGNED_INT, (const void *)(0 * sizeof(unsigned int)));

	glDisableVertexAttribArray(h_nor);
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();
}

void RigidBody::updatePosNor() {
	// update 
	for (int i = 0; i < numRB; i++) {
		bodies[i]->updatePosition();
	}

	for (int i = 0; i < numRB; i++) {
		for (int iface = 0; iface < nTriFaces; iface++) {
			Vector3d p1 = bodies[i]->nodes[out.trifacelist[3 * iface + 0]]->x;
			Vector3d p2 = bodies[i]->nodes[out.trifacelist[3 * iface + 1]]->x;
			Vector3d p3 = bodies[i]->nodes[out.trifacelist[3 * iface + 2]]->x;

			//Position
			Vector3d e1 = p2 - p1;
			Vector3d e2 = p3 - p1;
			Vector3d normal = e1.cross(e2);
			normal.normalize();

			for (int idx = 0; idx < 3; idx++) {
				posBuf[nTriFaces * 9 * i + 9 * iface + 0 + idx] = p1(idx);
				posBuf[nTriFaces * 9 * i + 9 * iface + 3 + idx] = p2(idx);
				posBuf[nTriFaces * 9 * i + 9 * iface + 6 + idx] = p3(idx);
				norBuf[nTriFaces * 9 * i + 9 * iface + 0 + idx] = normal(idx);
				norBuf[nTriFaces * 9 * i + 9 * iface + 3 + idx] = normal(idx);
				norBuf[nTriFaces * 9 * i + 9 * iface + 6 + idx] = normal(idx);
			}
		}
	}	
}

Matrix3d RigidBody::vec2crossmatrix(Vector3d a) {
	Matrix3d A;
	A << 0, -a(2), a(1),
		a(2), 0, -a(0),
		-a(1), a(0), 0;
	return A;
}