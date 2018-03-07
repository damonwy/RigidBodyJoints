#define TETLIBRARY
#include "RigidBody.h"

#include <iostream>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "GLSL.h"
#include "Joint.h"
#include "MatrixStack.h"
#include "odeBoxBox.h"
#include "MatlabDebug.h"
#include "Particle.h"
#include "Program.h"
#include "RBState.h"
#include "Spring.h"

#include <unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h>

using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> ETriplet;
double inf = numeric_limits<double>::infinity();

RigidBody::RigidBody() {
	this->numRB = 4;
	this->numJoints = 1;
	this->numFixed = 2;
	this->yfloor = -0.0;
	this->ynormal << 0.0, 1.0, 0.0;
	this->I.setIdentity();
	gamma.resize(3, 6);
	gamma.block<3, 3>(0, 3) = I;
	gamma_k.resize(3, 6);
	gamma_k.block<3, 3>(0, 3) = I;
	numVars = numRB * 6;
	numEqualities = 0;
	numInequalities = 0;

	joint_forces.resize(numVars);
	joint_forces.setZero();

	g << 0.0, -9.8, 0.0;

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

	// Initialize rigid bodies
	VectorXd init_v, init_w, init_p;
	init_v.resize(3 * numRB);
	init_w.resize(3 * numRB);
	init_p.resize(3 * numRB);

	init_v.setZero();
	init_v.segment<3>(3 * 1) << -10.0, 0.0, 0.0;
	
	init_w.setZero();
	//init_w << 0, 0, 0, 0, 0, 0;
	init_p << 
		0.0, 12.0, 0.0,
		0.0, 6.0, 0.0,
		6.0, 3.0, 1.0,
		15.0, 7.5, 2.0;		

	MatrixXd init_R, init_E;
	init_R.resize(3 * numRB, 3);
	init_R.block<3, 3>(0, 0) = I;
	for (int i = 0; i < numRB; i++) {
		init_R.block<3, 3>(i * 3, 0) = I;
	}
	
	/*init_R.block<3, 3>(9, 0) << 0, -1, 0,
		1, 0, 0,
		0, 0, 1;
	init_R.block<3, 3>(12, 0) << 0, -1, 0,
		1, 0, 0,
		0, 0, 1;
	init_R.block<3, 3>(15, 0) << 0, -1, 0,
		1, 0, 0,
		0, 0, 1;
	*/
	init_R.block<3, 3>(9, 0) << 0, -1, 0,
		1, 0, 0,
		0, 0, 1;

	init_R.block<3, 3>(6, 0) << sqrt(2)/2.0, 0, -sqrt(2) / 2.0,
		0, 1, 0,
		sqrt(2) / 2.0, 0, sqrt(2) / 2.0;

	init_fixed_rb.resize(numFixed);
	init_fixed_rb << 0, 3;
	convec.resize(6);
	conveck.resize(6);
	Vector3d dimensions;
	dimensions << 2 * abs(out.pointlist[0]), 2 * abs(out.pointlist[1]), 2 * abs(out.pointlist[2]);


	// Initialize Rigid Bodies
	for (int i = 0; i < numRB; i++) {
		
		bodies[i]->setDimensions(dimensions);
		bodies[i]->setMass(mass);
		bodies[i]->setLinearVelocity(init_v.segment<3>(3 * i));
		bodies[i]->setAngularVelocity(init_w.segment<3>(3 * i));
		bodies[i]->setRotational(init_R.block<3, 3>(3 * i, 0));
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
	for (int i = 0; i < numJoints; i++) {
		auto jt = make_shared<Joint>();
		joints.push_back(jt);
		// specify joint type here
		jt->type = HINGE_JOINT;
		if (jt->type == HINGE_JOINT) {
			jt->hinge_type = HINGE_JOINT_Z;
		}

		jt->i = i;
		jt->k = i + 1;
		jt->index = i;
		jt->Eij.block<3, 3>(0, 0) = I;
		jt->Eij.block<3, 1>(0, 3) << 0.0, -3.0, 0.0;
	}
	
	// Initialize Springs

	double stiffness = 1e2;
	springs.push_back(createSpring(1, 3, 0, 3, bodies, stiffness));

	// Compute the number of equality constraints
	for (int i = 0; i < numJoints; i++) {
		numEqualities += joints[i]->type;
	}
	numEqualities += numFixed * 6;

	// Initialize QP b vector and solution vector
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

	// Initialize A matrix
	A_.clear();
	for (int i = 0; i < numRB; i++) {
		for (int j = 0; j < 6; j++) {
			A_.push_back(ETriplet(6 * i + j, 6 * i + j, bodies[i]->M(j, j)));
		}
	}
	A.resize(numVars, numVars);
	A.setFromTriplets(A_.begin(), A_.end());
	program->setObjectiveMatrix(A);

	//sparse_to_file_as_dense(A, "A");

	// Initialize b vector

	// Set spring forces
	for (int i = 0; i < springs.size(); i++) {
		auto b0 = bodies[springs[i]->i];
		auto b1 = bodies[springs[i]->k];

		auto n0 = b0->nodes[springs[i]->in];
		auto n1 = b1->nodes[springs[i]->kn];

		Vector3d d = (n0->x - n1->x).normalized();
		double l = (n0->x - n1->x).norm();
		double f = springs[i]->E * (l / springs[i]->L - 1.0);

		Vector3d f0 = -f * d;
		Vector3d f1 = -f0;

		gamma.block<3, 3>(0, 0) = vec2crossmatrix(n0->x0).transpose();
		VectorXd wrench0 = (b0->R * gamma).transpose() * f0;
		gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(n1->x0).transpose();
		VectorXd wrench1 = (b1->R * gamma_k).transpose() * f1;

		b0->setBodyForce(wrench0);
		b1->setBodyForce(wrench1);
	}


	for (int i = 0; i < numRB; i++) {
		RHS.segment<6>(6 * i) = bodies[i]->computeForces(h);
	}
	program->setObjectiveVector(-RHS);

	MatrixXd Gi, Gk;
	Gi.resize(6, 6);
	Gk.resize(6, 6);

	G_.clear();

	// Initialize G matrix
	int currentrow = 0;

	// Push back joints constraints
	for (int i = 0; i < numJoints; i++) {
		Gi = joints[i]->computeAdjoint(joints[i]->Eij.inverse());
		joints[i]->computeEjk(bodies);
		Gk = -joints[i]->computeAdjoint(joints[i]->Ejk); // careful about the sign!

		if (joints[i]->type == BALL_JOINT) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 6; k++) {
					G_.push_back(ETriplet(currentrow , 6 * joints[i]->i + k, Gi(j + 3, k)));
					G_.push_back(ETriplet(currentrow , 6 * joints[i]->k + k, Gk(j + 3, k)));
				}
				currentrow += 1;
			}
		}

		if (joints[i]->type == HINGE_JOINT) {
			
			for (int j = 0; j < 6; j++) {
				if (j != joints[i]->hinge_type) {
					for (int k = 0; k < 6; k++) {
						G_.push_back(ETriplet(currentrow, 6 * joints[i]->i + k, Gi(j, k)));
						G_.push_back(ETriplet(currentrow, 6 * joints[i]->k + k, Gk(j, k)));
					}
					currentrow += 1;
				}
			}
		}
	}

	// Push back fixed constraints
	for (int i = 0; i < numFixed; i++) {
		for (int j = 0; j < 6; j++) {
			G_.push_back(ETriplet(currentrow + j, 6 * init_fixed_rb(i) + j, 1));
		}
		currentrow += 6;
	}

	GG.resize(numEqualities, numVars);
	GG.setFromTriplets(G_.begin(), G_.end());

	// Initialize equality vector
	equalvec.resize(numEqualities);
	equalvec.setZero();

	program->setNumberOfEqualities(numEqualities);
	program->setEqualityMatrix(GG);
	program->setEqualityVector(equalvec);

	//sparse_to_file_as_dense(GG, "GG");

	bool success = program->solve();
	sol = program->getPrimalSolution();
	joint_forces = -GG.transpose() * program->getDualEquality()/h;
	
	//vec_to_file(sol, "sol");


	for (int i = 0; i < numRB; i++) {
		bodies[i]->setAngularVelocity(sol.segment<3>(6 * i + 0));
		bodies[i]->setLinearVelocity(sol.segment<3>(6 * i + 3));
		bodies[i]->computeTempE(h);
	}

	// Update node positions to detect collisions with floor
	numInequalities = 0;
	numColFloor = 0;
	for (int i = 0; i < numRB; i++) {
		auto b = bodies[i];

		for (int j = 0; j < b->nodes.size(); j++) {
			auto n = b->nodes[j];
			n->xo = n->x;
			n->x = local2world(b->Etemp, n->x0);
			if (n->x(1) < yfloor) {
				colList.push_back(i);
				colList.push_back(j);
				numColFloor += 1;
			}
		}
	}
	numInequalities += numColFloor;
	numColBoxBox = 0;
	contacts.clear();

	// Detect collisions between rigid bodies
	for (int i = 0; i < numRB; i++) {

		if (i != 2) { //TODO: modify to the index of all the rigid bodies
			auto b = bodies[i];
			auto cts = make_shared<Contacts>(odeBoxBox(b->E, b->dimensions, bodies[2]->E, bodies[2]->dimensions));  // TODO : modify to the index of all the rigid bodies
			if (cts->count > 0) {
				cts->i = i;
				cts->k = 2; //TODO
				contacts.push_back(cts);
				numColBoxBox += cts->count;
			}
		}
	}

	numInequalities += numColBoxBox;
	
	// Add Inequalities constraints
	if (numInequalities == 0) {
		for (int i = 0; i < numRB; i++) {
			bodies[i]->updateTransformationMatrix(h);
		}
	}
	else {
		VectorXd inequalvec;

		// Initialize inequality vector
		inequalvec.resize(numInequalities);
		inequalvec.setZero();
		C.resize(numInequalities, numVars);
		program->setNumberOfInequalities(numInequalities);

		int currentrow = 0;

		// Floor collisions
		for (int i = 0; i < numColFloor; i++) {

			int ib = colList[2 * i + 0]; // the index of rigid body 
			int in = colList[2 * i + 1]; // the index of node 

			gamma.block<3, 3>(0, 0) = vec2crossmatrix(bodies[ib]->nodes[in]->x0).transpose();
			// R and p were not updated yet... were still the values before the collsion 

			//Matrix3d Rtemp = bodies[ib]->Etemp.block<3, 3>(0, 0);
			//convec = ynormal.transpose() * Rtemp * gamma;

			convec = ynormal.transpose() * bodies[ib]->R * gamma; // 1x6 vector that should be in the constraint matrix

			for (int t = 0; t < 6; t++) {
				C_.push_back(ETriplet(i, 6 * ib + t, -convec(t)));

			}
		}

		currentrow += numColFloor;
		
		// Push back constraints with rigid bodies
		for (int i = 0; i < contacts.size(); i++) {
			auto cts = contacts[i];
			
			for (int j = 0; j < cts->count; j++) {
				// Change the world position to local position
				Vector3d x0 = cts->positions[j];			// collision point in rb i in world frame
				Vector3d xi = local2world(bodies[cts->i]->Etemp.inverse(), x0);
				Vector3d xk = local2world(bodies[cts->k]->Etemp.inverse(), x0);


				gamma.block(0, 0, 3, 3) = vec2crossmatrix(xi).transpose();

				Matrix3d Ri = bodies[cts->i]->Etemp.block<3, 3>(0, 0);
				//convec = cts->normal.transpose() * bodies[cts->i]->R * gamma;
				
				convec = cts->normal.transpose() * Ri * gamma;
				Matrix3d Rk = bodies[cts->k]->Etemp.block<3, 3>(0, 0);

				gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(xk).transpose();
				//conveck = cts->normal.transpose() * bodies[cts->k]->R * gamma_k;

				conveck = cts->normal.transpose() * Rk * gamma_k;

			
				for (int t = 0; t < 6; t++) {
					C_.push_back(ETriplet(currentrow, 6 * cts->i + t, convec(t)));
					C_.push_back(ETriplet(currentrow, 6 * cts->k + t, -conveck(t)));
				}
				currentrow += 1;
				
			}
		}

		C.setFromTriplets(C_.begin(), C_.end());

		program->setInequalityMatrix(C);
		program->setInequalityVector(inequalvec);

		//sparse_to_file_as_dense(C, "C");
		//sparse_to_file_as_dense(GG, "GG");

		// Initialize equality vector
		equalvec.resize(numEqualities);
		equalvec.setZero();

		bool success = program->solve();
		sol = program->getPrimalSolution();
		//vec_to_file(sol, "sol");

		joint_forces = -GG.transpose() * program->getDualEquality() / h;

		for (int i = 0; i < numRB; i++) {
			bodies[i]->setAngularVelocity(sol.segment<3>(6 * i + 0));
			bodies[i]->setLinearVelocity(sol.segment<3>(6 * i + 3));
			bodies[i]->updateTransformationMatrix(h);
		}

		C_.clear();
		colList.clear();
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

void RigidBody::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p, const shared_ptr<Program> p2, shared_ptr<MatrixStack> P)const {

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
	
	for (int i = 0; i < numRB; i++) {
		glUniform3fv(p->getUniform("kdBack"), 1, Vector3f((150 + 20 * i) / 255.0, (100 + 10*i) / 255.0, (150.0 + 10 * i) / 255.0).data());
		glDrawElements(GL_TRIANGLES, 3 * nTriFaces, GL_UNSIGNED_INT, (const void *)(3 * nTriFaces * i * sizeof(unsigned int)));
	}
	
	glDisableVertexAttribArray(h_nor);
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();

	p->unbind();

	p2->bind();
	glUniformMatrix4fv(p2->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(p2->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->pushMatrix();
	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	glPointSize(5);

	glBegin(GL_POINTS);
	for (int i = 0; i < nVerts; i++) {
		switch (i) {
		case 0:
			glColor3f(1.0, 0.0, 0.0); // red
			break;
		case 1:
			glColor3f(1.0, 1.0, 0.0); // yellow
			break;
		case 2:
			glColor3f(0.0, 0.0, 1.0);  // blue
			break;
		case 3:
			glColor3f(0.0, 1.0, 0.0); // lime
			break;
		case 4:
			glColor3f(128.0/ 255.0, 0.0, 0.0);// maroon
			break;
		case 5:
			glColor3f(128.0 / 255.0, (128.0) / 255.0, 0.0);// olive
			break;
		case 6:
			glColor3f(0.0, 0.0, 128.0 / 255.0); // navy
			break;
		case 7:
			glColor3f(0.0, 128.0 / 255.0 ,0.0); // green
			break;
		}
		for (int j = 0; j < numRB; j++) {
			glVertex3f(float(bodies[j]->nodes[i]->x(0)), float(bodies[j]->nodes[i]->x(1)), float(bodies[j]->nodes[i]->x(2)));
		}
	}
	glEnd();

	glColor3f(1.0, 1.0, 0.0); // yellow
	if (numColBoxBox != 0) {
		glLineWidth(4);
		for (int i = 0; i < contacts.size(); i++) {
			for (int j = 0; j < contacts[i]->count; j++){
				Vector3d p1 = contacts[i]->positions[j];
				Vector3d p2 = 1.5 * contacts[i]->normal + contacts[i]->positions[j];
				glBegin(GL_LINES);
				glVertex3f(float(p1(0)), float(p1(1)), float(p1(2)));
				glVertex3f(float(p2(0)), float(p2(1)), float(p2(2)));
				glEnd();
				glPointSize(10);
				glBegin(GL_POINTS);
				glVertex3f(float(p1(0)), float(p1(1)), float(p1(2)));
				glEnd();
			}
		}
	}

	glColor3f(0.0, 0.0, 0.0); // black
	glLineWidth(5);
	for (int i = 0; i < springs.size(); i++) {
		Vector3d p0 = bodies[springs[i]->i]->nodes[springs[i]->in]->x;
		Vector3d p1 = bodies[springs[i]->k]->nodes[springs[i]->kn]->x;
		glBegin(GL_LINES);
		glVertex3f(float(p0(0)), float(p0(1)), float(p0(2)));
		glVertex3f(float(p1(0)), float(p1(1)), float(p1(2)));
		glEnd();
	}

	MV->popMatrix();
	p2->unbind();
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

Vector3d RigidBody::local2world(MatrixXd E, Vector3d x) {
	// homo coor
	VectorXd xh;
	xh.resize(4);
	xh.segment<3>(0) = x;
	xh(3) = 1.0;

	Vector3d xw = (E * xh).segment<3>(0);
	return xw;
}

Vector3d RigidBody::world2local(MatrixXd E, Vector3d x) {
	// homo coor
	VectorXd xh;
	xh.resize(4);
	xh.segment<3>(0) = x;
	xh(3) = 1.0;

	Vector3d xw = (E * xh).segment<3>(0);
	return xw;
}

shared_ptr<Spring> RigidBody::createSpring(int _i, int _k, int _in, int _kn, vector < shared_ptr<RBState> > bodies, double E)
{	
	auto s = make_shared<Spring>(bodies[_i]->nodes[_in], bodies[_k]->nodes[_kn]);
	s->i = _i;
	s->k = _k;
	s->in = _in;
	s->kn = _kn;
	s->E = E;
	Vector3d x0 = bodies[_i]->nodes[_in]->x0;
	Vector3d x1 = bodies[_k]->nodes[_kn]->x0;
	Vector3d dx = x1 - x0;
	s->L = dx.norm();
	return s;
}