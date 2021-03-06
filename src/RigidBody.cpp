#define TETLIBRARY
#include "RigidBody.h"

#include <iostream>
#include <fstream>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>


#include <math.h>       /* isnan, sqrt */

#include "GLSL.h"
#include "Joint.h"
#include "MatrixStack.h"
#include "odeBoxBox.h"
#include "MatlabDebug.h"
#include "Particle.h"
#include "Program.h"
#include "RBState.h"
#include "Spring.h"
#include "WrapSphere.h"
#include "WrapDoubleCylinder.h"
#include "WrapCylinder.h"

#include "Shape.h"
#include "Cylinder.h"
#include "DoubleCylinder.h"
#include "Helper.h"
#include "JsonEigen.h"

// for convenience
using json = nlohmann::json;

#include <unsupported/Eigen/MatrixFunctions> // TODO: avoid using this later, write a func instead

using namespace std;
using namespace Eigen;
using namespace Helper;

typedef Eigen::Triplet<double> ETriplet;
double inf = numeric_limits<double>::infinity();

RigidBody::RigidBody():
	numEqualities(0),
	numInequalities(0),
	steps(0),
	gamma(3, 6),
	gamma_k(3, 6),
	convec(6),
	conveck(6)
{
	//read a JSON file
	ifstream i("input.json");
	i >> js;
	i.close();

	// change when model changes
	this->numRB = js["numRB"];
	this->numJoints = js["numJoints"];
	this->numSprings = js["numSprings"];
	this->numCylinders = js["numCylinders"];
	this->numDoubleCylinders = js["numDoubleCylinders"];
	this->numWrapPoints = js["numWrapPoints"];
	this->numFinitePoints = js["numFinitePoints"];
	this->numFixed = js["numFixed"];
	this->yfloor = js["yfloor"];
	this->isBoxBoxCol = js["isBoxBoxCol"];
	this->isFloorCol = js["isFloorCol"];
	this->isFEM = js["isFEM"];
	this->debug_i = js["debug_i"];
	this->epsilon = js["epsilon"];
	this->muscle_mass = js["muscle_mass"];
	this->stiffness = js["stiffness"];
	this->ynormal << 0.0, 1.0, 0.0;
	this->g << 0.0, -9.8, 0.0;
	this->I.setIdentity();
	
	gamma.block<3, 3>(0, 3) = I;
	gamma_k.block<3, 3>(0, 3) = I;

	if (numRB != 0) {
		initAfterNumRB();
	}

	initRBShape();

	if (numFixed != 0) {
		init_fixed_rb.resize(numFixed);
		// Compute the number of equality constraints (part 2: fixed)
		numEqualities += numFixed * 6;
	}

	if (numJoints != 0) {
		initJoints();
	}

	if (numCylinders != 0) {
		initAfterNumCylinders();
	}

	if (numDoubleCylinders != 0) {
		initAfterNumDoubleCylinders();
	}

	test_pts.resize(3, numFinitePoints);
	test_pts_pert.resize(3, numFinitePoints);
	test_pts.setZero();
	test_pts_pert.setZero();
	
	//init_v.segment<3>(3 * 1) << -10.0, 0.0, 0.0;
	/*init_p << 0.0, 14.0, 0.0,
		3.0, 11.0, 0.0,
		3.0, 5.0, 0.0;*/

	//init_p << 0.0, 14.0, 0.0,
	//	3.0, 11.0, 0.0,
	//	6.0, 8.0, 0.0;

	init_p << 0.0, 14.0, 0.0,
		6.0, 14.0, 0.0;
		
	init_R.block<3, 3>(0, 0) << 0, -1, 0,
		1, 0, 0,
		0, 0, 1;
	//
	init_R.block<3, 3>(3, 0) << 0, -1, 0,
		1, 0, 0,
	0, 0, 1;

	init_fixed_rb << 0;

	/*init_cyl_P << 1.0, 0.0, 0.0;
	init_cyl_S << 1.0, -2.0, 0.0;
	init_cyl_O << 1.5, -1.5, -0.2;
	init_cyl_Z << 0.0, 0.0, 1.0;
	init_cyl_r << 0.4;
	init_cyl_O_rb << 1;
	init_cyl_P_rb << 0;
	init_cyl_S_rb << 1;*/
	
	//nlohmann::from_json(j["init_fixed_rb"], init_fixed_rb);

	//-------------------------------------------------------------
	
	initRBs();

	if (numCylinders != 0) {
		initCylinder();
	}

	if (numDoubleCylinders != 0) {
		initDoubleCylinder();
	}

	if (numSprings != 0) {
		// init after initRBs
		initSprings(stiffness);
	}
	initBuffers();
}

void RigidBody::initCylinder() {
	for (int i = 0; i < numCylinders; i++) {
		auto P = make_shared<Particle>();
		Ps.push_back(P);
		P->x0 = init_cyl_P.segment<3>(3 * i);
		P->rb_id = init_cyl_P_rb(i);
		P->x = transform(bodies[P->rb_id]->E, P->x0);

		auto S = make_shared<Particle>();
		Ss.push_back(S);
		S->x0 = init_cyl_S.segment<3>(3 * i);
		S->rb_id = init_cyl_S_rb(i);
		S->x = transform(bodies[S->rb_id]->E, S->x0);

		auto O = make_shared<Particle>();
		Os.push_back(O);
		O->x0 = init_cyl_O.segment<3>(3 * i);
		O->rb_id = init_cyl_O_rb(i);
		O->x = transform(bodies[O->rb_id]->E, O->x0);

		auto cylinder = make_shared<Cylinder>(P, S, O);
		cylinders.push_back(cylinder);
		cylinder->r = init_cyl_r(i);
		cylinder->Z = init_cyl_Z.segment<3>(3 * i);
		cylinder->mass = muscle_mass;
	}
}

void RigidBody::initDoubleCylinder() {
	for (int i = 0; i < numDoubleCylinders; i++) {
		auto P = make_shared<Particle>();
		dPs.push_back(P);
		P->x0 = init_dcyl_P.segment<3>(3 * i);
		P->rb_id = init_dcyl_P_rb(i);
		P->x = transform(bodies[P->rb_id]->E, P->x0);

		auto S = make_shared<Particle>();
		dSs.push_back(S);
		S->x0 = init_dcyl_S.segment<3>(3 * i);
		S->rb_id = init_dcyl_S_rb(i);
		S->x = transform(bodies[S->rb_id]->E, S->x0);

		auto U = make_shared<Particle>();
		dUs.push_back(U);
		U->x0 = init_dcyl_U.segment<3>(3 * i);
		U->rb_id = init_dcyl_U_rb(i);
		U->x = transform(bodies[U->rb_id]->E, U->x0);

		auto V = make_shared<Particle>();
		dVs.push_back(V);
		V->x0 = init_dcyl_V.segment<3>(3 * i);
		V->rb_id = init_dcyl_V_rb(i);
		V->x = transform(bodies[V->rb_id]->E, V->x0);

		auto doublecylinder = make_shared<DoubleCylinder>(P, S, U, V);

		doublecylinders.push_back(doublecylinder);
		doublecylinder->Ur = init_dcyl_Ur(i);
		doublecylinder->Zu = init_dcyl_Uz.segment<3>(3 * i);
		doublecylinder->Vr = init_dcyl_Vr(i);
		doublecylinder->Zv = init_dcyl_Vz.segment<3>(3 * i);
	}
}

void RigidBody::initRBs() {
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

void RigidBody::initAfterNumRB() {
	numVars = numRB * 6;
	// Initialize QP b vector and solution vector

	RHS.resize(numVars);
	RHS.setZero();
	sol.resize(numVars);
	xl.resize(numVars);
	xu.resize(numVars);
	for (int i = 0; i < numVars; i++) {
		xl(i) = -inf;
		xu(i) = inf;
	}

	// Create rigid body structs
	for (int i = 0; i < numRB; i++) {
		bodies.push_back(make_shared<RBState>());
	}

	// Initialize rigid bodies
	init_v.resize(3 * numRB);
	init_w.resize(3 * numRB);
	init_p.resize(3 * numRB);
	init_R.resize(3 * numRB, 3);

	init_v.setZero();
	init_w.setZero();
	init_p.setZero();
	for (int i = 0; i < numRB; i++) {
		init_R.block<3, 3>(i * 3, 0) = I;
	}

	Adense.resize(6 * numRB, 6 * numRB);
	Adense.setZero();
}

void RigidBody::initAfterNumCylinders() {

	// init cylinders wrapper
	init_cyl_P.resize(numCylinders * 3);
	init_cyl_S.resize(numCylinders * 3);
	init_cyl_Z.resize(numCylinders * 3);
	init_cyl_O.resize(numCylinders * 3);
	init_cyl_r.resize(numCylinders);
	init_cyl_O_rb.resize(numCylinders);
	init_cyl_P_rb.resize(numCylinders);
	init_cyl_S_rb.resize(numCylinders);
	wpc.resize(3 * numCylinders, numWrapPoints + 1);
	wpc_stat.resize(numCylinders);
	wpc_length.resize(numCylinders);

	wpc_pert.resize(3 * 12, numWrapPoints + 1 + 2);
	wpc_stat_pert.resize(12);
	wpc_length_pert.resize(12);
}

void RigidBody::initAfterNumDoubleCylinders() {
	// init doublecylinders wrapper
	init_dcyl_P.resize(numDoubleCylinders * 3);
	init_dcyl_P_rb.resize(numDoubleCylinders);

	init_dcyl_S.resize(numDoubleCylinders * 3);
	init_dcyl_S_rb.resize(numDoubleCylinders);

	init_dcyl_Uz.resize(numDoubleCylinders * 3);
	init_dcyl_U.resize(numDoubleCylinders * 3);
	init_dcyl_Ur.resize(numDoubleCylinders);
	init_dcyl_U_rb.resize(numDoubleCylinders);

	init_dcyl_Vz.resize(numDoubleCylinders * 3);
	init_dcyl_V.resize(numDoubleCylinders * 3);
	init_dcyl_Vr.resize(numDoubleCylinders);
	init_dcyl_V_rb.resize(numDoubleCylinders);


	wpdc.resize(3 * numDoubleCylinders, 3 * numWrapPoints + 1);
	wpdc_stat.resize(numDoubleCylinders);
	wpdc_length.resize(numDoubleCylinders);
}

void RigidBody::initRBShape() {
	// init cubes
	in.load_ply("rectcube");
	tetrahedralize("pqz", &in, &out);
	this->nVerts = out.numberofpoints;
	this->nTriFaces = out.numberoftrifaces;
	this->nEdges = out.numberofedges;
	mass = nVerts * 1.0;
	dimensions << 2 * abs(out.pointlist[0]), 2 * abs(out.pointlist[1]), 2 * abs(out.pointlist[2]);
}

void RigidBody::initJoints() {
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
		jt->Eij.block<3, 1>(0, 3) << 0.0, -dimensions(1)*0.5, 0.0; // use y direction 

		// Compute the number of equality constraints (part 1: joints)
		numEqualities += jt->type;
	}
}

void RigidBody::initSprings(double stiffness) {
	springs.push_back(createSpring2RB(0, 1, 4, 1, bodies, stiffness, muscle_mass));
}

void RigidBody::initBuffers() {
	// Build buffers
	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();

	posBuf.resize(nTriFaces * 9 * numRB);
	norBuf.resize(nTriFaces * 9 * numRB);
	eleBuf.resize(nTriFaces * 3 * numRB);

	updatePosNor();
	updateWrapCylinders();
	updateDoubleWrapCylinders();

	for (int i = 0; i < numRB * nTriFaces; i++) {
		for (int j = 0; j < 3; j++) {
			eleBuf[3 * i + j] = 3 * i + j;
		}
	}
}

void RigidBody::computeSpringForces() {
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
		VectorXd wrench0 = (b0->R * gamma).transpose() * (f0);
		gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(n1->x0).transpose();
		VectorXd wrench1 = (b1->R * gamma_k).transpose() * (f1);

		b0->setBodyForce(wrench0);
		b1->setBodyForce(wrench1);
	}
}

void RigidBody::computeWrapCylinderForces() {
	for (int i = 0; i < numCylinders; i++) {
		auto c = cylinders[i];
		double f = c->stiffness * (c->l / c->L - 1);

		if (wpc_stat[i] == wrap) {
			// if the status is wrapped, then we need to consider four forces

			// P->C1 path:
			if (c->P->rb_id == c->O->rb_id) {
				// if they are in the same rigid body, forces cancel out
				c->fp.setZero();
				c->fpc.setZero();
			}
			else {
				// if they are in different rigid body, compute the forces and apply them
				Vector3d fp = -f * c->pdir; // apply to the rb that p is on
				Vector3d fpc = -fp; // apply to the rb that cylinder is on

				c->fp = fp;
				c->fpc = fpc;

				gamma.block<3, 3>(0, 0) = vec2crossmatrix(c->P->x0).transpose();
				VectorXd wrench0 = (bodies[c->P->rb_id]->R * gamma).transpose() * fp;
				bodies[c->P->rb_id]->setBodyForce(wrench0);

				// compute the local coordinate of the contact point in cylinder	
				Vector3d c1 = transform(bodies[c->O->rb_id]->E.inverse(), wpc.block<3, 1>(3 * i, numWrapPoints));

				gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(c1).transpose();
				VectorXd wrench1 = (bodies[c->O->rb_id]->R * gamma_k).transpose() * fpc;
				bodies[c->O->rb_id]->setBodyForce(wrench1);
			}

			// C2->S path:
			if (c->S->rb_id == c->O->rb_id) {
				c->fs.setZero();
				c->fsc.setZero();
			}
			else {

				Vector3d fs = -f * c->sdir; // apply to the rb that s is on
				Vector3d fsc = -fs; // apply to the rb that cylinder is on

				c->fs = fs;
				c->fsc = fsc;

				gamma.block<3, 3>(0, 0) = vec2crossmatrix(c->S->x0).transpose();
				VectorXd wrench0 = (bodies[c->S->rb_id]->R * gamma).transpose() * fs;
				bodies[c->S->rb_id]->setBodyForce(wrench0);

				// compute the local coordinate of the contact point in cylinder	
				Vector3d c2 = transform(bodies[c->O->rb_id]->E.inverse(), wpc.block<3, 1>(3 * i, 0));

				gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(c2).transpose();
				VectorXd wrench1 = (bodies[c->O->rb_id]->R * gamma_k).transpose() * fsc;
				bodies[c->O->rb_id]->setBodyForce(wrench1);
			}
		}
		else {
			// if the status is not wrapped, only consider two forces, applied on P, S seperately
			Vector3d fp = -f * c->pdir; // apply to the rb that p is on
			c->fp = fp;
			Vector3d fs = -f * c->sdir; // apply to the rb that s is on
			c->fs = fs;
			c->fpc.setZero();
			c->fsc.setZero();

			gamma.block<3, 3>(0, 0) = vec2crossmatrix(c->P->x0).transpose();
			VectorXd wrench0 = (bodies[c->P->rb_id]->R * gamma).transpose() * fp;

			gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(c->S->x0).transpose();
			VectorXd wrench1 = (bodies[c->S->rb_id]->R * gamma_k).transpose() * fs;

			bodies[c->P->rb_id]->setBodyForce(wrench0);
			bodies[c->S->rb_id]->setBodyForce(wrench1);
		}
	}
}

void RigidBody::computeWCGravityForces() {
	for (int i = 0; i < numCylinders; i++) {
		auto c = cylinders[i];

	}

}

void RigidBody::computeWrapDoubleCylinderForces() {
	for (int i = 0; i < numDoubleCylinders; i++) {
		auto dc = doublecylinders[i];
		double f = dc->stiffness * (dc->l / dc->L - 1);

		if (wpdc_stat[i] == 1) {
			// if the status is wrapped, then we need to consider six forces

			// P->U1 path:
			if (dc->P->rb_id == dc->U->rb_id) {
				// if they are in the same rigid body, forces cancel out
				dc->fpu.setZero();
				dc->fup.setZero();
			}
			else {
				// if they are in different rigid body, compute the forces and apply them
				Vector3d fpu = -f * dc->pdir;	// apply to the rb that p is on
				Vector3d fup = -fpu;			// apply to the rb that cylinder U is on

												// save this for drawing purpose
				dc->fpu = fpu;
				dc->fup = fup;

				gamma.block<3, 3>(0, 0) = vec2crossmatrix(dc->P->x0).transpose();
				VectorXd wrench0 = (bodies[dc->P->rb_id]->R * gamma).transpose() * fpu;
				bodies[dc->P->rb_id]->setBodyForce(wrench0);

				// compute the local coordinate of the contact point U1 in U cylinder	
				Vector3d u1 = transform(bodies[dc->U->rb_id]->E.inverse(), dc->u1);

				gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(u1).transpose();
				VectorXd wrench1 = (bodies[dc->U->rb_id]->R * gamma_k).transpose() * fup;
				bodies[dc->U->rb_id]->setBodyForce(wrench1);
			}

			// U2->V1 path:
			if (dc->U->rb_id == dc->V->rb_id) {
				dc->fuv.setZero();
				dc->fvu.setZero();
			}
			else {

				Vector3d fuv = -f * dc->uvdir;
				Vector3d fvu = -fuv;

				dc->fuv = fuv;
				dc->fvu = fvu;
				// compute the local coordinate of the contact point u2 in cylinder	U
				Vector3d u2 = transform(bodies[dc->U->rb_id]->E.inverse(), dc->u2);
				gamma.block<3, 3>(0, 0) = vec2crossmatrix(u2).transpose();
				VectorXd wrench0 = (bodies[dc->U->rb_id]->R * gamma).transpose() * fuv;
				bodies[dc->U->rb_id]->setBodyForce(wrench0);

				// compute the local coordinate of the contact point v1 in cylinder	V
				Vector3d v1 = transform(bodies[dc->V->rb_id]->E.inverse(), dc->v1);

				gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(v1).transpose();
				VectorXd wrench1 = (bodies[dc->V->rb_id]->R * gamma_k).transpose() * fvu;
				bodies[dc->V->rb_id]->setBodyForce(wrench1);
			}

			// V2->S path:
			if (dc->V->rb_id == dc->S->rb_id) {
				dc->fvs.setZero();
				dc->fsv.setZero();
			}
			else {

				Vector3d fvs = -f * (-dc->uvdir);
				Vector3d fsv = -fvs;
				dc->fvs = fvs;
				dc->fsv = fsv;

				// compute the local coordinate of the contact point v2 in cylinder	V
				Vector3d v2 = transform(bodies[dc->V->rb_id]->E.inverse(), dc->v2);
				gamma.block<3, 3>(0, 0) = vec2crossmatrix(v2).transpose();
				VectorXd wrench0 = (bodies[dc->V->rb_id]->R * gamma).transpose() * fvs;
				bodies[dc->V->rb_id]->setBodyForce(wrench0);

				// compute the local coordinate of S

				gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(dc->S->x0).transpose();
				VectorXd wrench1 = (bodies[dc->S->rb_id]->R * gamma_k).transpose() * fsv;
				bodies[dc->S->rb_id]->setBodyForce(wrench1);
			}
		}
	}
}

void RigidBody::setJointConstraints(int &currentrow, int type) {
	Matrix6d Gi, Gk;

	for (int i = 0; i < numJoints; i++) {
		Gi = joints[i]->computeAdjoint(joints[i]->Eij.inverse());
		if (type == 1) {
			joints[i]->computeEjk(bodies);
		}
		if (type == 0) {
			joints[i]->computeEjkTemp(bodies);
		}

		Gk = -joints[i]->computeAdjoint(joints[i]->Ejk); // careful about the sign!

		if (joints[i]->type == BALL_JOINT) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 6; k++) {
					G_.push_back(ETriplet(currentrow, 6 * joints[i]->i + k, Gi(j + 3, k)));
					G_.push_back(ETriplet(currentrow, 6 * joints[i]->k + k, Gk(j + 3, k)));
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
}

void RigidBody::setFixedConstraints(int &currentrow) {
	// Push back fixed constraints
	for (int i = 0; i < numFixed; i++) {
		for (int j = 0; j < 6; j++) {
			G_.push_back(ETriplet(currentrow + j, 6 * init_fixed_rb(i) + j, 1));
		}
		currentrow += 6;
	}
}

void RigidBody::detectFloorCol() {
	numColFloor = 0;

	for (int i = 0; i < numRB; i++) {
		auto b = bodies[i];

		for (int j = 0; j < b->nodes.size(); j++) {
			auto n = b->nodes[j];
			n->xo = n->x;
			n->x = transform(b->Etemp, n->x0);
			if (n->x(1) < yfloor) {
				colList.push_back(i);
				colList.push_back(j);
				numColFloor += 1;
			}
		}
	}
	numInequalities += numColFloor;
}

void RigidBody::detectBoxBoxCol() {
	contacts.clear();

	numColBoxBox = 0;
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
}

void RigidBody::setInequality(double h) {
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
				Vector3d xi = transform(bodies[cts->i]->Etemp.inverse(), x0);
				Vector3d xk = transform(bodies[cts->k]->Etemp.inverse(), x0);

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

		// Initialize equality vector
		equalvec.resize(numEqualities);
		equalvec.setZero();

		bool success = program->solve();
		sol = program->getPrimalSolution();

		for (int i = 0; i < numRB; i++) {
			bodies[i]->setAngularVelocity(sol.segment<3>(6 * i + 0));
			bodies[i]->setLinearVelocity(sol.segment<3>(6 * i + 3));
			bodies[i]->updateTransformationMatrix(h);
		}
		C_.clear();
		colList.clear();
	}
}

void RigidBody::setEquality() {
	// Initialize G matrix
	int currentrow = 0;
	G_.clear();

	if (numJoints != 0) {
		setJointConstraints(currentrow, 1);
	}

	if (numFixed != 0) {
		setFixedConstraints(currentrow);
	}

	GG.resize(numEqualities, numVars);
	GG.setFromTriplets(G_.begin(), G_.end());

	// Initialize equality vector
	equalvec.resize(numEqualities);
	equalvec.setZero();

	program->setNumberOfEqualities(numEqualities);
	program->setEqualityMatrix(GG);
	program->setEqualityVector(equalvec);
}

void RigidBody::updateInertia() {

	if (isFEM == false) {
		// Update the inertia matrix of two rigid bodies connected by a spring 
		// Using gamma matrix (3x6) material Jacobian
		for (int i = 0; i < springs.size(); i++) {
			int ia = springs[i]->i;
			int ib = springs[i]->k;
			muscle_density = springs[i]->mass / (springs[i]->p1->x - springs[i]->p0->x).norm();

			auto b0 = bodies[ia];
			auto b1 = bodies[ib];

			auto n0 = b0->nodes[springs[i]->in];
			auto n1 = b1->nodes[springs[i]->kn];

			gamma.block<3, 3>(0, 0) = vec2crossmatrix(n0->x0).transpose();
			MatrixXd gamma_a = (b0->R * gamma);
			gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(n1->x0).transpose();
			MatrixXd gamma_b = (b1->R * gamma_k);

			Matrix12d IM;
			IM.setZero();

			// Add gravity force 
				Vector3d fg = springs[i]->mass * g/numFinitePoints;
				Matrix3x12d Jt;
				Jt.block<3, 6>(0, 0) = gamma_a;
				Jt.block<3, 6>(0, 6) = gamma_b;

				VectorXd wrench = Jt.transpose() * fg; // 12x1
				//bodies[ia]->setBodyForce(wrench.segment<6>(0));
				//bodies[ib]->setBodyForce(wrench.segment<6>(6));
				double dm = springs[i]->mass / numFinitePoints;

				Vector12d phi12;
				phi12.segment<6>(0) = b0->Phi;
				phi12.segment<6>(6) = b1->Phi;

				for (int t = 0; t < numFinitePoints; t++) {

					MatrixXd It = Jt.transpose() * Jt * dm; // It:12x12
					IM += It;
					bodies[ia]->setBodyForce(wrench.segment<6>(0));
					bodies[ib]->setBodyForce(wrench.segment<6>(6));
					//double ke = 0.5 * phi12.transpose() * It * phi12;
					//Vector3d pt = (len - t * dl) / len * springs[i]->p0->x + t * dl / len * springs[i]->p1->x;

					//double pe = dm * g.transpose() * pt;

					//spring_pe += pe;
					//spring_ke += ke;
				}


			/*MatrixXd I11 = gamma_a.transpose() * gamma_a * 1.0 / 3.0 * muscle_density;
			MatrixXd I12 = gamma_a.transpose() * gamma_b * 1.0 / 6.0 * muscle_density;
			MatrixXd I21 = gamma_b.transpose() * gamma_a * 1.0 / 6.0 * muscle_density;
			MatrixXd I22 = gamma_b.transpose() * gamma_b * 1.0 / 3.0 * muscle_density;*/

			MatrixXd I11 = IM.block<6, 6>(0, 0);
			MatrixXd I12 = IM.block<6, 6>(0, 6);
			MatrixXd I21 = IM.block<6, 6>(6, 0);
			MatrixXd I22 = IM.block<6, 6>(6, 6);

			for (int ii = 0; ii < 6; ii++) {
				for (int jj = 0; jj < 6; jj++) {
					A_.push_back(ETriplet(6 * ia + ii, 6 * ia + jj, I11(ii, jj)));
					A_.push_back(ETriplet(6 * ia + ii, 6 * ia + jj, I12(ii, jj)));
					A_.push_back(ETriplet(6 * ib + ii, 6 * ia + jj, I21(ii, jj)));
					A_.push_back(ETriplet(6 * ib + ii, 6 * ib + jj, I22(ii, jj)));
				}
			}
		}
	} else {
		// Integrate inertia matrix and compute material Jacobian using finite difference with mass samples
		// This method approximates the material Jacobian

		for (int i = 0; i < springs.size(); i++) {
			if (springs[i]->type == two_end_rbs) {
				int ia = springs[i]->i;
				int ib = springs[i]->k;
				auto b0 = bodies[ia];
				auto b1 = bodies[ib];
				double len = (springs[i]->p1->x - springs[i]->p0->x).norm();
				muscle_density = springs[i]->mass;

				// The approximated material Jacobian matrix of rb0, rb1
				Matrix3x6d J0, J1;
				J0.setZero();
				J1.setZero();

				Vector6d pert;
				// for each component of phi(i = 0, 1, 2..,11) add a relative small perturbation
				for (int k = 0; k < 6; k++) {

					pert.setZero();
					pert(k) = 1.0 * epsilon; // change kth component
					MatrixXd E_pert = b0->E * vec2crossmatrix(pert).exp();
					// Compute position of the end point p0 of spring
					Vector3d p_pert = transform(E_pert, springs[i]->p0->x0);
					J0.block(0, k, 3, 1) = 1.0 / epsilon * (p_pert - springs[i]->p0->x);
				}

				for (int k = 0; k < 6; k++) {

					pert.setZero();
					pert(k) = 1.0 * epsilon;

					MatrixXd E_pert = b1->E * vec2crossmatrix(pert).exp();

					Vector3d p_pert = transform(E_pert, springs[i]->p1->x0);
					J1.block(0, k, 3, 1) = 1.0 / epsilon * (p_pert - springs[i]->p1->x);
				}

				Matrix3x12d Jt;
				//Jt.block<3, 6>(0, 0) = J0;
				//Jt.block<3, 6>(0, 6) = J1;

			//	MatrixXd I11 = J0.transpose() * J0 * 1.0 / 3.0 * muscle_density;
			//	MatrixXd I12 = J0.transpose() * J1 * 1.0 / 6.0 * muscle_density;
			//	MatrixXd I21 = J1.transpose() * J0 * 1.0 / 6.0 * muscle_density;
			//	MatrixXd I22 = J1.transpose() * J1 * 1.0 / 3.0 * muscle_density;

			//	// Update the inertia matrix
			//for (int ii = 0; ii < 6; ii++) {
			//		for (int jj = 0; jj < 6; jj++) {
			//			A_.push_back(ETriplet(6 * ia + ii, 6 * ia + jj, I11(ii, jj)));
			//			A_.push_back(ETriplet(6 * ia + ii, 6 * ib + jj, I12(ii, jj)));
			//			A_.push_back(ETriplet(6 * ib + ii, 6 * ia + jj, I21(ii, jj)));
			//			A_.push_back(ETriplet(6 * ib + ii, 6 * ib + jj, I22(ii, jj)));
			//		}
			//	}
		
				// Add gravity force 
				Vector3d fg = springs[i]->mass * g;
				Vector12d wrench = Jt.transpose() * fg; 

				Jt.setZero();
				double dl = len / (numFinitePoints - 1);
				double ss;
				double dm = springs[i]->mass / numFinitePoints;
				Vector3d dfg = dm * g;
				Matrix12d IM; // add to body0 and body1
				IM.setZero();
				auto n0 = b0->nodes[springs[i]->in];
				auto n1 = b1->nodes[springs[i]->kn];

				for (int z = 0; z < numFinitePoints; z++) {
					ss = dl * z / len;

					gamma.block<3, 3>(0, 0) = vec2crossmatrix(n0->x0).transpose();
					MatrixXd gamma_a = (b0->R * gamma);
					gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(n1->x0).transpose();
					MatrixXd gamma_b = (b1->R * gamma_k);
					Jt.block<3, 6>(0, 0) = (1 - ss) * gamma_a;
					Jt.block<3, 6>(0, 6) = ss * gamma_b;

					MatrixXd It = Jt.transpose() * Jt * dm; // It:12x12
					IM += It;

					// compute the pos of z th point
					Vector3d zpoint = springs[i]->p0->x * (1 - ss) + springs[i]->p1->x * ss;
					Vector12d wr = Jt.transpose() * dfg;
					bodies[ia]->setBodyForce(wr.segment<6>(0));
					bodies[ib]->setBodyForce(wr.segment<6>(6));
				}
				
				Matrix6d I11 = IM.block<6, 6>(0, 0);
				Matrix6d I12 = IM.block<6, 6>(0, 6);
				Matrix6d I21 = IM.block<6, 6>(6, 0);
				Matrix6d I22 = IM.block<6, 6>(6, 6);

				for (int ii = 0; ii < 6; ii++) {
					for (int jj = 0; jj < 6; jj++) {
						A_.push_back(ETriplet(6 * ia + ii, 6 * ia + jj, I11(ii, jj)));
						A_.push_back(ETriplet(6 * ia + ii, 6 * ib + jj, I12(ii, jj)));
						A_.push_back(ETriplet(6 * ib + ii, 6 * ia + jj, I21(ii, jj)));
						A_.push_back(ETriplet(6 * ib + ii, 6 * ib + jj, I22(ii, jj)));
					}
				}
			}
		}

		for (int i = 0; i < numCylinders; i++) {
			auto c = cylinders[i];
			int ia = c->P->rb_id;
			int ib = c->S->rb_id;

			auto b0 = bodies[ia];
			auto b1 = bodies[ib];

			double total_length = c->l;
			muscle_density = c->mass / total_length;

			Matrix12d IM; // add to body0 and body1
			IM.setZero();

			if (wpc_stat[i] == wrap) {

				int numElements = numFinitePoints;
				double dm = c->mass / numElements;
				double dl = total_length / (numElements - 1);

				Vector3d length;
				length[0] = c->l_pc1;
				length[1] = c->l_c1c2;
				length[2] = c->l_c2s;

				Vector3i count;
				count[0] = floor(length[0] / dl);					// the number of elements in path P to C1
				count[2] = floor(length[2] / dl);					// the number of elements in path C2 to S
				count[1] = numElements - count[0] - count[2];		// the number of elements in path C1 to C2

				Vector6d pert;

				wpc_pert.setZero();

				// Jacobians for all the points
				MatrixXd J;
				J.resize(3 * numElements, 12);
				J.setZero();

				for (int k = 0; k < 12; k++) {

					pert.setZero();
					MatrixXd E_nopert, E_pert;
					// Compute the position after perturbation using wrapping
					Vector3d O_pert, P_pert, S_pert, Z_pert;

					if (k < 6) {
						// if move b0
						pert(k) = 1.0 * epsilon;	// change kth component
						E_nopert = b0->E;			// wiggle b0 when k < 6
						E_pert = E_nopert * vec2crossmatrix(pert).exp();
						P_pert = transform(E_pert, c->P->x0);
						S_pert = c->S->x;

						if (c->O->rb_id == c->P->rb_id) {
							// if cylinder is attached to b0, pos/vec changes
							O_pert = transform(E_pert, c->O->x0);
							Z_pert = transformVector(E_pert, transformVector(b0->E.inverse(), c->Z));
						}
						else {
							// if cylinder is attached to b1, pos/vec doesn't change
							O_pert = c->O->x;
							Z_pert = c->Z;
						}
					}
					else {
						// if move b1
						pert(k - 6) = 1.0 * epsilon;
						E_nopert = b1->E;
						E_pert = E_nopert * vec2crossmatrix(pert).exp();
						P_pert = c->P->x;
						S_pert = transform(E_pert, c->S->x0);

						if (c->O->rb_id == c->P->rb_id) {
							// if cylinder is attached to b0, pos/vec doesn't changes
							O_pert = c->O->x;
							Z_pert = c->Z;
						}
						else {
							// if cylinder is attached to b1, pos/vec changes
							O_pert = transform(E_pert, c->O->x0);
							Z_pert = transformVector(E_pert, transformVector(b1->E.inverse(), c->Z));
						}
					}

					// test if after perturbation I get the correct wrapping thing
					if (k == debug_i) {
						test_o = O_pert;
						test_z = Z_pert;
					}

					WrapCylinder wc(P_pert, S_pert, O_pert, Z_pert, c->r);
					wc.compute();

					Vector3d length_pert;
					double theta_s_pert, theta_e_pert;	// used to compute angle of the arc
					double total_length_pert;			// 
					Matrix3d _M_pert;					// some matrix used to compute the arc
					Vector3d C1_pert, C2_pert;

					if (wc.getStatus() == wrap) {
						wpc_stat_pert(k) = wrap;
						wpc_pert.block(k * 3, 0, 3, numWrapPoints + 1) = wc.getPoints(numWrapPoints, theta_s_pert, theta_e_pert, _M_pert);
						wpc_pert.block(k * 3, numWrapPoints + 1, 3, 1) = P_pert;
						wpc_pert.block(k * 3, numWrapPoints + 2, 3, 1) = S_pert;

						C2_pert = wpc_pert.block(k * 3, 0, 3, 1);

						for (int j = 0; j < numWrapPoints + 1; j++) {
							if (isnan(wpc_pert(k * 3, j))) {
								C1_pert = wpc_pert.block(k * 3, j - 1, 3, 1);
								break;
							}
							else {
								C1_pert = wpc_pert.block(k * 3, j, 3, 1);
							}
						}

						length_pert[0] = (C1_pert - P_pert).norm();
						length_pert[1] = wc.getLength();
						length_pert[2] = (C2_pert - S_pert).norm();

						total_length_pert = length_pert[0] + length_pert[1] + length_pert[2];
					}
					else {
						wpc_stat_pert(k) = no_wrap;
					}

					double dl_pert = total_length_pert / (numElements - 1);

					Vector3i count_pert;
					count_pert[0] = floor(length_pert[0] / dl_pert);				// the number of elements in path p to c1
					count_pert[2] = floor(length_pert[2] / dl_pert);				// the number of elements in path c2 to s
					count_pert[1] = numElements - count_pert[0] - count_pert[2];	// the number of elements in path C1 to C2

					for (int id = 0; id < numElements; id++) {
						Vector3d p_nopert;
						Vector3d p_pert;

						// Compute the position after perturbation
						if (id <= count_pert[0]) {
							double s = dl_pert * id / length_pert[0];
							p_pert = (1 - s) * P_pert + s * C1_pert;

						}
						else if (id > count_pert[0] && id < count_pert[0] + count_pert[1] - 1) {
							double dtheta = (dl_pert * id - length_pert[0]) / length_pert[1] * (theta_e_pert - theta_s_pert);
							double theta = theta_e_pert - dtheta;

							p_pert = _M_pert.transpose() * Vector3d(c->r * cos(theta), c->r * sin(theta), 0.0) + O_pert;
							p_pert(2) = c->P->x(2);
						}
						else {
							double s = (dl_pert * id - length_pert[1] - length_pert[0]) / length_pert[2];
							p_pert = (1 - s) * C2_pert + s * S_pert;
						}

						// Compute the position before perturbation
						// TODO: In fact, it should be moved outside the loop . Put it here for now
						if (id <= count[0]) {
							double s = dl * id / length[0];
							p_nopert = (1 - s) * c->P->x + s * c->c1;
						}
						else if (id > count[0] && id < count[0] + count[1] - 1) {
							double dtheta = (dl * id - length[0]) / length[1] * (c->theta_e - c->theta_s);
							double theta = c->theta_e - dtheta;

							p_nopert = c->M.transpose() * Vector3d(c->r * cos(theta), c->r * sin(theta), 0.0) + c->O->x;
							p_nopert(2) = c->P->x(2);

						}
						else {
							double s = (dl * id - length[0] - length[1]) / length[2];
							p_nopert = (1 - s) * c->c2 + s * c->S->x;
						}

						// to test if I get every element right
						if (k == debug_i) {
							test_pts.block<3, 1>(0, id) = p_nopert;
							test_pts_pert.block<3, 1>(0, id) = p_pert;
						}

						J.block(3 * id, k, 3, 1) = 1.0 / epsilon * (p_pert - p_nopert);
					}
				}

				Vector12d phi12;
				phi12.segment<6>(0) = b0->Phi;
				phi12.segment<6>(6) = b1->Phi;

				// Compute Inertia matrix by summing up all the JTJ*mu
				// Compute the total kinetic energy of a muscle
				for (int t = 0; t < numElements; t++) {
					
					Matrix3x12d Jt = J.block(3 * t, 0, 3, 12);
					Matrix12d It = Jt.transpose() * Jt * dm; 
					IM += It;

					// Add gravity force for each sample
					Vector3d fg =  dm * g;
					Vector12d wrench = Jt.transpose() * fg; 
					bodies[ia]->setBodyForce(wrench.segment<6>(0));
					bodies[ib]->setBodyForce(wrench.segment<6>(6));
				}

				Matrix6d I11 = IM.block<6, 6>(0, 0);
				Matrix6d I12 = IM.block<6, 6>(0, 6);
				Matrix6d I21 = IM.block<6, 6>(6, 0);
				Matrix6d I22 = IM.block<6, 6>(6, 6);

				for (int ii = 0; ii < 6; ii++) {
					for (int jj = 0; jj < 6; jj++) {
						A_.push_back(ETriplet(6 * ia + ii, 6 * ia + jj, I11(ii, jj)));
						A_.push_back(ETriplet(6 * ia + ii, 6 * ib + jj, I12(ii, jj)));
						A_.push_back(ETriplet(6 * ib + ii, 6 * ia + jj, I21(ii, jj)));
						A_.push_back(ETriplet(6 * ib + ii, 6 * ib + jj, I22(ii, jj)));
					}
				}
			}
			else {

				// if the status is not wrapped, update inertia the same way as springs
				gamma.block<3, 3>(0, 0) = vec2crossmatrix(c->P->x0).transpose();
				Matrix3d gamma_a = (b0->R * gamma);
				gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(c->S->x0).transpose();
				Matrix3d gamma_b = (b1->R * gamma_k);

				Matrix6d I11 = gamma_a.transpose() * gamma_a * 1.0 / 3.0 * muscle_density;
				Matrix6d I12 = gamma_a.transpose() * gamma_b * 1.0 / 6.0 * muscle_density;
				Matrix6d I21 = gamma_b.transpose() * gamma_a * 1.0 / 6.0 * muscle_density;
				Matrix6d I22 = gamma_b.transpose() * gamma_b * 1.0 / 3.0 * muscle_density;

				for (int ii = 0; ii < 6; ii++) {
					for (int jj = 0; jj < 6; jj++) {
						A_.push_back(ETriplet(6 * ia + ii, 6 * ia + jj, I11(ii, jj)));
						A_.push_back(ETriplet(6 * ia + ii, 6 * ib + jj, I12(ii, jj)));
						A_.push_back(ETriplet(6 * ib + ii, 6 * ia + jj, I21(ii, jj)));
						A_.push_back(ETriplet(6 * ib + ii, 6 * ib + jj, I22(ii, jj)));
					}
				}
			}
		}
	}
}

void RigidBody::setObjective(double h) {
	// Initialize A matrix
	A_.clear();

	for (int i = 0; i < numRB; i++) {
		bodies[i]->clearBodyForce(); // only gravity
		for (int j = 0; j < 6; j++) {
			A_.push_back(ETriplet(6 * i + j, 6 * i + j, bodies[i]->M(j, j)));
		}
	}
	
	Ao.resize(numVars, numVars);
	Ao.setFromTriplets(A_.begin(), A_.end());
	
	Aold = MatrixXd(Ao);

	updateInertia();

	A.resize(numVars, numVars);
	A.setFromTriplets(A_.begin(), A_.end());
	program->setObjectiveMatrix(A);
	
	// Initialize b vector
	if (numSprings != 0) {
		//computeSpringForces();
	}

	if (numCylinders != 0) {
		//computeWrapCylinderForces();
	}

	if (numDoubleCylinders != 0) {
		//computeWrapDoubleCylinderForces();
	}

	Adense = MatrixXd(A);
	MatrixXd PhiT;
	VectorXd twists, bodyforces;

	PhiT.resize(6 * numRB, 6 * numRB);
	twists.resize(6 * numRB);
	bodyforces.resize(6 * numRB);

	PhiT.setZero();
	bodyforces.setZero();
	twists.setZero();

	for (int i = 0; i < numRB; i++) {
		twists.segment<6>(i * 6) = bodies[i]->Phi;
		bodyforces.segment<6>(i * 6) = bodies[i]->B;
		PhiT.block(6 * i, 6 * i, 6, 6) = bodies[i]->PhiT;
	}

	RHS = Adense * twists + h * (PhiT * Aold * twists + bodyforces);
	program->setObjectiveVector(-RHS);
}

void RigidBody::postStabilization(int &currentrow) {
	equalvec.resize(numEqualities);
	equalvec.setZero();
	G_.clear();
	Matrix6d Gi, Gk;
	Matrix4d j0toB;
	j0toB.setZero();
	j0toB.block<3, 3>(0, 0) = I;
	j0toB.block<3, 1>(0, 3) << 0.0, dimensions(1)*0.5, 0.0;
	j0toB(3, 3) = 1;

	for (int i = 0; i < numJoints; i++) {
		Gi = joints[i]->computeAdjoint(joints[i]->Eij.inverse());

		joints[i]->computeEjkTemp(bodies);//btoJ

		Gk = -joints[i]->computeAdjoint(joints[i]->Ejk); // careful about the sign!

		if (joints[i]->type == BALL_JOINT) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 6; k++) {
					G_.push_back(ETriplet(currentrow, 6 * joints[i]->i + k, Gi(j + 3, k)));
					G_.push_back(ETriplet(currentrow, 6 * joints[i]->k + k, Gk(j + 3, k)));
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

		MatrixXd j0toJ0 = joints[i]->Ejk * j0toB;
		MatrixXd err = j0toJ0.log();

		if (isnan(err(0, 0))) {
			err.setZero();
		}

		VectorXd correctedg = crossmatrix2vec(err);
		equalvec.segment(i * 5, 2) = -correctedg.segment(0, 2);
		equalvec.segment(i * 5 + 2, 3) = -correctedg.segment(3, 3);
	}

	GG.resize(numEqualities, numVars);
	GG.setFromTriplets(G_.begin(), G_.end());

	program->setNumberOfEqualities(numEqualities);
	program->setEqualityMatrix(GG);
	program->setEqualityVector(equalvec);

	RHS.setZero();
	program->setObjectiveVector(RHS);
	bool success = program->solve();
	sol2 = program->getPrimalSolution();
}

void RigidBody::step(double h) {

	steps++;

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

	setObjective(h);
	setEquality();

	bool success = program->solve();
	sol = program->getPrimalSolution();

	for (int i = 0; i < numRB; i++) {	// Update node positions to detect collisions with floor
		bodies[i]->setAngularVelocity(sol.segment<3>(6 * i + 0));
		bodies[i]->setLinearVelocity(sol.segment<3>(6 * i + 3));
		bodies[i]->computeTempE(h);
	}
	int currentrow = 0;

	// To soft the joint drifting problem
	if (js["isPostStable"]) {
		postStabilization(currentrow);
	}
	
	numInequalities = 0;
	numColBoxBox = 0;
	numColFloor = 0;

	if (isFloorCol == true) {
		detectFloorCol();
	}

	if (isBoxBoxCol == true) {
		detectBoxBoxCol();
	}

	setInequality(h);
	updatePosNor();

	// use parameter baum = 0.1 to avoid some extreme situations, 1 is too large for this case
	if (js["isPostStable"]) {
		double baum = js["baum"];
		for (int i = 0; i < numRB; i++) {

			if (init_fixed_rb.size() >= 1 && i != init_fixed_rb[0]) {
				bodies[i]->correctPosition(-baum * sol2.segment<6>(6 * i + 0), baum);
			}
		}
		updatePosNor();
	}

	updateWrapCylinders();
	updateDoubleWrapCylinders();

	spring_pe = 0.0;
	spring_ke = 0.0;
	computeSpringEnergy();

	rb_ke = 0.0;
	rb_pe = 0.0;
	computeRigidBodyEnergy();

	cylinder_ke = 0.0;
	cylinder_pe = 0.0;
	computeCylinderEnergy();	
	
	saveData(js["save_step"]);
}

void RigidBody::saveData(int save_step) {
	// save data and plot comparison in MATLAB

	K.push_back(rb_ke + spring_ke + cylinder_ke);
	V.push_back(rb_pe + spring_pe + cylinder_pe);
	T.push_back(steps);

	if (steps == save_step) {
		cout << "finished" << endl;
		VectorXd Kvec;
		VectorXd Vvec;
		VectorXd Tvec;
		Kvec.resize(K.size());
		Vvec.resize(V.size());
		Tvec.resize(T.size());

		for (int i = 0; i < K.size(); i++) {
			Kvec(i) = K[i];
			Vvec(i) = V[i];
			Tvec(i) = T[i];
		}

		vec_to_file(Kvec, "K");
		vec_to_file(Vvec, "V");
		vec_to_file(Tvec, "T");
	}
}

void RigidBody::computeRigidBodyEnergy() {
	
	for (int i = 0; i < numRB; i++) {
		double pei = bodies[i]->mass * g.transpose() * bodies[i]->p;
		rb_pe += pei;

		double kei = 0.5 * bodies[i]->Phi.transpose() * bodies[i]->M * bodies[i]->Phi;
		rb_ke += kei;
	}
}

void RigidBody::computeSpringEnergy() {
	
	for (int i = 0; i < springs.size(); i++) {
		if (springs[i]->type == two_end_rbs) {
			int ia = springs[i]->i;
			int ib = springs[i]->k;
			auto b0 = bodies[ia];
			auto b1 = bodies[ib];
			double len = (springs[i]->p0->x - springs[i]->p1->x).norm();
			muscle_density = springs[i]->mass ;

			// The approximated material Jacobian matrix of rb0, rb1
			Matrix3x6d J0, J1;	
			J0.setZero();
			J1.setZero();

			Vector6d pert;
			// for each component of phi(i = 0, 1, 2..,11) add a relative small perturbation
			for (int k = 0; k < 6; k++) {

				pert.setZero();
				pert(k) = 1.0 * epsilon; // change kth component
				Matrix4d E_pert = b0->E * vec2crossmatrix(pert).exp();
				// Compute position of the end point p0 of spring
				Vector3d p_pert = transform(E_pert, springs[i]->p0->x0);
				J0.block(0, k, 3, 1) = 1.0 / epsilon * (p_pert - springs[i]->p0->x);
			}

			for (int k = 0; k < 6; k++) {
				pert.setZero();
				pert(k) = 1.0 * epsilon;
				Matrix4d E_pert = b1->E * vec2crossmatrix(pert).exp();
				Vector3d p_pert = transform(E_pert, springs[i]->p1->x0);
				J1.block(0, k, 3, 1) = 1.0 / epsilon * (p_pert - springs[i]->p1->x);
			}

			Matrix3x12d Jt;
			Jt.block<3, 6>(0, 0) = J0;
			Jt.block<3, 6>(0, 6) = J1;

			Matrix6d I11 = J0.transpose() * J0 * 1.0 / 3.0 * muscle_density;
			Matrix6d I12 = J0.transpose() * J1 * 1.0 / 6.0 * muscle_density;
			Matrix6d I21 = J1.transpose() * J0 * 1.0 / 6.0 * muscle_density;
			Matrix6d I22 = J1.transpose() * J1 * 1.0 / 3.0 * muscle_density;

			Matrix12d IM; // add to body0 and body1
			// Compute the total kinetic energy of a muscle
			IM.block<6, 6>(0, 0) = I11;
			IM.block<6, 6>(0, 6) = I12;
			IM.block<6, 6>(6, 0) = I21;
			IM.block<6, 6>(6, 6) = I22;

			Vector12d PHIM;
			PHIM.setZero();
			PHIM.segment<6>(0) = b0->Phi;
			PHIM.segment<6>(6) = b1->Phi;

			// Compute the total potential energy of a muscle
			double dl = len / (numFinitePoints - 1);
			double ss;
			double dm = springs[i]->mass / numFinitePoints;
			Vector3d dfg = dm * g;
			auto n0 = b0->nodes[springs[i]->in];
			auto n1 = b1->nodes[springs[i]->kn];
	
			IM.setZero();
			Jt.setZero();

			for (int z = 0; z < numFinitePoints; z++) {
				MatrixXd It = Jt.transpose() * Jt * dm; // It:12x12
				IM += It;
				// compute the pos of z th point
				ss = dl * z / len;
				Vector3d zpoint = springs[i]->p0->x * (1 - ss) + springs[i]->p1->x * ss;	
				double pe = dm * g.transpose() * zpoint;
				gamma.block<3, 3>(0, 0) = vec2crossmatrix(n0->x0).transpose();
				MatrixXd gamma_a = (b0->R * gamma);
				gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(n1->x0).transpose();
				MatrixXd gamma_b = (b1->R * gamma_k);
				Jt.block<3, 6>(0, 0) = (1 - ss)*gamma_a;
				Jt.block<3, 6>(0, 6) = ss * gamma_b;
				Vector3d zvel = Jt * PHIM;
				spring_ke += 0.5 * dm * zvel.transpose() * zvel;
				spring_pe += pe;
			}
		}
	}
}

void RigidBody::computeCylinderEnergy() {
	for (int i = 0; i < numCylinders; i++) {
		auto c = cylinders[i];
		int ia = c->P->rb_id;
		int ib = c->S->rb_id;

		auto b0 = bodies[ia];
		auto b1 = bodies[ib];

		double total_length = c->l;
		muscle_density = c->mass / total_length;

		Matrix12d IM; // add to body0 and body1
		IM.setZero();

		if (wpc_stat[i] == wrap) {

			int numElements = numFinitePoints;
			double dm = c->mass / numElements;
			double dl = total_length / (numElements - 1);

			Vector3d length;
			length[0] = c->l_pc1;
			length[1] = c->l_c1c2;
			length[2] = c->l_c2s;

			Vector3i count;
			count[0] = floor(length[0] / dl);					// the number of elements in path P to C1
			count[2] = floor(length[2] / dl);					// the number of elements in path C2 to S
			count[1] = numElements - count[0] - count[2];		// the number of elements in path C1 to C2

			Vector6d pert;
			wpc_pert.setZero();

			// Jacobians for all the points
			MatrixXd J;
			J.resize(3 * numElements, 12);
			J.setZero();

			for (int k = 0; k < 12; k++) {

				pert.setZero();
				Matrix4d E_nopert, E_pert;
				// Compute the position after perturbation using wrapping
				Vector3d O_pert, P_pert, S_pert, Z_pert;

				if (k < 6) {
					// if move b0
					pert(k) = 1.0 * epsilon;	// change kth component
					E_nopert = b0->E;			// wiggle b0 when k < 6
					E_pert = E_nopert * vec2crossmatrix(pert).exp();
					P_pert = transform(E_pert, c->P->x0);
					S_pert = c->S->x;

					if (c->O->rb_id == c->P->rb_id) {
						// if cylinder is attached to b0, pos/vec changes
						O_pert = transform(E_pert, c->O->x0);
						Z_pert = transformVector(E_pert, transformVector(b0->E.inverse(), c->Z));
					}
					else {
						// if cylinder is attached to b1, pos/vec doesn't change
						O_pert = c->O->x;
						Z_pert = c->Z;
					}
				}
				else {
					// if move b1
					pert(k - 6) = 1.0 * epsilon;
					E_nopert = b1->E;
					E_pert = E_nopert * vec2crossmatrix(pert).exp();
					P_pert = c->P->x;
					S_pert = transform(E_pert, c->S->x0);

					if (c->O->rb_id == c->P->rb_id) {
						// if cylinder is attached to b0, pos/vec doesn't changes
						O_pert = c->O->x;
						Z_pert = c->Z;
					}
					else {
						// if cylinder is attached to b1, pos/vec changes
						O_pert = transform(E_pert, c->O->x0);
						Z_pert = transformVector(E_pert, transformVector(b1->E.inverse(), c->Z));
					}
				}

				WrapCylinder wc(P_pert, S_pert, O_pert, Z_pert, c->r);
				wc.compute();

				Vector3d length_pert;
				double theta_s_pert, theta_e_pert;	// used to compute angle of the arc
				double total_length_pert;			// 
				Matrix3d _M_pert;					// some matrix used to compute the arc
				Vector3d C1_pert, C2_pert;

				if (wc.getStatus() == wrap) {
					wpc_stat_pert(k) = wrap;
					wpc_pert.block(k * 3, 0, 3, numWrapPoints + 1) = wc.getPoints(numWrapPoints, theta_s_pert, theta_e_pert, _M_pert);
					wpc_pert.block(k * 3, numWrapPoints + 1, 3, 1) = P_pert;
					wpc_pert.block(k * 3, numWrapPoints + 2, 3, 1) = S_pert;

					C2_pert = wpc_pert.block(k * 3, 0, 3, 1);

					for (int j = 0; j < numWrapPoints + 1; j++) {
						if (isnan(wpc_pert(k * 3, j))) {
							C1_pert = wpc_pert.block(k * 3, j - 1, 3, 1);
							break;
						}
						else {
							C1_pert = wpc_pert.block(k * 3, j, 3, 1);
						}
					}

					length_pert[0] = (C1_pert - P_pert).norm();
					length_pert[1] = wc.getLength();
					length_pert[2] = (C2_pert - S_pert).norm();

					total_length_pert = length_pert[0] + length_pert[1] + length_pert[2];
				}
				else {
					wpc_stat_pert(k) = no_wrap;
				}

				double dl_pert = total_length_pert / (numElements - 1);

				Vector3i count_pert;
				count_pert[0] = floor(length_pert[0] / dl_pert);				// the number of elements in path p to c1
				count_pert[2] = floor(length_pert[2] / dl_pert);				// the number of elements in path c2 to s
				count_pert[1] = numElements - count_pert[0] - count_pert[2];	// the number of elements in path C1 to C2

				for (int id = 0; id < numElements; id++) {
					Vector3d p_nopert;
					Vector3d p_pert;

					// Compute the position after perturbation
					if (id <= count_pert[0]) {
						double s = dl_pert * id / length_pert[0];
						p_pert = (1 - s) * P_pert + s * C1_pert;

					}
					else if (id > count_pert[0] && id < count_pert[0] + count_pert[1] - 1) {
						double dtheta = (dl_pert * id - length_pert[0]) / length_pert[1] * (theta_e_pert - theta_s_pert);
						double theta = theta_e_pert - dtheta;

						p_pert = _M_pert.transpose() * Vector3d(c->r * cos(theta), c->r * sin(theta), 0.0) + O_pert;
						p_pert(2) = c->P->x(2);
					}
					else {
						double s = (dl_pert * id - length_pert[1] - length_pert[0]) / length_pert[2];
						p_pert = (1 - s) * C2_pert + s * S_pert;
					}

					// Compute the position before perturbation
					// TODO: In fact, it should be moved outside the loop . Put it here for now
					if (id <= count[0]) {
						double s = dl * id / length[0];
						p_nopert = (1 - s) * c->P->x + s * c->c1;
					}
					else if (id > count[0] && id < count[0] + count[1] - 1) {
						double dtheta = (dl * id - length[0]) / length[1] * (c->theta_e - c->theta_s);
						double theta = c->theta_e - dtheta;

						p_nopert = c->M.transpose() * Vector3d(c->r * cos(theta), c->r * sin(theta), 0.0) + c->O->x;
						p_nopert(2) = c->P->x(2);

					}
					else {
						double s = (dl * id - length[0] - length[1]) / length[2];
						p_nopert = (1 - s) * c->c2 + s * c->S->x;
					}
					test_pts.block<3, 1>(0, id) = p_nopert;
					J.block(3 * id, k, 3, 1) = 1.0 / epsilon * (p_pert - p_nopert);
				}
			}

			Vector12d phi12;
			phi12.segment<6>(0) = b0->Phi;
			phi12.segment<6>(6) = b1->Phi;

			// Compute Inertia matrix by summing up all the JTJ*mu
			// Compute the total kinetic energy of a muscle

			for (int t = 0; t < numElements; t++) {
				Matrix3x12d Jt = J.block(3 * t, 0, 3, 12);
				Matrix12d It = Jt.transpose() * Jt * dm; 
				IM += It;

				double ke = 0.5 * phi12.transpose() * It * phi12;
				double pe = dm * g.transpose() * test_pts.block<3, 1>(0, t);
	
				cylinder_pe += pe;
				cylinder_ke += ke;
			}
		}
		else {
			// if the status is not wrapped, update inertia the same way as springs
			gamma.block<3, 3>(0, 0) = vec2crossmatrix(c->P->x0).transpose();
			Matrix3d gamma_a = (b0->R * gamma);
			gamma_k.block<3, 3>(0, 0) = vec2crossmatrix(c->S->x0).transpose();
			Matrix3d gamma_b = (b1->R * gamma_k);
			Matrix6d I11 = gamma_a.transpose() * gamma_a * 1.0 / 3.0 * muscle_density;
			Matrix6d I12 = gamma_a.transpose() * gamma_b * 1.0 / 6.0 * muscle_density;
			Matrix6d I21 = gamma_b.transpose() * gamma_a * 1.0 / 6.0 * muscle_density;
			Matrix6d I22 = gamma_b.transpose() * gamma_b * 1.0 / 3.0 * muscle_density;
		}
	}
}

void RigidBody::updateWrapCylinders() {

	for (int i = 0; i < numCylinders; i++) {
		auto c = cylinders[i];

		c->O->x = transform(bodies[c->O->rb_id]->E, c->O->x0);
		c->P->x = transform(bodies[c->P->rb_id]->E, c->P->x0);
		c->S->x = transform(bodies[c->S->rb_id]->E, c->S->x0);

		WrapCylinder wc(c->P->x, c->S->x, c->O->x, c->Z, c->r);
		wc.compute();

		if (wc.getStatus() == wrap) {
			double theta_s, theta_e;
			Matrix3d _M;

			wpc.block(int(3 * i), 0, 3, numWrapPoints + 1) = wc.getPoints(numWrapPoints, theta_s, theta_e, _M);

			wpc_stat(i) = wrap;
			wpc_length(i) = wc.getLength();

			Vector3d end = wpc.block<3, 1>(3 * i, 0);
			Vector3d start;
			for (int j = 0; j < numWrapPoints + 1; j++) {
				if (isnan(wpc(3 * i, j))) {
					start = wpc.block<3, 1>(3 * i, j - 1);
					break;
				}
				else {
					start = wpc.block<3, 1>(3 * i, j);
				}
			}
			c->c1 = start;
			c->c2 = end;

			c->theta_e = theta_e;
			c->theta_s = theta_s;
			c->M = _M;

			c->l_pc1 = (start - c->P->x).norm();
			c->l_c1c2 = wpc_length(i);
			c->l_c2s = (end - c->S->x).norm();

			c->l = c->l_c1c2 + c->l_c2s + c->l_pc1;


			c->pdir = (c->P->x - start).normalized();
			c->sdir = (c->S->x - end).normalized();

			if (cylinders[i]->L < 0.0) {
				cylinders[i]->L = cylinders[i]->l;
			}
		}
		else {
			wpc.block(3 * i, 0, 3, numWrapPoints + 1).setZero();
			wpc_stat(i) = no_wrap;
			wpc_length(i) = 0.0;

			cylinders[i]->l = wpc_length(i) + (cylinders[i]->S->x - cylinders[i]->P->x).norm();
			c->pdir = (c->P->x - c->S->x).normalized();
			c->sdir = -c->pdir;

			// also init L if it's the first time
			if (cylinders[i]->L < 0.0) {
				cylinders[i]->L = cylinders[i]->l;
			}
		}
	}
}

void RigidBody::updateDoubleWrapCylinders() {
	for (int i = 0; i < numDoubleCylinders; i++) {
		auto dc = doublecylinders[i];

		dc->U->x = transform(bodies[dc->U->rb_id]->E, dc->U->x0);
		dc->V->x = transform(bodies[dc->V->rb_id]->E, dc->V->x0);
		dc->P->x = transform(bodies[dc->P->rb_id]->E, dc->P->x0);
		dc->S->x = transform(bodies[dc->S->rb_id]->E, dc->S->x0);

		WrapDoubleCylinder wdc(dc->P->x,
			dc->S->x,
			dc->U->x,
			dc->Zu,
			dc->Ur,
			dc->V->x,
			dc->Zv,
			dc->Vr);

		wdc.compute();

		if (wdc.getStatus() == wrap) {


			wpdc.block(3 * i, 0, 3, 3 * numWrapPoints + 1) = wdc.getPoints(numWrapPoints);

			wpdc_stat(i) = 1;
			wpdc_length(i) = double(wdc.getLength());

			Vector3d start = wpdc.block<3, 1>(3 * i, 0);
			Vector3d end;

			for (int j = 0; j < 3 * numWrapPoints + 1; j++) {
				if (isnan(wpdc(3 * i, j))) {
					end = wpdc.block<3, 1>(3 * i, j - 1);
					break;
				}
				else {
					end = wpdc.block<3, 1>(3 * i, j);
				}
			}
			// set u1, u2, v1, v2
			dc->u1 = start;
			dc->u2 = wpdc.block<3, 1>(3 * i, numWrapPoints);
			dc->v2 = wpdc.block<3, 1>(3 * i, 2 * numWrapPoints);
			dc->v1 = end;

			dc->l = wpdc_length(i) + (end - dc->S->x).norm() + (start - dc->P->x).norm();

			dc->pdir = (dc->P->x - start).normalized();
			dc->sdir = (dc->S->x - end).normalized();
			dc->uvdir = (dc->u2 - dc->v1).normalized();

			if (dc->L < 0.0) {
				dc->L = dc->l;
			}
		}
		else {

			wpdc.block(3 * i, 0, 3, 3 * numWrapPoints + 1).setZero();
			wpdc_stat(i) = 0;
			wpdc_length(i) = 0.0;

			dc->l = wpdc_length(i) + (dc->S->x - dc->P->x).norm();

			dc->pdir = (dc->P->x - dc->S->x).normalized();
			dc->sdir = -dc->pdir;

			// also init L if it's the first time
			if (dc->L < 0.0) {
				dc->L = dc->l;
			}
		}
	}
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

void RigidBody::drawRBnodes()const {
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
			glColor3f(128.0 / 255.0, 0.0, 0.0);// maroon
			break;
		case 5:
			glColor3f(128.0 / 255.0, (128.0) / 255.0, 0.0);// olive
			break;
		case 6:
			glColor3f(0.0, 0.0, 128.0 / 255.0); // navy
			break;
		case 7:
			glColor3f(0.0, 128.0 / 255.0, 0.0); // green
			break;
		}
		for (int j = 0; j < numRB; j++) {
			glVertex3f(float(bodies[j]->nodes[i]->x(0)), float(bodies[j]->nodes[i]->x(1)), float(bodies[j]->nodes[i]->x(2)));
		}
	}
	glEnd();
}

void RigidBody::drawSprings()const {
	glColor3f(0.0, 0.0, 0.0); // black
	glLineWidth(3);
	for (int i = 0; i < springs.size(); i++) {
		Vector3d p0 = bodies[springs[i]->i]->nodes[springs[i]->in]->x;
		Vector3d p1 = bodies[springs[i]->k]->nodes[springs[i]->kn]->x;
		glBegin(GL_LINES);
		glVertex3f(float(p0(0)), float(p0(1)), float(p0(2)));
		glVertex3f(float(p1(0)), float(p1(1)), float(p1(2)));
		glEnd();
	}
}

void RigidBody::drawWrapCylinders()const {
	for (int t = 0; t < numCylinders; t++) {
		auto c = cylinders[t];
		// Draw Wrapper
		glColor3f(0.0, 0.0, 0.0); // black
		glLineWidth(3);
		glBegin(GL_LINE_STRIP);
		Vector3f end = c->S->x.cast <float>();
		Vector3f start = c->P->x.cast <float>();

		glVertex3f(end(0), end(1), end(2));
		if (wpc_stat(t) == wrap) {

			for (int i = 0; i < numWrapPoints + 1; i++) {
				Vector3f p = wpc.block<3, 1>(3 * t, i).cast<float>();
				glVertex3f(p(0), p(1), p(2));
			}
		}

		glVertex3f(start(0), start(1), start(2));
		glEnd();

		// Draw Wrapper Forces
		glColor3f(0.5, 0.5, 0.0);
		glLineWidth(5);
		glBegin(GL_LINES);

		Vector3d e = c->fp + c->P->x;
		Vector3d f = c->fs + c->S->x;
		Vector3d g = c->fpc + c->c1;
		Vector3d h = c->fsc + c->c2;

		glVertex3f(float(c->P->x(0)), float(c->P->x(1)), float(c->P->x(2)));
		glVertex3f(float(e(0)), float(e(1)), float(e(2)));

		glVertex3f(float(c->S->x(0)), float(c->S->x(1)), float(c->S->x(2)));
		glVertex3f(float(f(0)), float(f(1)), float(f(2)));

		glVertex3f(float(c->c1(0)), float(c->c1(1)), float(c->c1(2)));
		glVertex3f(float(g(0)), float(g(1)), float(g(2)));

		glVertex3f(float(c->c2(0)), float(c->c2(1)), float(c->c2(2)));
		glVertex3f(float(h(0)), float(h(1)), float(h(2)));

		glEnd();
	}
}

void RigidBody::drawWrapCylindersPerturbed()const {
	for (int t = 0; t < numCylinders; t++) {
		int k = 6;
		auto c = cylinders[t];
		// Draw Wrapper
		glColor3f(1.0, 0.0, 0.0); // black
		glLineWidth(2);
		glBegin(GL_LINE_STRIP);
		Vector3f end = wpc_pert.block(k * 3, numWrapPoints + 2, 3, 1).cast <float>();
		Vector3f start = wpc_pert.block(k * 3, numWrapPoints + 1, 3, 1).cast <float>();

		glVertex3f(end(0), end(1), end(2));
		if (wpc_stat_pert(k) == wrap) {

			for (int i = 0; i < numWrapPoints + 1; i++) {
				Vector3f p = wpc_pert.block<3, 1>(k * 3, i).cast<float>();
				glVertex3f(p(0), p(1), p(2));
			}
		}

		glVertex3f(start(0), start(1), start(2));
		glEnd();

		// Draw Cylinder after perturbation
		Vector3f op = test_o.cast<float>();
		Vector3f zaxis = (test_z + test_o).cast<float>();
		glBegin(GL_LINES);
		glVertex3f(op(0), op(1), op(2));
		glVertex3f(zaxis(0), zaxis(1), zaxis(2));
		glEnd();

		// Draw finite points before perturbation
		glColor3f(0.0, 1.0, 1.0);
		glPointSize(9.0f);
		glBegin(GL_POINTS);
		for (int i = 0; i < numFinitePoints; i++) {
			Vector3f pt = test_pts.block<3, 1>(0, i).cast<float>();
			glVertex3f(pt(0), pt(1), pt(2));
		}
		glEnd();

		// Draw finite points after perturbation

		glColor3f(0.0, 0.0, 1.0);
		glBegin(GL_POINTS);
		for (int i = 0; i < numFinitePoints; i++) {
			Vector3f ptp = test_pts_pert.block<3, 1>(0, i).cast<float>();
			glVertex3f(ptp(0), ptp(1), ptp(2));
		}
		glEnd();

	}
}

void RigidBody::drawDoubleWrapCylinders()const {


	for (int t = 0; t < numDoubleCylinders; t++) {
		auto dc = doublecylinders[t];
		Vector3f end = dc->S->x.cast <float>();
		Vector3f start = dc->P->x.cast <float>();

		// Draw Double Wrapper Points
		glColor3f(0.0, 0.0, 0.0);
		glPointSize(5);
		glBegin(GL_POINTS);
		glVertex3f(end(0), end(1), end(2));
		glVertex3f(start(0), start(1), start(2));
		glEnd();

		// Draw Double Wrapper
		glColor3f(0.0, 0.0, 0.0); // black
		glLineWidth(3);
		glBegin(GL_LINE_STRIP);

		glVertex3f(start(0), start(1), start(2));
		if (wpdc_stat(t) == 1) {

			for (int i = 0; i < numWrapPoints; i++) {
				if (isnan(wpdc(3 * t, i))) {
					break;
				}
				else {
					Vector3f p = wpdc.block<3, 1>(3 * t, i).cast<float>();
					glVertex3f(p(0), p(1), p(2));
				}
			}

			for (int i = 3 * numWrapPoints; i > 2 * numWrapPoints + 1; i--) {
				if (isnan(wpdc(3 * t, i))) {
					break;
				}
				else {
					Vector3f p = wpdc.block<3, 1>(3 * t, i).cast<float>();
					glVertex3f(p(0), p(1), p(2));
				}
			}

		}
		glVertex3f(end(0), end(1), end(2));

		glEnd();

		// Draw Wrapper Forces
		glColor3f(1.0, 0.0, 0.0); // red
		glLineWidth(5);
		glBegin(GL_LINES);
		Vector3d e = dc->fpu + dc->P->x;
		Vector3d f = dc->fsv + dc->S->x;
		Vector3d g = dc->fup + dc->u1;

		Vector3d h = dc->fuv + dc->u2;
		Vector3d q = dc->fvu + dc->v1;
		Vector3d p = dc->fvs + dc->v2;
		glVertex3f(float(dc->P->x(0)), float(dc->P->x(1)), float(dc->P->x(2)));
		glVertex3f(float(e(0)), float(e(1)), float(e(2)));

		glVertex3f(float(dc->S->x(0)), float(dc->S->x(1)), float(dc->S->x(2)));
		glVertex3f(float(f(0)), float(f(1)), float(f(2)));

		glVertex3f(float(dc->u1(0)), float(dc->u1(1)), float(dc->u1(2)));
		glVertex3f(float(g(0)), float(g(1)), float(g(2)));

		glVertex3f(float(dc->u2(0)), float(dc->u2(1)), float(dc->u2(2)));
		glVertex3f(float(h(0)), float(h(1)), float(h(2)));

		glVertex3f(float(dc->v1(0)), float(dc->v1(1)), float(dc->v1(2)));
		glVertex3f(float(q(0)), float(q(1)), float(q(2)));

		glVertex3f(float(dc->v2(0)), float(dc->v2(1)), float(dc->v2(2)));
		glVertex3f(float(p(0)), float(p(1)), float(p(2)));

		glEnd();
	}
}

void RigidBody::drawBoxBoxCol()const {
	glColor3f(1.0, 1.0, 0.0); // yellow
	if (numColBoxBox != 0) {
		glLineWidth(4);
		for (int i = 0; i < contacts.size(); i++) {
			for (int j = 0; j < contacts[i]->count; j++) {
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
		glUniform3fv(p->getUniform("kdBack"), 1, Vector3f((150 + 20 * i) / 255.0, (100 + 10 * i) / 255.0, (150.0 + 10 * i) / 255.0).data());
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

	drawRBnodes();
	drawSprings();
	drawBoxBoxCol();
	drawWrapCylinders();
	drawWrapCylindersPerturbed();
	drawDoubleWrapCylinders();

	MV->popMatrix();
	p2->unbind();
}

shared_ptr<Spring> RigidBody::createSpring2RB(int _i, int _k, int _in, int _kn, vector < shared_ptr<RBState> > bodies, double E, double _muscle_mass)
{
	auto s = make_shared<Spring>(bodies[_i]->nodes[_in], bodies[_k]->nodes[_kn]);
	s->i = _i;
	s->k = _k;
	s->in = _in;
	s->kn = _kn;
	s->E = E;
	s->type = two_end_rbs;
	s->mass = _muscle_mass;
	Vector3d x0 = transform(bodies[_i]->E, bodies[_i]->nodes[_in]->x0);
	Vector3d x1 = transform(bodies[_k]->E, bodies[_k]->nodes[_kn]->x0);
	Vector3d dx = x1 - x0;
	s->L = dx.norm();
	return s;
}

shared_ptr<Spring> RigidBody::createSpring1RB(int _i, int _in, Vector3d pos, vector <shared_ptr<RBState> > bodies, double E, double _muscle_mass) {
	auto p0 = make_shared<Particle>();
	p0->x0 = pos;
	p0->x = pos;

	auto s = make_shared<Spring>(p0, bodies[_i]->nodes[_in]);
	s->i = -1; // in world frame
	s->k = _i;
	s->in = -1; // no index
	s->kn = _in;
	s->E = E;
	s->type = one_end_fixed;
	s->mass = _muscle_mass;
	Vector3d x1 = transform(bodies[_i]->E, bodies[_i]->nodes[_in]->x0);
	Vector3d dx = x1 - pos;
	s->L = dx.norm();
	return s;
}