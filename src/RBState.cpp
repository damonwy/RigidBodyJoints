#define TETLIBRARY
#include "RBState.h"
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

void RBState::setDimensions(Eigen::Vector3d _dimensions) {
	this->dimensions = _dimensions;
}

void RBState::setTransformationMatrix(MatrixXd _E) {
	this->E = _E;
	R = E.block(0, 0, 3, 3);
	p = E.block(0, 3, 3, 1);
	VectorXd force;
	force.resize(6);
	force.setZero();

	setBodyForce(force);
}

void RBState::setRotational(Matrix3d _R) {
	this->R = _R;
	this->E.block<3, 3>(0, 0) = R; // update E
	VectorXd force;
	force.resize(6);
	force.setZero();
	setBodyForce(force);
}

void RBState::setSpatialInertiaMatrix() {
	M.resize(6, 6);
	this->M << mass / 12.0 * (dimensions(1)*dimensions(1) + dimensions(2)+dimensions(2)), 0, 0, 0, 0, 0,
			  0, mass / 12.0 * (dimensions(0)*dimensions(0) + dimensions(2) + dimensions(2)), 0, 0, 0, 0,
			  0, 0, mass / 12.0 * (dimensions(1)*dimensions(1) + dimensions(0) + dimensions(0)), 0, 0, 0,
		0, 0, 0, mass, 0, 0,
		0, 0, 0, 0, mass, 0,
		0, 0, 0, 0, 0, mass;
}

void RBState::setLocalFrameOrigin(Vector3d _p) {
	this->p = _p;
	E.block<3, 1>(0, 3) = p;
}

void RBState::setMass(double _mass) {
	mass = _mass;
}

void RBState::setLinearVelocity(Vector3d _V) {
	this->V = _V;
	Phi.segment<3>(3) = V;
	EE.block<3, 1>(0, 3) = V;
}

void RBState::setAngularVelocity(Eigen::Vector3d _Omega) {
	this->Omega = _Omega;
	Phi.segment<3>(0) = Omega;
	OMEGA = vec2crossmatrix(Omega);
	EE.block<3, 3>(0, 0) = OMEGA;
	PHI.block(0, 0, 3, 3) = OMEGA;
	PHI.block(3, 3, 3, 3) = OMEGA;
	PhiT = PHI.transpose();
}

void RBState::setBodyForce(VectorXd force) {
	// if only gravity is involved
	Eigen::Vector3d g;
	g << 0.0, -9.8, 0.0;
	B.setZero();
	B.segment<3>(3) = R.transpose() * mass * g;
	B = B + force;
}


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

void RBState::computeTempE(double h) {
	Etemp = E * (h * EE).exp();
}

void RBState::updateTransformationMatrix(double h) {
	E = E * (h * EE).exp();
	R = E.block(0, 0, 3, 3);
	p = E.block(0, 3, 3, 1);
	VectorXd force;
	force.resize(6);
	force.setZero();

	setBodyForce(force);
}

