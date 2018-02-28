#define TETLIBRARY
#include "Joint.h"
#include <iostream>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Particle.h"
#include "MatrixStack.h"
#include "Program.h"
#include "GLSL.h"
#include "RBState.h"

#include <unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h>

using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> ETriplet;


Joint::Joint() {
	Eij.resize(4, 4);
	Ejk.resize(4, 4);
	Eij.setZero();
	Eij(3, 3) = 1;
	Gi.resize(6, 6);
	Gi.setZero();
	Gk.resize(6, 6);
	Gk.setZero();
};

Matrix3d Joint::vec2crossmatrix(Vector3d a) {
	Matrix3d A;
	A << 0, -a(2), a(1),
		a(2), 0, -a(0),
		-a(1), a(0), 0;
	return A;
}

void Joint::computeEjk(const vector< shared_ptr<RBState> > bodies) {
	Ejk = Eij.inverse() * bodies[i]->E.inverse() * bodies[k]->E;
}

MatrixXd Joint::computeAdjoint(Matrix4d E) {
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