#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "Helper.h"

using namespace std;
using namespace Eigen;


MatrixXd Helper::vec2crossmatrix(VectorXd a) {
	MatrixXd A;
	if (a.size() == 3) {
		// dim = 3
		A.resize(3, 3);
		A.setZero();
		A << 0, -a(2), a(1),
			a(2), 0, -a(0),
			-a(1), a(0), 0;
	}
	else {
		// dim = 6
		A.resize(4, 4);
		A.setZero();
		A.block<3, 1>(0, 3) = a.segment<3>(3);
		A.block<3, 3>(0, 0) = vec2crossmatrix(a.segment<3>(0));
	}
	return A;
}

VectorXd Helper::crossmatrix2vec(MatrixXd A) {
	VectorXd x;
	if (A.cols() < 4) {
		x.resize(3);
		x << A(2, 1), A(0, 2), A(1, 0);
	}
	else {
		x.resize(6);
		x << A(2, 1), A(0, 2), A(1, 0), A(0, 3), A(1, 3), A(2, 3);
	}
	return x;
}

Vector3d Helper::transform(MatrixXd E, Vector3d x) {
	// homo coor
	VectorXd xh;
	xh.resize(4);
	xh.segment<3>(0) = x;
	xh(3) = 1.0;

	Vector3d xw = (E * xh).segment<3>(0);
	return xw;
}

Vector3d Helper::transformVector(MatrixXd E, Vector3d vec) {
	VectorXd xh;
	xh.resize(4);
	xh.segment<3>(0) = vec;
	xh(3) = 0.0;
	Vector3d xw = (E * xh).segment<3>(0);
	return xw;
}