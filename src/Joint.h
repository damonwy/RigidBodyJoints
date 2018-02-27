#pragma once

#ifndef __Joint__
#define __Joint__

#include <vector>
#include <memory>
#include <tetgen.h>
#include "QuadProgMosek.h"

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class MatrixStack;
class Program;
class Particle;
struct RBState;

typedef Eigen::Triplet<double> ETriplet;

struct Joint {
	static const int BALL_JOINT = 3; // number is the dof
	static const int HINGE_JOINT = 5;

	// fixed
	Eigen::MatrixXd Eij;   // j w.r.t i
	Eigen::MatrixXd Gi;	   // Used to constraint Phi_i


	int i; // index of the first rigid body in RigidBody::bodies
	int k; // index of the second rigid body
	int index; // index of joint
	int type; // the type of joint

			  // computed
	Eigen::MatrixXd Ejk;   // k w.r.t j
	Eigen::MatrixXd Gk;	   // Used to constraint Phi_k

	Joint() {
		Eij.resize(4, 4);
		Eij.setZero();
		Eij(3, 3) = 1;
		Gi.resize(6, 6);
		Gi.setZero();
		Gk.resize(6, 6);
		Gk.setZero();
	};

	Eigen::MatrixXd computeAdjoint(Eigen::MatrixXd E);
	Eigen::Matrix3d vec2crossmatrix(Eigen::Vector3d a);
	void computeEjk(const std::vector< std::shared_ptr<RBState> > bodies);
};

#endif