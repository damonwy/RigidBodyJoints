#pragma once
#ifndef RigidBodyJoints_SRC_JOINT_H_
#define RigidBodyJoints_SRC_JOINT_H_

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
	static const int HINGE_JOINT_X = 0;
	static const int HINGE_JOINT_Y = 1;
	static const int HINGE_JOINT_Z = 2;

	// fixed
	Eigen::MatrixXd Eij;   // j w.r.t i
	Eigen::MatrixXd Gi;	   // Used to constraint Phi_i


	int i;		// index of the first rigid body in RigidBody::bodies
	int k;		// index of the second rigid body
	int index;	// index of joint
	int type;	// the type of joint
	int hinge_type;

	// computed
	Eigen::MatrixXd Ejk;   // k w.r.t j
	Eigen::MatrixXd Gk;	   // Used to constraint Phi_k


	Eigen::MatrixXd computeAdjoint(Eigen::Matrix4d E);
	Eigen::Matrix3d vec2crossmatrix(Eigen::Vector3d a);
	void computeEjk(const std::vector< std::shared_ptr<RBState> > bodies);
	void computeEjkTemp(const std::vector< std::shared_ptr<RBState> > bodies);
	Joint();
};

#endif // RigidBodyJoints_SRC_JOINT_H_