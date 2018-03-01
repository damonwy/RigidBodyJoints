#pragma once

#ifndef __RBState__
#define __RBState__

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

typedef Eigen::Triplet<double> ETriplet;

struct RBState
{
	std::vector <std::shared_ptr<Particle>>nodes;
	// constants
	double mass;
	Eigen::MatrixXd M;

	// states
	Eigen::VectorXd Phi;	// 6x1 angular and linear velocity repackeaged into 6-vector
	Eigen::MatrixXd PHI;	// 6x6
	Eigen::MatrixXd PhiT;	// 6x6 spatial cross product matrix cosisting of [omega_i]

	Eigen::Vector3d Omega;	// 3x1 angular velocity
	Eigen::Matrix3d OMEGA;	// 3X3 the crossproduct matrix of Omega
	Eigen::Vector3d V;		// 3x1 linear velocity
	Eigen::Matrix3d VC;		// 3x3 the crossproduct matrix of V

	Eigen::MatrixXd E;      // 4x4 transformation matrix consisting of rotational and translational components
	Eigen::MatrixXd Etemp;  // 4x4 transformation matrix (temporary)
	Eigen::MatrixXd EE;		// 4x4 spatial velocity matrix	

	Eigen::Vector3d p;		// 3x1 position of the local frame's origin in world coordinates
	Eigen::Matrix3d R;		// 3x3 rotation matrix of which each colum corresponds to the frame's basis vectors e_k, expressed in world coordinates

							// computed
	Eigen::VectorXd B; // 6x1 forces
	Eigen::Vector3d dimensions;  // width, height, length

	// set functions
	void setMass(double _mass);
	void setDimensions(Eigen::Vector3d _dimensions);
	void setLinearVelocity(Eigen::Vector3d _V);
	void setAngularVelocity(Eigen::Vector3d _Omega);
	void setSpatialInertiaMatrix();
	void setTransformationMatrix(Eigen::MatrixXd _E);
	void setRotational(Eigen::Matrix3d _R);
	void setLocalFrameOrigin(Eigen::Vector3d _p);
	void setBodyForce();

	Eigen::VectorXd computeForces(double h);
	void computeTempE(double h);

	// update
	void updatePosition();
	void updateTransformationMatrix(double h);

	Eigen::Vector3d local2world(Eigen::MatrixXd E, Eigen::Vector3d x);	// compute the world position given a local position x on a rigid body
	Eigen::Matrix3d vec2crossmatrix(Eigen::Vector3d a);					// repackage a vector into a cross-product matrix

	RBState() {
		E.resize(4, 4);
		EE.resize(4, 4);
		Etemp.resize(4, 4);
		E.setZero();
		E(3, 3) = 1.0;
		EE.setZero();
		B.resize(6);
		B.setZero();
		PHI.resize(6, 6);
		PHI.setZero();
		Phi.resize(6);
		PhiT.resize(6, 6);
		M.resize(6, 6);
	}
};

#endif