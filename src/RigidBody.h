#pragma once
#ifndef __RigidBody__
#define __RigidBody__

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
	Eigen::MatrixXd M; // fixed

	// states
	
	Eigen::VectorXd Phi;	// 6x1 angular and linear velocity repackeaged into 6-vector
	Eigen::MatrixXd PHI;	// 6x6
	Eigen::MatrixXd PhiT;	// 6x6 spatial cross product matrix cosisting of [omega_i]
	
	Eigen::Vector3d Omega;	// 3x1 angular velocity
	Eigen::Matrix3d OMEGA;	// 3X3 the crossproduct matrix of Omega
	Eigen::Vector3d V;		// 3x1 linear velocity
	Eigen::Matrix3d VC;		// 3x3 the crossproduct matrix of V

	Eigen::MatrixXd E;      // 4x4 transformation matrix consisting of rotational and translational components
	Eigen::MatrixXd EE;		// 4x4 spatial velocity matrix	
	
	Eigen::Vector3d p;		// 3x1 position of the local frame's origin in world coordinates
	Eigen::Matrix3d R;		// 3x3 rotation matrix of which each colum corresponds to the frame's basis vectors e_k, expressed in world coordinates
	
	// computed
	Eigen::VectorXd B; // 6x1 forces

	// get functions
	double getMass()
	{
		return mass;
	}
	Eigen::Vector3d getLinearVelocity() {
		return V;
	}
	Eigen::Vector3d getAngularVelocity() {
		return Omega;
	}
	Eigen::Matrix3d getRotational() {
		return R;
	}
	Eigen::Vector3d getLocalFrameOrigin() {
		return p;
	}


	// set functions
	void setMass(double _mass) {
		mass = _mass;
	}
	void setLinearVelocity(Eigen::Vector3d _V) {
		this->V = _V;
		Phi.segment<3>(3) = V;
		EE.block<3,1>(0, 3) = V;
	}
	void setAngularVelocity(Eigen::Vector3d _Omega) {
		this->Omega = _Omega;
		Phi.segment<3>(0) = Omega;
		OMEGA = vec2crossmatrix(Omega);
		EE.block<3,3>(0, 0) = OMEGA;
		PHI.block(0, 0, 3, 3) = OMEGA;
		PHI.block(3, 3, 3, 3) = OMEGA;
		PhiT = PHI.transpose();

	}
	void setSpatialInertiaMatrix() {
		M.resize(6, 6);
		this->M << mass / 3.0 * 10, 0, 0, 0, 0, 0,
					0, mass / 3.0 * 2.0, 0, 0, 0, 0,
					0, 0, mass / 3.0 * 10.0, 0, 0, 0,
					0, 0, 0, mass, 0, 0,
					0, 0, 0, 0, mass, 0,
					0, 0, 0, 0, 0, mass;

	}
	void setTransformationMatrix(Eigen::MatrixXd _E) {
		this->E = _E;
		R = E.block(0, 0, 3, 3);
		p = E.block(0, 3, 3, 1);
		setBodyForce();
	}
	void setRotational(Eigen::Matrix3d _R) {
		this->R = _R;
		this->E.block<3, 3>(0, 0) = R; // update E
		setBodyForce();
	}
	void setLocalFrameOrigin(Eigen::Vector3d _p) {
		this->p = _p;
		E.block<3, 1>(0, 3) = p;
	}
	void setBodyForce() {
		// if only gravity is involved
		Eigen::Vector3d g;
		g << 0.0, -9.8, 0.0;
		B.segment<3>(3) = R.transpose() * mass * g;
	}

	Eigen::VectorXd computeForces(double h);

	// update
	void updatePosition();
	void updateTransformationMatrix(double h);

	Eigen::Vector3d local2world(Eigen::MatrixXd E, Eigen::Vector3d x);	// compute the world position given a local position x on a rigid body
	Eigen::Matrix3d vec2crossmatrix(Eigen::Vector3d a);		// repackage a vector into a cross-product matrix

	RBState() {
		E.resize(4, 4);
		E.setZero();
		E(3, 3) = 1.0;
		EE.resize(4, 4);
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


class RigidBody {
public:
	RigidBody();
	void init();
	void updatePosNor();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p)const;
	void step(double h);

	Eigen::MatrixXd computeAdjoint(Eigen::MatrixXd E);
	Eigen::Matrix3d vec2crossmatrix(Eigen::Vector3d a);	// repackage a vector into a cross-product matrix

	tetgenio in, out;
	
	int nVerts;
	int nTriFaces;
	int nEdges;

	double mass;
	int numCol;
	int numJoints;
	int numVars;
	
	Eigen::Vector3d ynormal;
	std::shared_ptr<QuadProg> program;
	std::vector<ETriplet> A_;
	Eigen::SparseMatrix<double> A;
	
	Eigen::VectorXd xl;
	Eigen::VectorXd xu;
	Eigen::VectorXd b;
	Eigen::Matrix3d I;
	Eigen::Vector3d g; // gravity
	Eigen::VectorXd RHS;
	Eigen::VectorXd sol;
	
	int numRB;
	Eigen::MatrixXd Eij;
	std::vector <std::shared_ptr<RBState> > bodies;

private:
	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;
	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;
};

#endif