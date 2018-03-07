#pragma once
#ifndef RigidBodyJoints_SRC_RIGIDBODY_H_
#define RigidBodyJoints_SRC_RIGIDBODY_H_

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
struct Joint;
struct Contacts;
class Spring;

typedef Eigen::Triplet<double> ETriplet;

class RigidBody {
public:
	RigidBody();
	void init();
	void updatePosNor();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p, const std::shared_ptr<Program> p2, std::shared_ptr<MatrixStack> P)const;
	void step(double h);
	std::shared_ptr<Spring> createSpring(int _i, int _k, int _in, int _kn, std::vector < std::shared_ptr<RBState> > bodies, double E);
	Eigen::MatrixXd computeAdjoint(Eigen::MatrixXd E);
	Eigen::Matrix3d vec2crossmatrix(Eigen::Vector3d a);	// repackage a vector into a cross-product matrix
	Eigen::Vector3d local2world(Eigen::MatrixXd E, Eigen::Vector3d x);	// compute the world position given a local position x on a rigid body
	Eigen::Vector3d world2local(Eigen::MatrixXd E, Eigen::Vector3d x);
	tetgenio in, out;

	int nVerts;
	int nTriFaces;
	int nEdges;
	Eigen::VectorXi init_fixed_rb;

	double mass;
	double yfloor;
	int numFixed;
	int numColFloor;
	int numColBoxBox;
	int numJoints;
	int numVars;
	int numEqualities;
	int numInequalities;

	Eigen::Vector3d ynormal;
	std::shared_ptr<QuadProg> program;
	std::vector<ETriplet> A_;
	Eigen::SparseMatrix<double> A;
	std::vector<ETriplet> G_;
	Eigen::SparseMatrix<double> GG;
	std::vector<ETriplet> C_;
	Eigen::SparseMatrix<double> C;


	Eigen::VectorXd xl;
	Eigen::VectorXd xu;
	Eigen::VectorXd b;
	Eigen::Matrix3d I;
	Eigen::Vector3d g; // gravity
	Eigen::VectorXd RHS;
	Eigen::VectorXd sol;
	Eigen::VectorXd equalvec;
	Eigen::VectorXd convec;
	Eigen::VectorXd conveck;
	Eigen::MatrixXd gamma;
	Eigen::MatrixXd gamma_k;
	Eigen::VectorXd joint_forces;

	int numRB;
	Eigen::MatrixXd Eij;
	std::vector < std::shared_ptr<Joint> > joints;
	std::vector < std::shared_ptr<RBState> > bodies;
	std::vector < std::shared_ptr<Contacts> > contacts;
	std::vector < std::shared_ptr<Spring> > springs;

	std::vector < int > colList;

	static const int BALL_JOINT = 3;
	static const int HINGE_JOINT = 5;
	static const int HINGE_JOINT_X = 0;
	static const int HINGE_JOINT_Y = 1;
	static const int HINGE_JOINT_Z = 2;

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

#endif  // RigidBodyJoints_SRC_RIGIDBODY_H_