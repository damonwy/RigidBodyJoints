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
struct RBState;
struct Joint;

typedef Eigen::Triplet<double> ETriplet;

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
	Eigen::VectorXi init_fixed_rb;

	double mass;
	int numFixed;
	int numCol;
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

	
	Eigen::VectorXd xl;
	Eigen::VectorXd xu;
	Eigen::VectorXd b;
	Eigen::Matrix3d I;
	Eigen::Vector3d g; // gravity
	Eigen::VectorXd RHS;
	Eigen::VectorXd sol;
	Eigen::VectorXd equalvec;
	
	int numRB;
	Eigen::MatrixXd Eij;
	std::vector < std::shared_ptr<Joint> > joints;
	std::vector <std::shared_ptr<RBState> > bodies;

	static const int BALL_JOINT = 3;
	static const int HINGE_JOINT = 5;

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