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
#include <Eigen/SparseLU>


class MatrixStack;
class Program;
class Particle;
struct RBState;
struct Joint;
struct Contacts;
class Spring;
class Cylinder;
class DoubleCylinder;
class WrapCylinder;
class WrapDoubleCylinder;

typedef Eigen::Triplet<double> ETriplet;

class RigidBody {
public:
	RigidBody();
	void init();
	void initAfterNumRB();
	void initAfterNumCylinders();
	void initAfterNumDoubleCylinders();
	void initConstant();
	void initRBShape();
	void initJoints();
	void initSprings(double stiffness);
	void initCylinder();
	void initDoubleCylinder();

	void initRBs();
	void initBuffers();
	void updatePosNor();
	void updateWrapCylinders();
	void updateDoubleWrapCylinders();
	void postStabilization(int &currentrow);
	void updateInertia(double h);

	void computeSpringForces();
	void computeWrapCylinderForces();
	void computeWrapDoubleCylinderForces();
	void setJointConstraints(int &currentrow, int type);
	void setFixedConstraints(int &currentrow);

	void setEquality();
	void setInequality(double h);
	void setObjective(double h);

	void detectFloorCol();
	void detectBoxBoxCol();

	void drawRBnodes()const;
	void drawSprings()const;
	void drawWrapCylinders()const;
	void drawWrapCylindersPerturbed()const;
	void drawDoubleWrapCylinders()const;
	void drawBoxBoxCol()const;
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p, const std::shared_ptr<Program> p2, std::shared_ptr<MatrixStack> P)const;
	

	void step(double h);
	std::shared_ptr<Spring> createSpring2RB(int _i, int _k, int _in, int _kn, std::vector < std::shared_ptr<RBState> > bodies, double E, double _muscle_mass);
	std::shared_ptr<Spring> createSpring1RB(int _i, int _in, Eigen::Vector3d pos, std::vector <std::shared_ptr<RBState> > bodies, double E, double _muscle_mass);
	Eigen::MatrixXd computeAdjoint(Eigen::MatrixXd E);
	Eigen::MatrixXd vec2crossmatrix(Eigen::VectorXd a);					// repackage a vector into a cross-product matrix
	Eigen::VectorXd crossmatrix2vec(Eigen::MatrixXd A);
	Eigen::Vector3d transform(Eigen::MatrixXd E, Eigen::Vector3d x);	// compute the world position given a local position x on a rigid body
	Eigen::Vector3d transformVector(Eigen::MatrixXd E, Eigen::Vector3d vec);
	tetgenio in, out;

	double stiffness;
	double epsilon;
	double mass;
	double muscle_mass;
	double yfloor;
	double muscle_density;

	int nVerts;
	int nTriFaces;
	int nEdges;

	int numRB;
	int numFixed;
	int numColFloor;
	int numColBoxBox;
	int numJoints;
	int numSprings;
	int numCylinders;
	int numDoubleCylinders;
	int numWrapPoints;   
	int numFinitePoints;

	int numVars;
	int numEqualities;
	int numInequalities;

	int debug_i;

	bool isBoxBoxCol;	// is box to box collision on?
	bool isFloorCol;	// is box to floor collision on?
	bool isFEM;			// use FEM to compute Jacobian Matrix?

	Eigen::Vector3d dimensions;

	Eigen::Vector3d ynormal;
	std::shared_ptr<QuadProg> program;
	std::vector<ETriplet> A_;
	Eigen::SparseMatrix<double> A;
	std::vector<ETriplet> G_;
	Eigen::SparseMatrix<double> GG;
	std::vector<ETriplet> C_;
	Eigen::SparseMatrix<double> C;

	Eigen::VectorXd init_v; // used to init the linear velocity of all the rigid bodies
	Eigen::VectorXd init_w; // ................ angular velocity ......................
	Eigen::VectorXd init_p; // ................ position ..............................
	Eigen::MatrixXd init_R; // ................ rotation ..............................
	Eigen::VectorXi init_fixed_rb;
	Eigen::VectorXd init_cyl_P; // in local frame
	Eigen::VectorXd init_cyl_S; // ...
	Eigen::VectorXd init_cyl_O; // ...
	Eigen::VectorXd init_cyl_Z; 
	Eigen::VectorXd init_cyl_r; 
	Eigen::VectorXi init_cyl_O_rb;
	Eigen::VectorXi init_cyl_P_rb;
	Eigen::VectorXi init_cyl_S_rb;

	Eigen::VectorXd init_dcyl_P; // in local frame
	Eigen::VectorXd init_dcyl_S; // ...
	Eigen::VectorXi init_dcyl_P_rb;
	Eigen::VectorXi init_dcyl_S_rb;

	Eigen::VectorXd init_dcyl_U; 
	Eigen::VectorXd init_dcyl_Uz;
	Eigen::VectorXd init_dcyl_Ur;
	Eigen::VectorXi init_dcyl_U_rb;

	Eigen::VectorXd init_dcyl_Vz;
	Eigen::VectorXd init_dcyl_Vr;
	Eigen::VectorXd init_dcyl_V;
	Eigen::VectorXi init_dcyl_V_rb;

	Eigen::VectorXd xl;
	Eigen::VectorXd xu;
	Eigen::VectorXd b;
	Eigen::Matrix3d I;
	Eigen::Vector3d g; // gravity
	Eigen::VectorXd RHS;
	Eigen::VectorXd sol;
	Eigen::VectorXd sol2;
	Eigen::VectorXd equalvec;
	Eigen::VectorXd convec;
	Eigen::VectorXd conveck;
	Eigen::MatrixXd gamma;
	Eigen::MatrixXd gamma_k;

	Eigen::MatrixXd wpc;
	Eigen::VectorXi wpc_stat;
	Eigen::VectorXd wpc_length;

	Eigen::MatrixXd wpc_pert;
	Eigen::VectorXi wpc_stat_pert;
	Eigen::VectorXd wpc_length_pert;

	Eigen::MatrixXd wpdc;
	Eigen::VectorXi wpdc_stat;
	Eigen::VectorXd wpdc_length;
	
	Eigen::MatrixXd Eij;

	std::vector < std::shared_ptr<Joint> > joints;
	std::vector < std::shared_ptr<RBState> > bodies;
	std::vector < std::shared_ptr<Contacts> > contacts;
	std::vector < std::shared_ptr<Spring> > springs;
	std::vector < std::shared_ptr<Cylinder> > cylinders;
	std::vector < std::shared_ptr<DoubleCylinder> > doublecylinders;

	std::vector < std::shared_ptr<Particle> > Ps;
	std::vector < std::shared_ptr<Particle> > Ss;
	std::vector < std::shared_ptr<Particle> > Os;

	std::vector < std::shared_ptr<Particle> > dPs;
	std::vector < std::shared_ptr<Particle> > dSs;
	std::vector < std::shared_ptr<Particle> > dUs;
	std::vector < std::shared_ptr<Particle> > dVs;

	std::vector < int > colList;

	static const int BALL_JOINT = 3;
	static const int HINGE_JOINT = 5;
	static const int HINGE_JOINT_X = 0;
	static const int HINGE_JOINT_Y = 1;
	static const int HINGE_JOINT_Z = 2;

	Eigen::Vector3d test_o;
	Eigen::Vector3d test_z;
	Eigen::MatrixXd test_pts;		// 3 x numFinitePoints
	Eigen::MatrixXd test_pts_pert;	// 3 x numFinitePoints

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