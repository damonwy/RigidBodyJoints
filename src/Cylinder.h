#pragma once
#ifndef RIGIDBODYJOINTS_SRC_CYLINDER_H_
#define RIGIDBODYJOINTS_SRC_CYLINDER_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Shape;
class Program;
class MatrixStack;
class Particle;

class Cylinder
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Cylinder(std::shared_ptr<Particle> _P, std::shared_ptr<Particle> _S, std::shared_ptr<Particle> _O);
	
	virtual ~Cylinder();

	double r; // radius
	double L; // initial length
	double l; // current length

	double l_pc1;	// length of pc1
	double l_c1c2;	// length of c1c2
	double l_c2s;	// length of c2s

	double stiffness;
	std::shared_ptr<Particle> P;	// start point
	std::shared_ptr<Particle> S;	// end point
	std::shared_ptr<Particle> O;	// the origin 
	Eigen::Vector3d Z;				// the z axis

	Eigen::Vector3d pdir;			// the direction of p
	Eigen::Vector3d sdir;			// the direction of s

	Eigen::Vector3d fp;
	Eigen::Vector3d fs;
	Eigen::Vector3d fpc;
	Eigen::Vector3d fsc;

	Eigen::Vector3d c1;				// contact point 1
	Eigen::Vector3d c2;

	int status;  // 1: wrap 0: no wrap
	double mass;

	// used to compute the arc
	double theta_s;
	double theta_e;
	Eigen::Matrix3d M;
};

#endif // RIGIDBODYJOINTS_SRC_CYLINDER_H_