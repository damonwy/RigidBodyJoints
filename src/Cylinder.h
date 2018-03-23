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

	std::shared_ptr<Particle> P;	// start point
	std::shared_ptr<Particle> S;	// end point
	std::shared_ptr<Particle> O;	// the origin 
	Eigen::Vector3d Z;				// the z axis
};

#endif // RIGIDBODYJOINTS_SRC_CYLINDER_H_