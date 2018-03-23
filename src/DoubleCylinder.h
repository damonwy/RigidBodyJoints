#pragma once
#ifndef RIGIDBODYJOINTS_SRC_DOUBLECYLINDER_H_
#define RIGIDBODYJOINTS_SRC_DOUBLECYLINDER_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Shape;
class Program;
class MatrixStack;
class Particle;

class DoubleCylinder
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		DoubleCylinder(std::shared_ptr<Particle> _P, std::shared_ptr<Particle> _S, std::shared_ptr<Particle> _U, std::shared_ptr<Particle> _V);

	virtual ~DoubleCylinder();

	double Ur; // radius
	double Vr;
	double L; // initial length
	double l; // current length
	double stiffness;
	std::shared_ptr<Particle> P;	// start point
	std::shared_ptr<Particle> S;	// end point
	std::shared_ptr<Particle> U;	// the origin of U
	std::shared_ptr<Particle> V;	// the origin of V 
	Eigen::Vector3d Zu;				// the z axis of U
	Eigen::Vector3d Zv;				// the z axis of V

	Eigen::Vector3d pdir;			// the direction of p
	Eigen::Vector3d sdir;			// the direction of s

	Eigen::Vector3d fp;				// forces along P
	Eigen::Vector3d fs;				// forces along S
};

#endif // RIGIDBODYJOINTS_SRC_DOUBLECYLINDER_H_