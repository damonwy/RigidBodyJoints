#pragma once
#ifndef RIGIDBODYJOINTS_SRC_PARTICLE_H_
#define RIGIDBODYJOINTS_SRC_PARTICLE_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Shape;
class Program;
class MatrixStack;

class Particle
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Particle();
	Particle(const std::shared_ptr<Shape> shape);
	virtual ~Particle();
	void tare();
	void reset();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	
	double r; // radius
	double m; // mass
	int i;  // starting index
	Eigen::Vector3d x0; // initial position
	Eigen::Vector3d v0; // initial velocity
	Eigen::Vector3d x;  // world position
	Eigen::Vector3d xi; // local position
	Eigen::Vector3d xo; 
	Eigen::Vector3d v;  // velocity

	Eigen::Vector3f P;
	Eigen::Vector3f S;
	Eigen::Vector3f O;
	Eigen::Vector3f Z;
	
	int rb_id;

	bool fixed;
	
private:
	const std::shared_ptr<Shape> sphere;
};

#endif // RIGIDBODYJOINTS_SRC_PARTICLE_H_
