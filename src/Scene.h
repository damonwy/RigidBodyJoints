#pragma once
#ifndef RIGIDBODYJOINTS_SRC_SCENE_H_
#define RIGIDBODYJOINTS_SRC_SCENE_H_

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Particle;
class MatrixStack;
class Program;
class Shape;
class RigidBody;
class Cylinder;
class DoubleCylinder;

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Scene();
	virtual ~Scene();
	
	void load(const std::string &RESOURCE_DIR);
	void init();
	void tare();
	void reset();
	void step();
	
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;
	
	double getTime() const { return t; }
	
private:
	double t;
	double h;
	Eigen::Vector3d grav;

	int numCylinders;
	int numDoubleCylinders;
	
	std::shared_ptr<RigidBody> rigidbody;

	std::vector< std::shared_ptr<Shape> > cylinderShapes;
	std::vector< std::shared_ptr<Shape> > dcShapes;
	std::vector< std::shared_ptr<Particle> > cylinders;
	std::vector < std::shared_ptr<Particle> > doublecylinders;
};

#endif // RIGIDBODYJOINTS_SRC_SCENE_H_
