#include <iostream>

#include "Scene.h"
#include "Particle.h"
//#include "Cloth.h"
#include "Shape.h"
#include "Program.h"
#include "RigidBody.h"

using namespace std;
using namespace Eigen;

Scene::Scene() :
	t(0.0),
	h(1e-3),
	grav(0.0, 0.0, 0.0)
{
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{
	// Units: meters, kilograms, seconds
	h = 1e-2;
	
	grav << 0.0, -9.8, 0.0;
	
	rigidbody = make_shared<RigidBody>();
	numCylinders = rigidbody->numCylinders;

	cylinderShape = make_shared<Shape>();
	cylinderShape->loadMesh(RESOURCE_DIR + "cylinder.obj");
	cylinderShape->resize(4.0);

	Matrix3d R;
	R << 0, 0, 1, 1, 0, 0, 0, 1, 0;
	cylinderShape->rotate(R);

	VectorXd init_cyl_x;
	init_cyl_x.resize(numCylinders * 3);
	init_cyl_x.setZero();
	init_cyl_x.segment<3>(0) = Vector3d(1.0, 7.2, -0.2);

	for (int i = 0; i < numCylinders; i++) {
		auto cylinder = make_shared<Particle>(cylinderShape);

		cylinders.push_back(cylinder);
		cylinder->r = 0.1;
		cylinder->x = init_cyl_x.segment<3>(3 * i);
	}
	

}

void Scene::init()
{
	rigidbody->init();
	cylinderShape->init();
}

void Scene::tare()
{
	
	for (int i = 0; i < numCylinders; i++) {
		cylinders[i]->tare();
	}
}

void Scene::reset()
{
	t = 0.0;
	for (int i = 0; i < numCylinders; i++) {
		cylinders[i]->reset();
	}
}

void Scene::step()
{
	t += h;

	rigidbody->step(h);

	// move the cylinders
	//if (!cylinders.empty()) {
	//	auto c = cylinders.front();
	//	c->x = // TODO
	//}
	//TODO
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> prog2, shared_ptr<MatrixStack> P) const
{
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
	
	for (int i = 0; i < numCylinders; i++) {
		
		cylinders[i]->draw(MV, prog);
	}

	rigidbody->draw(MV, prog, prog2, P);


}
