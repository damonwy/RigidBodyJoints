#include <iostream>

#include "Scene.h"
#include "Particle.h"
//#include "Cloth.h"
#include "Shape.h"
#include "Program.h"
#include "RigidBody.h"
#include "Cylinder.h"

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
	

	Matrix3d R;
	R << 0, 0, 1, 1, 0, 0, 0, 1, 0;

	for (int i = 0; i < numCylinders; i++) {

		auto cylinderShape = make_shared<Shape>();
		cylinderShapes.push_back(cylinderShape);

		cylinderShape->loadMesh(RESOURCE_DIR + "cylinder.obj");
		//cylinderShape->resize(1.0);
		cylinderShape->resize(rigidbody->cylinders[i]->r/0.4);
		cylinderShape->rotate(R);
		

		auto cylinder = make_shared<Particle>(cylinderShapes[i]);

		cylinders.push_back(cylinder);
		cylinder->r = rigidbody->cylinders[i]->r;
		cylinder->x = rigidbody->cylinders[i]->O->x;
	}
}

void Scene::init()
{
	rigidbody->init();
	for (int i = 0; i < numCylinders; i++) {
		cylinderShapes[i]->init();
	}
	
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
	if (!cylinders.empty()) {
		for (int i = 0; i < numCylinders; i++) {
			cylinders[i]->x = rigidbody->cylinders[i]->O->x;
		}
	}

}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> prog2, shared_ptr<MatrixStack> P) const
{
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
	
	for (int i = 0; i < numCylinders; i++) {
		
		cylinders[i]->draw(MV, prog);
	}

	rigidbody->draw(MV, prog, prog2, P);


}
