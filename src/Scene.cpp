#include <iostream>

#include "Scene.h"
#include "Particle.h"
//#include "Cloth.h"
#include "Shape.h"
#include "Program.h"
#include "RigidBody.h"
#include "Cylinder.h"
#include "DoubleCylinder.h"

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
	numDoubleCylinders = rigidbody->numDoubleCylinders;

	Matrix3d R;
	R << 0, 0, 1, 1, 0, 0, 0, 1, 0;

	for (int i = 0; i < numCylinders; i++) {

		auto cylinderShape = make_shared<Shape>();
		cylinderShapes.push_back(cylinderShape);

		cylinderShape->loadMesh(RESOURCE_DIR + "cylinder.obj");
		cylinderShape->rotate(R);

		auto cylinder = make_shared<Particle>(cylinderShapes[i]);

		cylinders.push_back(cylinder);
		cylinder->r = rigidbody->cylinders[i]->r;
		cylinder->x = rigidbody->cylinders[i]->O->x;
	}


	for (int i = 0; i < numDoubleCylinders; i++) {

		auto dcShape0 = make_shared<Shape>();
		auto dcShape1 = make_shared<Shape>();
		dcShapes.push_back(dcShape0);
		dcShapes.push_back(dcShape1);

		dcShape0->loadMesh(RESOURCE_DIR + "cylinder.obj");
		dcShape0->rotate(R);
		dcShape1->loadMesh(RESOURCE_DIR + "cylinder.obj");
		dcShape1->rotate(R);
		auto dc0 = make_shared<Particle>(dcShapes[2*i + 0]);
		auto dc1 = make_shared<Particle>(dcShapes[2*i + 1]);
		doublecylinders.push_back(dc0);
		doublecylinders.push_back(dc1);
		dc0->r = rigidbody->doublecylinders[i]->Ur;
		dc0->x = rigidbody->doublecylinders[i]->U->x;
		dc1->r = rigidbody->doublecylinders[i]->Vr;
		dc1->x = rigidbody->doublecylinders[i]->V->x;
	}
}

void Scene::init()
{
	rigidbody->init();
	for (int i = 0; i < numCylinders; i++) {
		cylinderShapes[i]->init();
	}
	for (int i = 0; i < 2*numDoubleCylinders; i++) {
		dcShapes[i]->init();
	}
	
}

void Scene::tare()
{
	
	for (int i = 0; i < numCylinders; i++) {
		cylinders[i]->tare();
	}

	for (int i = 0; i < 2 * numCylinders; i++) {
		doublecylinders[i]->tare();
	}
}

void Scene::reset()
{
	t = 0.0;
	for (int i = 0; i < numCylinders; i++) {
		cylinders[i]->reset();
	}

	for (int i = 0; i < 2*numCylinders; i++) {
		doublecylinders[i]->reset();
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
	if (!doublecylinders.empty()) {
		for (int i = 0; i < numDoubleCylinders; i++) {
			doublecylinders[2*i + 0]->x = rigidbody->doublecylinders[i]->U->x;
			doublecylinders[2*i + 1]->x = rigidbody->doublecylinders[i]->V->x;
		}
	}
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> prog2, shared_ptr<MatrixStack> P) const
{
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
	
	for (int i = 0; i < numCylinders; i++) {
		
		cylinders[i]->draw(MV, prog);
	}

	for (int i = 0; i < 2*numDoubleCylinders; i++) {

		doublecylinders[i]->draw(MV, prog);
	}

	rigidbody->draw(MV, prog, prog2, P);


}
