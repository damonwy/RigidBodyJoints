#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "DoubleCylinder.h"
#include "Particle.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"

using namespace std;

DoubleCylinder::DoubleCylinder(std::shared_ptr<Particle> _P, std::shared_ptr<Particle> _S, std::shared_ptr<Particle> _U, std::shared_ptr<Particle> _V) {
	this->P = _P;
	this->S = _S;
	this->U = _U;
	this->V = _V;
	this->L = -1.0;
	this->stiffness = 20;
}

DoubleCylinder::~DoubleCylinder() {

}