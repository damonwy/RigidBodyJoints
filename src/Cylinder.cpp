#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Cylinder.h"
#include "Particle.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"

using namespace std;

Cylinder::Cylinder(shared_ptr<Particle> _P, shared_ptr<Particle> _S, shared_ptr<Particle> _O)
{
	this->P = _P;
	this->S = _S;
	this->O = _O;
	this->L = -1.0;
	this->stiffness = 1e2;
}

Cylinder::~Cylinder()
{
}

