#include "Spring.h"
#include "Particle.h"
#include "Helper.h"
#include "RBState.h"
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;
using namespace Eigen;
using namespace Helper;

Spring::Spring(shared_ptr<Particle> p0, shared_ptr<Particle> p1) :
	E(1.0)
{
	assert(p0);
	assert(p1);
	assert(p0 != p1);
	this->p0 = p0;
	this->p1 = p1;
	Vector3d x0 = p0->x;
	Vector3d x1 = p1->x;
	Vector3d dx = x1 - x0;
	L = dx.norm();
	updateDensity();
	//assert(L > 0.0);

	IM.resize(12, 12);
	I11.resize(6, 6);
	I12.resize(6, 6);
	I21.resize(6, 6);
	I22.resize(6, 6);

	J0.resize(3, 6);
	J1.resize(3, 6);

	IM.setZero();
	I11.setZero();
	I12.setZero();
	I21.setZero();
	I22.setZero();

	J0.setZero();
	J1.setZero();

	wrench.resize(12);
}

void Spring::updateInertiaBlocks() {
	IM.block<6, 6>(0, 0) = I11;
	IM.block<6, 6>(0, 6) = I12;
	IM.block<6, 6>(6, 0) = I21;
	IM.block<6, 6>(6, 6) = I22;
}

void Spring::updateDensity() {
	muscle_density = mass / (p1->x - p0->x).norm();
}

void Spring::computeJacobians(double epsilon, const vector<shared_ptr<RBState> > bodies, Vector3d g) {
	auto b0 = bodies[i];
	auto b1 = bodies[k];
	updateDensity();

	VectorXd pert;
	pert.resize(6);

	// for each component of phi(i = 0, 1, 2..,11) add a relative small perturbation
	for (int k = 0; k < 6; k++) {

		pert.setZero();
		pert(k) = 1.0 * epsilon; // change kth component

		MatrixXd E_pert = b0->E * vec2crossmatrix(pert).exp();

		// Compute position of the end point p0 of spring
		Vector3d p_pert = transform(E_pert, p0->x0);
		J0.block(0, k, 3, 1) = 1.0 / epsilon * (p_pert - p0->x);
	}

	for (int k = 0; k < 6; k++) {

		pert.setZero();
		pert(k) = 1.0 * epsilon;

		MatrixXd E_pert = b1->E * vec2crossmatrix(pert).exp();

		Vector3d p_pert = transform(E_pert, p1->x0);
		J1.block(0, k, 3, 1) = 1.0 / epsilon * (p_pert - p1->x);
	}

	I11 = J0.transpose() * J0 * 1.0 / 3.0 * muscle_density;
	I12 = J0.transpose() * J1 * 1.0 / 6.0 * muscle_density;
	I21 = J1.transpose() * J0 * 1.0 / 6.0 * muscle_density;
	I22 = J1.transpose() * J1 * 1.0 / 3.0 * muscle_density;

	updateInertiaBlocks();
	computeGravityWrench(g);

}

double Spring::getPotentialEnergy(Vector3d g) {
	
	PE = 0.5 * mass * g.transpose() * (p0->x + p1->x);
	return PE;
}

double Spring::getKineticEnergy(const vector<shared_ptr<RBState> > bodies) {
	// Compute the kinetic energy of a muscle
	auto b0 = bodies[i];
	auto b1 = bodies[k];

	VectorXd PHIM;
	PHIM.resize(12);
	PHIM.setZero();
	PHIM.segment<6>(0) = b0->Phi;
	PHIM.segment<6>(6) = b1->Phi;

	KE = 0.5 * PHIM.transpose() * IM * PHIM;
	return KE;
}


void Spring::computeGravityWrench(Vector3d g) {
	// Add gravity force 
	Vector3d fg = mass * g;
	MatrixXd Jt;
	Jt.resize(3, 12);
	Jt.block<3, 6>(0, 0) = J0;
	Jt.block<3, 6>(0, 6) = J1;
	wrench = Jt.transpose() * fg; // 12x1
	
}

Spring::~Spring()
{
	
}

