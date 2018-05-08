#pragma once
#ifndef RIGIDBODYJOINTS_SRC_SPRING_H_
#define RIGIDBODYJOINTS_SRC_SPRING_H_

#include <memory>
#include <vector>
#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Particle;
class RBState;

enum SpringType { one_end_fixed, two_end_fixed, two_end_rbs };

class Spring
{
public:
	Spring(std::shared_ptr<Particle> p0, std::shared_ptr<Particle> p1);
	void updateInertiaBlocks();
	void updateDensity();
	void computeJacobians(double epsilon, const std::vector <std::shared_ptr<RBState> > bodies, Eigen::Vector3d g);
	double getPotentialEnergy(Eigen::Vector3d g);
	double getKineticEnergy(const std::vector <std::shared_ptr<RBState> > bodies);
	void computeGravityWrench(Eigen::Vector3d g);

	virtual ~Spring();
	int i;
	int k;
	int in;
	int kn;
	std::shared_ptr<Particle> p0;
	std::shared_ptr<Particle> p1;
	double E;
	double L;
	SpringType type;
	double mass;
	double muscle_density;

	Eigen::MatrixXd IM;  // 12x12  inertia matrix
	Eigen::MatrixXd I11; // 6x6
	Eigen::MatrixXd I12; 
	Eigen::MatrixXd I21;
	Eigen::MatrixXd I22;

	// The approximated material Jacobian matrix of rb0, rb1
	Eigen::MatrixXd J0; // 
	Eigen::MatrixXd J1;

	double KE;
	double PE;
	Eigen::VectorXd wrench;//12x1

};

#endif // RIGIDBODYJOINTS_SRC_SPRING_H_
