#ifndef WrapCylinder_H
#define WrapCylinder_H

/*
* WrapCylinder.hpp
*
* Obstacle Set Algorithm Simulation for Cylinder Obstacles
*
*/

#include "WrapObst.h"

class WrapCylinder : public WrapObst
{
private:
	Eigen::Vector3d vec_z;      // Cylinder Positive z axis

public:
	// default constructor
	WrapCylinder()
	{
		vec_z = Eigen::Vector3d(0.0, 0.0, 0.0);
		type = cylinder;
	}

	void setCylinderConfig(const Eigen::Vector3d &O,
		const Eigen::Vector3d &Z)
	{
		this->point_O = O;
		this->vec_z = Z;
	}

	// constructor
	WrapCylinder(const Eigen::Vector3d &P,
		const Eigen::Vector3d &S,
		const Eigen::Vector3d &O,
		const Eigen::Vector3d &Z,
		const double R)
		: WrapObst(P, S, O, R), vec_z(Z)
	{
		type = cylinder;
	}

	using WrapObst::compute;
	void compute();

	using WrapObst::getPoints;
	Eigen::MatrixXd getPoints(int num_points, double &theta_s, double &theta_e, Eigen::Matrix3d &_M);
};

#endif