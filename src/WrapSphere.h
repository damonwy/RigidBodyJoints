#ifndef WrapSphere_H
#define WrapSphere_H

/*
* WrapSphere.hpp
*
* Obstacle Set Algorithm Simulation for Sphere Obstacles
*
*/

#include "WrapObst.h"

class WrapSphere : public WrapObst
{
public:
	// default constructor
	WrapSphere()
	{
		type = sphere;
	}

	void setSphereConfig(const Eigen::Vector3d &O)
	{
		this->point_O = O;
	}

	// constructor
	WrapSphere(const Eigen::Vector3d &P,
		const Eigen::Vector3d &S,
		const Eigen::Vector3d &O,
		const double R)
		: WrapObst(P, S, O, R)
	{
		type = sphere;
	}

	using WrapObst::compute;
	void compute();

	using WrapObst::getPoints;
	Eigen::MatrixXd getPoints(int num_points);
};

#endif