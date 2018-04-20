#include "WrapCylinder.h"

void WrapCylinder::compute()
{
	Eigen::Vector3d OP = this->point_P - this->point_O;
	OP = OP / OP.norm();
	Eigen::Vector3d vec_Z = vec_z / vec_z.norm();
	Eigen::Vector3d vec_X = vec_Z.cross(OP);
	vec_X = vec_X / vec_X.norm();
	Eigen::Vector3d vec_Y = vec_Z.cross(vec_X);
	vec_Y = vec_Y / vec_Y.norm();

	this->M << vec_X.transpose(), vec_Y.transpose(), vec_Z.transpose();

	Eigen::Vector3d p = this->M * (this->point_P - this->point_O);
	Eigen::Vector3d s = this->M * (this->point_S - this->point_O);

	double denom_q = p(0)*p(0) + p(1)*p(1);
	double denom_t = s(0)*s(0) + s(1)*s(1);
	double R = this->radius;

	if ((denom_q - R*R < 0.0f) || (denom_t - R*R < 0.0))
	{
		this->status = inside_radius;
	}

	double root_q = sqrt(denom_q - R*R);
	double root_t = sqrt(denom_t - R*R);

	Eigen::Vector3d q(0.0, 0.0, 0.0);
	Eigen::Vector3d t(0.0, 0.0, 0.0);
	q(0) = (p(0) * R*R + R * p(1) * root_q) / denom_q;
	q(1) = (p(1) * R*R - R * p(0) * root_q) / denom_q;
	t(0) = (s(0) * R*R - R * s(1) * root_t) / denom_t;
	t(1) = (s(1) * R*R + R * s(0) * root_t) / denom_t;

	if (R * (q(0) * t(1) - q(1) * t(0)) > 0.0)
	{
		this->status = no_wrap;
	}

	this->status = wrap;

	std::complex<double> qt_i = 1.0 - 0.5 *
		((q(0) - t(0)) * (q(0) - t(0))
			+ (q(1) - t(1)) * (q(1) - t(1))) / (R*R);
	double qt_xy = abs(R * acos(qt_i));// changed
	this->path_length = qt_xy;

	double pq_xy = sqrt((p(0) - q(0)) * (p(0) - q(0)) +
		(p(1) - q(1)) * (p(1) - q(1)));
	double ts_xy = sqrt((t(0) - s(0)) * (t(0) - s(0)) +
		(t(1) - s(1)) * (t(1) - s(1)));
	q(2) = p(2) + (s(2) - p(2)) * pq_xy / (pq_xy + qt_xy + ts_xy);
	t(2) = s(2) - (s(2) - p(2)) * ts_xy / (pq_xy + qt_xy + ts_xy);

	this->point_q = q;
	this->point_t = t;

	Eigen::Vector3d Q = this->M.transpose() * q + this->point_O;
	Eigen::Vector3d T = this->M.transpose() * t + this->point_O;

	// std::cout << Q.transpose() << std::endl << T.transpose() << std::endl;
}

Eigen::MatrixXd WrapCylinder::getPoints(int num_points, double &theta_s, double &theta_e, Eigen::Matrix3d &_M)
{
	double theta_q = atan(this->point_q(1) / this->point_q(0));
	if (this->point_q(0) < 0.0)
		theta_q += PI;

	double theta_t = atan(this->point_t(1) / this->point_t(0));
	if (this->point_t(0) < 0.0)
		theta_t += PI;

	Eigen::MatrixXd points(3, num_points + 1);

	double z_s, z_e;
	//theta_s, theta_e
	_M = this->M;

	if (theta_q < theta_t)
	{
		theta_s = theta_q; theta_e = theta_t;
		z_s = this->point_q(2); z_e = this->point_t(2);
	}
	else
	{
		theta_s = theta_t; theta_e = theta_q;
		z_s = this->point_t(2); z_e = this->point_q(2);
	}

	if (theta_e - theta_s > theta_s + 2 * PI - theta_e)
	{
		double tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2 * PI;
		tmp = z_s; z_s = z_e; z_e = tmp;
	}

	int col = 0;
	double z_i = z_s, dz = (z_e - z_s) / num_points;
	for (double i = theta_s; i <= theta_e + 0.001;
		i += (theta_e - theta_s) / num_points)
	{	
		if (col == num_points + 1) {
			break;
		}

		Eigen::Vector3d point = this->M.transpose() *
			Eigen::Vector3d(this->radius * cos(i), this->radius * sin(i), z_i) +
			this->point_O;
		z_i += dz;
		points.col(col++) = point;
	}

	return points;
}