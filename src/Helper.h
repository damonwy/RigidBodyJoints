#pragma once
#ifndef RIGIDBODYJOINTS_SRC_HELPER_H_
#define RIGIDBODYJOINTS_SRC_HELPER_H_
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

namespace Helper
{	
	Eigen::MatrixXd vec2crossmatrix(Eigen::VectorXd a); // repackage a vector into a cross-product matrix
	Eigen::VectorXd crossmatrix2vec(Eigen::MatrixXd A);
	Eigen::Vector3d transform(Eigen::MatrixXd E, Eigen::Vector3d x);
	Eigen::Vector3d transformVector(Eigen::MatrixXd E, Eigen::Vector3d vec);
	Eigen::MatrixXd computeAdjoint(Eigen::MatrixXd E);
}
	

#endif // RIGIDBODYJOINTS_SRC_HELPER_H_
