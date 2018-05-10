#pragma once
#ifndef RigidBodyJoints_SRC_JSONEIGEN_H_
#define RigidBodyJoints_SRC_JSONEIGEN_H_

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

#include <json.hpp>



namespace nlohmann {
	template <>
	struct adl_serializer<Eigen::ArrayXd> {
		static void to_json(json& j, const Eigen::ArrayXd& v) {
			for (unsigned int i = 0; i < v.size(); i++) {
				j[i] = v[i];
			}
		}

		static void from_json(const json& j, Eigen::ArrayXd& v) {
			unsigned int s = j.size();

			v.resize(s);
			for (unsigned int i = 0; i < s; i++) {
				v[i] = j[i];
			}
		}
	};

	template <>
	struct adl_serializer<Eigen::VectorXi> {
		static void to_json(json& j, const Eigen::VectorXi& v) {
			for (unsigned int i = 0; i < v.size(); i++) {
				j[i] = v[i];
			}
		}

		static void from_json(const json& j, Eigen::VectorXi& v) {
			unsigned int s = j.size();

			v.resize(s);
			for (unsigned int i = 0; i < s; i++) {
				v[i] = j[i];
			}
		}
	};


	//template <>
	//struct adl_serializer<MatrixXd> {
	//	static void to_json(json& j, const MatrixXd& m) {
	//		for (unsigned int i = 0; i < m.rows(); i++) {
	//			j[i] = json();
	//			for (unsigned int j = 0; j < m.cols(); j++) {
	//				j[i][j] = m(i, j);
	//			}
	//		}
	//	}
	//	//static void from_json(const json&j,)
	//};
}




#endif  // RigidBodyJoints_SRC_JSONEIGEN_H_