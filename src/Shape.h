#pragma once
#ifndef RIGIDBODYJOINTS_SRC_SHAPE_H_
#define RIGIDBODYJOINTS_SRC_SHAPE_H_

#include <string>
#include <vector>
#include <memory>
#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Program;

/**
 * A shape defined by a list of triangles
 * - posBuf should be of length 3*ntris
 * - norBuf should be of length 3*ntris (if normals are available)
 * - texBuf should be of length 2*ntris (if texture coords are available)
 * posBufID, norBufID, and texBufID are OpenGL buffer identifiers.
 */
class Shape
{
public:
	Shape();
	virtual ~Shape();
	void loadMesh(const std::string &meshName);
	void init();
	void draw(const std::shared_ptr<Program> prog) const;
	void resize(double ratio);
	void rotate(Eigen::Matrix3d R);
	
private:
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;
};

#endif // RIGIDBODYJOINTS_SRC_SHAPE_H_
