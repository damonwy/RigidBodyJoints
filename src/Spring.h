#pragma once
#ifndef RIGIDBODYJOINTS_SRC_SPRING_H_
#define RIGIDBODYJOINTS_SRC_SPRING_H_

#include <memory>

class Particle;

class Spring
{
public:
	Spring(std::shared_ptr<Particle> p0, std::shared_ptr<Particle> p1);
	virtual ~Spring();
	int i;
	int k;
	int in;
	int kn;
	std::shared_ptr<Particle> p0;
	std::shared_ptr<Particle> p1;
	double E;
	double L;
};

#endif // RIGIDBODYJOINTS_SRC_SPRING_H_
