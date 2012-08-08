/*
 * Particle.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "../Utilitaries/Pointcpp.h"


class Particle {
public:
	pointcpp<double> M;	// current position of the particle
	int mesh_index;	// mesh index of the current location of the particle
	double t;	// current time
	int no;	// identifier of current particle
public:
	Particle();
	virtual ~Particle();
};

#endif /* PARTICLE_H_ */
