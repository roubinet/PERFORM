/*
 * Parameters.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

enum{INFINITE_MATRIX,FINITE_MATRIX};

#include <string>

class Parameters{
public:
	double Lx;	// domain size in the longitudinal direction
	double Ly;	// domain size in the transversal direction
	double Dm;	// matrix diffusion coefficient
	double porosity;	// matrix porosity
	int nb_part;	// number of particles for transport simulation
	double proba_transfer;	// determine segment discretization for particle transport
	int simu_option;	// 0/1 = INFINITE_MATRIX/FINITE_MATRIX
	int Nt;	// number of time steps for result post-processing
	std::string code_path;
	std::string file_name_param;
	std::string file_name_domain;

public:
	Parameters();
	~Parameters(){};
	void read_param();
};



#endif /* PARAMETERS_H_ */
