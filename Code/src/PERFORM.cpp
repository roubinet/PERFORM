//============================================================================
// Name        : PERFORM.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : C++ code for solute transport in heterogeneous fractured porous media
//============================================================================

#include "Input_Output/Parameters.h"
#include "Domain_Definition/Domain.h"
#include "Domain_Definition/NetworkMeshes.h"
#include "Transport/Transport.h"
#include "Input_Output/Results.h"
#include "Utilitaries/Structures.h"
#include <iostream>

using namespace std;

int main() {

	cout << "PERFORM begin" << endl;

	// 1. Parameter and domain definition reading
	// 1.1. Parameter reading
	Parameters param;
	// 1.2. Domain definition
	Domain domain(param.Lx,param.Ly);
	NetworkMeshes net_mesh(param.code_path,param.file_name_domain);

	// 2. Transport modeling
	RngStream_a rng_tracker;
	Transport transport(rng_tracker,net_mesh,domain,param);
	map<int,double> arrival_times = transport.Particles_Transport();

	// 3. Results post-processing
	Results results(arrival_times,param);
	results.post_processing();
	results.writing();

	cout << "PERFORM end" << endl;

	return 0;
}
