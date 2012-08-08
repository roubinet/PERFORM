/*
 * Results.h
 *
 *  Created on: Jun 29, 2012
 *      Author: delphine
 */

#ifndef RESULTS_H_
#define RESULTS_H_

#include "../Input_Output/Parameters.h"
#include <map>
#include <string>


class Results{
public:
	int Nt;	// number of time steps for results post-processing
	std::map<int,double> arrival_times;	//<particle identifier,arrival time>
	std::map<double,int> cum_dist_times;	//<arrival time,number of particles arrived before this time>
	std::string output_file;	// file name to write results
public:
	Results(){};
	virtual ~Results(){};
	Results(std::map<int,double>,Parameters);
	void post_processing();
	void writing();
};


#endif /* RESULTS_H_ */
