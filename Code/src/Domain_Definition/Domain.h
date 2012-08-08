/*
 * Domain.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef DOMAIN_H_
#define DOMAIN_H_

#include "../Utilitaries/Pointcpp.h"

class Domain {
public:
	pointcpp<double> min_pt;
	pointcpp<double> max_pt;
public:
	Domain();
	virtual ~Domain();
	Domain(double,double);
	bool on_input_limit(pointcpp<double>);
	bool on_output_limit(pointcpp<double>);
	double domain_size_x(){return max_pt.i-min_pt.i;}
	double domain_size_y(){return max_pt.j-min_pt.j;}
	void get_extremities(pointcpp<double> & pmin, pointcpp<double> & pmax){pmin=min_pt; pmax=max_pt;}
};

#endif /* DOMAIN_H_ */
