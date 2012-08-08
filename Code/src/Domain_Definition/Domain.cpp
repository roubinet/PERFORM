/*
 * Domain.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#include "Domain.h"
#include "../Utilitaries/Constantes.h"

using namespace std;

Domain::Domain(){}
Domain::~Domain(){}

Domain::Domain(double Lx, double Ly){
	min_pt = pointcpp<double>(-0.5*Lx,-0.5*Ly);
	max_pt = pointcpp<double>(0.5*Lx,0.5*Ly);
}

bool Domain::on_input_limit(pointcpp<double> pt){
	if ((pt.j>max_pt.j)||(pt.j<min_pt.j)){
		cout << "WARNING in on_input_limit (Domain.cpp): point out of the domain" << endl;
		return false;
	}
	if (abs(pt.i-min_pt.i)<EPSILON){return true;}
	return false;
}

bool Domain::on_output_limit(pointcpp<double> pt){
	if ((pt.j>max_pt.j)||(pt.j<min_pt.j)){
		cout << "WARNING in on_input_limit (Domain.cpp): point out of the domain" << endl;
		return false;
	}
	if (abs(pt.i-max_pt.i)<EPSILON){return true;}
	return false;
}
