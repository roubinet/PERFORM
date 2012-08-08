/*
 * Constantes.h
 *
 *  Created on: Jul 2, 2012
 *      Author: delphine
 */

#ifndef CONSTANTES_H_
#define CONSTANTES_H_

#include "Pointcpp.h"
#include <math.h>

#ifndef EPSILON
#define EPSILON 1E-6		//epsilon value
#endif

static const double PI=4.*atan(1.0);

typedef std::map<pointcpp<double>,std::pair<int,double> > Projection_Map;


const int SCALE_LIN = 0;
const int SCALE_LOG = 1;

#endif /* CONTANTES_H_ */
