/*
 * FractureMesh.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef FRACTUREMESH_H_
#define FRACTUREMESH_H_

#include "../Utilitaries/FluxPoint2D.h"
#include "../Utilitaries/Pointcpp.h"
#include "../Utilitaries/RandomNumber.h"
#include "../Utilitaries/Segment.h"
#include "../Utilitaries/Constantes.h"


// definition of a mesh
class FractureMesh{
public:
	FluxPoint2D p_ori; /**< origin of the mesh. */
	FluxPoint2D p_tar;	/**< target of the mesh. */
	double velocity;	/**< the velocty in the mesh */
	double aperture;	// mesh aperture
	int mesh_index; /**< the index of the mesh .*/
	std::vector<int> neigh_meshes;	/**< the neighboring meshes.*/
public:
	FractureMesh(){};
	virtual ~FractureMesh(){};
	FractureMesh(FluxPoint2D,FluxPoint2D,double,double,int,std::vector<int>);
	double distance_from_M_to_extremity(pointcpp<double>, pointcpp<double> &);
	double Mesh_Time_Advection(pointcpp<double>,pointcpp<double>);
	double Get_Total_Time_From_Advec_Time(double,double,double,RngStream_a &);
	double Advection_Time_Computation(const pointcpp<double> &,const Projection_Map&,const double,const double,const double,const double);
	int Is_Collinear(CgalPoint2D);
	int Is_Intersected(Segment2D&,CgalPoint2D&);
	pointcpp<double> velocity_interpolation(const pointcpp<double>&);
	pointcpp<double> Mesh_Scale_Advection(pointcpp<double>,double,bool opposite=false);
};


#endif /* FRACTUREMESH_H_ */
