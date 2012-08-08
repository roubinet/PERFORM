/*
 * FractureMesh.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */


#include "FractureMesh.h"
#include <boost/math/special_functions/erf.hpp>

using namespace std;


FractureMesh::FractureMesh(FluxPoint2D p_ori_,FluxPoint2D p_tar_,double velocity_,double aperture_,int mesh_index_,std::vector<int> neigh_meshes_){
	p_ori = p_ori_;
	p_tar = p_tar_;
	velocity = velocity_;
	aperture = aperture_;
	mesh_index = mesh_index_;
	neigh_meshes = neigh_meshes_;
}

// return the distance from the current position M to the extremity for a displacement with the velocity
double FractureMesh::distance_from_M_to_extremity(pointcpp<double> M, pointcpp<double> & M_out){

	double distance = 0;
	CgalPoint2D M_(M.i,M.j);
	if (velocity>0){
		distance = sqrt(CGAL::squared_distance(M_,p_tar.p));
		M_out.define((double)p_tar.p.x(),(double)p_tar.p.y());
	}
	else if (velocity<0){
		distance = sqrt(CGAL::squared_distance(M_,p_ori.p));
		M_out.define((double)(p_ori.p.x()),(double)(p_ori.p.y()));
	}
	return distance;

}

// return the time required to travel by advection from M_in to M_out in the current mesh
double FractureMesh::Mesh_Time_Advection(pointcpp<double> M_in, pointcpp<double> M_out){

	if (velocity == 0){
		cout << "WARNING in Mesh_Time_Advection (FractureMesh.cpp): velocity null" << endl;
	}

	CgalPoint2D M_in_(M_in.i,M_in.j), M_out_(M_out.i,M_out.j);
	double distance = sqrt(CGAL::squared_distance(M_in_,M_out_));
	return distance/velocity;
}

/**returns the total time spent in a segment from the advection time and considering diffusion in the infinite surrounding matrix*/
// Dm and porosity are matrix diffusion coefficient and matrix porosity
double FractureMesh::Get_Total_Time_From_Advec_Time(double advection_time, double Dm, double porosity, RngStream_a & rng_tracker){

	double B = sqrt(Dm)*porosity*advection_time/(this->aperture*0.5);
	double u = rng_tracker.uniform(0,1);
	double value = boost::math::erfc_inv(u);
	double time = advection_time+0.25*B/value*B/value;
	return time;

}

/** returns the advection time required for the scheme stability (transfer probability lower than transfer_proba_max)*/
double FractureMesh::Advection_Time_Computation(const pointcpp<double> & current_position, const Projection_Map & trans_positions, const double porosity, const double diffusion, const double transfer_proba_max, const double t_advec_end){

	//case without restriction
	if (transfer_proba_max == -1) return t_advec_end;
	//0. PARAMETERS
	if ((trans_positions.size()<0)||(trans_positions.size()>2)){
		std::cout<<"\n PB IN Advection_Time_Computation (Utilitary_Transfer.cpp) : res not defined";
		return -1;
	}

	//1. NO SURROUNDING FRACTURES
	if (trans_positions.size()==0) return t_advec_end;

	//2. ONE OR TWO SURROUNDING FRACTURES
	//barriers parameters
	double distance1 = 0, distance2 = 0;
	if (trans_positions.size()==1){
		distance1 = trans_positions.begin()->first.Point_Distance(current_position);
		distance2 = distance1;
	}
	else if (trans_positions.size()==2){
		Projection_Map::const_iterator it = trans_positions.begin();
		distance1 = it->first.Point_Distance(current_position);
		it++;
		distance2 = it->first.Point_Distance(current_position);
	}
	double mean_transfer_time = distance1*distance2/(2*diffusion);
	double advection_time = aperture*std::sqrt(mean_transfer_time)*boost::math::erfc_inv(1-transfer_proba_max)/(porosity*std::sqrt(diffusion));
	return std::min(advection_time, t_advec_end);
}

pointcpp<double> FractureMesh::Mesh_Scale_Advection(pointcpp<double> M_in, double dt, bool opposite){
	pointcpp<double> v_inter = velocity_interpolation(M_in);
	if (opposite) v_inter = pointcpp<double>(-v_inter.i,-v_inter.j,-v_inter.k);
	return M_in + v_inter * dt;
}

/// Gets the interpolated velocity at point M
pointcpp<double> FractureMesh::velocity_interpolation(const pointcpp<double> & M){
	pointcpp<double> res(0.,0.,0.);
	CgalVector2D vect(p_ori.p, p_tar.p);
	res.i = vect.x()/CGAL::sqrt(vect.squared_length()) * velocity;
	res.j = vect.y()/CGAL::sqrt(vect.squared_length()) * velocity;
	return res;
}

int FractureMesh::Is_Collinear(CgalPoint2D pt){
	if (CGAL::collinear(p_ori.p,p_tar.p,pt)) return 1;
	return 0;
}

//returns 1 if the current mesh intersects the segment seg
//inter is the intersection point
int FractureMesh::Is_Intersected(Segment2D & seg, CgalPoint2D & inter){
	Segment2D seg_mesh(p_ori.p,p_tar.p);
	if (seg_mesh.intersection(seg,inter)) return 1;
	return 0;
}


