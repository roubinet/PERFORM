/*
 * Projection.cpp
 *
 *  Created on: Jul 2, 2012
 *      Author: delphine
 */


#include "Transport.h"


using namespace std;

// Projection of a point within the domain
void Transport::M_project_in_system(pointcpp<double> & M_out, const double & EPS){
	// if M is out the system, project it inside the system
	pointcpp<double> pmin,pmax;
	domain.get_extremities(pmin,pmax);
	if(M_out.i < pmin.i-EPS) M_out.i = pmin.i + EPS;
	if(M_out.i > pmax.i+EPS) M_out.i = pmax.i - EPS;
	if(M_out.j < pmin.j-EPS) M_out.j = pmin.j + EPS;
	if(M_out.j > pmax.j+EPS) M_out.j = pmax.j - EPS;
}

// Projection of a point within the domain
void Transport::M_project_in_system(CgalPoint2D & M_out, const double & EPS){
	// conversion from CgalPoint2D to pointcpp
	pointcpp<double> point(M_out);
	// projection in the system
	M_project_in_system(point,EPS);
	// conversion from pointcpp to CgalPoint2D
	M_out = point.CgalPoint();
}

/**
*returns the projection through the domain boundaries (consider periodic boundaries)
*when there is no projection on seg_ortho1 -> we use the projection by periodicity on seg_ortho2
*/
bool Transport::Projection_Periodic_Boundaries(CgalPoint2D init_point,Segment2D seg_ortho1,Segment2D seg_ortho2,CgalPoint2D M_end1,CgalPoint2D M_end2,CgalPoint2D & proj_point,int & proj_mesh,double & proj_dist){
	// 0. Parameters
	bool success = false;
	vector<FractureMesh> meshes = this->net_mesh.meshes;
	// 1. Verification of the closest projection
	if (Closest_Orthogonal_Projection(init_point,seg_ortho1,proj_point,proj_mesh,proj_dist)){
		std::cout<<"\n PB IN Projection_Periodic_Boundaries (Frac_Tracker.cpp): closest projection exists";
	}
	// 2. Determination of the projection through the periodic boundary
	else{
		// 2.1. Determination of the projection point -> the further from the initial point (closest to the second boundary)
		CgalPoint2D current_point;
		double max_dist = -1;
		int meshes_size = meshes.size();
		for (int i = 0; i<meshes_size; i++){
			if ((!meshes[i].Is_Collinear(init_point))&&(meshes[i].Is_Intersected(seg_ortho2, current_point))&&(!identic(init_point,current_point,EPSILON))&&(distance_2D(init_point,current_point)>max_dist)){
				proj_mesh = meshes[i].mesh_index;
				proj_point = current_point;
				max_dist = distance_2D(init_point,proj_point);
				success = true;
			}
		}
		// 2.2. Determination of the distance between the initial point to the projected point through the domain periodic boundaries
		// 2.2.1. Distance from the initial point to the closest boundary
		M_project_in_system(M_end1,EPSILON);
		double dist1 = distance_2D(init_point,M_end1);
		// 2.2.2. Distance from the second boundary to the projected point
		M_project_in_system(M_end2,EPSILON);
		double dist2 = distance_2D(M_end2,proj_point);
		// 2.2.3. total distance
		proj_dist = dist1+dist2;
	}
	return success;
}

/**returns the definition (point, mesh index and distance) of the closest projection of M orthogonal to seg_ortho*/
bool Transport::Closest_Orthogonal_Projection(CgalPoint2D init_point,Segment2D seg_ortho,CgalPoint2D & proj_point,int & proj_mesh,double & proj_dist){
	// 0. Parameters
	bool succes = false;
	vector<FractureMesh> meshes = this->net_mesh.meshes;
	int meshes_size = meshes.size();
	// 1. loop on the meshes to determine the closest intersection
	proj_dist = 1e9;
	CgalPoint2D current_point;
	for (int i = 0; i<meshes_size; i++){
		if ((!meshes[i].Is_Collinear(init_point))&&(meshes[i].Is_Intersected(seg_ortho, current_point))&&(!identic(init_point,current_point,EPSILON))&&(distance_2D(init_point,current_point)<proj_dist)){
			proj_mesh = meshes[i].mesh_index;
			proj_point = current_point;
			proj_dist = distance_2D(init_point,proj_point);
			succes = true;
		}
	}
	return succes;
}

/**
* returns the two closest projections of M1 orthogonal to the segment (M1,M2)
* map<projected point,pair<mesh index,distance from initial point to projected point>>
*/
Projection_Map Transport::Orthogonal_Projection_M(pointcpp<double> M, pointcpp<double> M2){
	//0. PARAMETERS
	Projection_Map proj_points;
	double L = max(this->domain.domain_size_x(),this->domain.domain_size_y());
	CgalPoint2D inter;

	//1. CONSTRUCTION OF THE SEGMENTS (ONE BY SIDE) ORTHOGONAL TO SEG AND STARTING FROM THE MIDDLE OF THE SEG
	//1.1. middle and angle of the initial segment
	CgalPoint2D M1 = M.CgalPoint(), M2_ = M2.CgalPoint();
	Segment2D seg(M1,M2_);
	double alpha = seg.get_orientation();	//angle of the current mesh
	//1.2. determination of the extremities of the orthogonal segments
	alpha += PI/2;
	CgalPoint2D M_end1 = CgalPoint2D(M1.x()+2*L*std::cos(alpha),M1.y()+2*L*std::sin(alpha));
	CgalPoint2D M_end2 = CgalPoint2D(M1.x()-2*L*std::cos(alpha),M1.y()-2*L*std::sin(alpha));
	Segment2D seg_ortho1(M1,M_end1);
	Segment2D seg_ortho2(M1,M_end2);

	//2. DETERMINATION OF THE INTERSECTIONS WITH THE FRACTURE NETWORK
	//2.0. Parameters to describe the selected projection
	double min_dist = 1e9;
	int mesh_index = -1;
	//2.1. determination of the projection : the closest one or the one through periodic boundary
	//2.1.1. first direction: from M to Mend1
	if (!Closest_Orthogonal_Projection(M1,seg_ortho1,inter,mesh_index,min_dist)){
		if (!Projection_Periodic_Boundaries(M1,seg_ortho1,seg_ortho2,M_end1,M_end2,inter,mesh_index,min_dist)){
			std::cout<<"\n PB IN Orthogonal_Projection_M (Projection.cpp): orthogonal projection not determined" << endl;
		}
	}
	proj_points[pointcpp<double>(inter)] = std::make_pair<int,double>(mesh_index,min_dist);
	//2.1.2. second direction: from M to Mend2
	if (!Closest_Orthogonal_Projection(M1,seg_ortho2,inter,mesh_index,min_dist)){
		if (!Projection_Periodic_Boundaries(M1,seg_ortho2,seg_ortho1,M_end2,M_end1,inter,mesh_index,min_dist)){
			std::cout<<"\n PB IN Orthogonal_Projection_M (Projection.cpp): orthogonal projection not determined" << endl;
		}
	}
	proj_points[pointcpp<double>(inter)] = std::make_pair<int,double>(mesh_index,min_dist);

	if (proj_points.size()!=2){
		cout << "\n PB IN Orthogonal_Projection_M (Projection.cpp) : method implemented only for 2 intersections" << endl;
		cout << "Number of intersections is " << proj_points.size();
	}

	return proj_points;
}
