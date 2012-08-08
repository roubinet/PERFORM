/*
 * Transport.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#include "Transport.h"
#include <math.h>

using namespace std;


Transport::Transport(RngStream_a rng_tracker_,NetworkMeshes net_mesh_,Domain domain_,Parameters param){
	rng_tracker = rng_tracker_;
	net_mesh = net_mesh_;
	domain = domain_;
	phys_param = PhysicsParam(param.Dm,param.porosity);
	num_param = NumericalParam(param.nb_part,param.proba_transfer,param.simu_option);
};

// injection of particle on the left border of the domain
vector<Particle> Transport::Particles_Injection(){

	// 0. Local variables
	pointcpp<double> pt1, pt2;
	vector<Particle> part_vect(this->num_param.nb_part);

	// 1. Define positions where particles are injected (on the left border)
	map<int,pointcpp<double> > init_pos;	//<mesh index,position>

	// 2. Search segments (fracture meshes) intersecting the left border
	for (vector<FractureMesh>::iterator it_mesh=net_mesh.meshes.begin();it_mesh!=net_mesh.meshes.end();it_mesh++){
		pt1 = pointcpp<double>(it_mesh->p_ori.p);
		pt2 = pointcpp<double>(it_mesh->p_tar.p);
		if (domain.on_input_limit(pt1)){
			init_pos[it_mesh->mesh_index] = pt1;
		}
		else if (domain.on_input_limit(pt2)){
			init_pos[it_mesh->mesh_index] = pt2;
		}
	}

	// 3. Initial positions affectation
	int nb_part_tot = part_vect.size();
	int nb_pos = init_pos.size();
	int nb_part_pos = (int)(nb_part_tot/(double)nb_pos);	// number of particles per position
	int part_count = 0;	// particle counter
	map<int,pointcpp<double> >::iterator it_pos = init_pos.begin();	// iterator on positions

	for (int i=0;i<nb_part_tot;i++){
		// update particle counter
		if (part_count<=nb_part_pos){part_count++;}
		// update position counter
		else if (part_count>nb_part_pos){
			it_pos++;
			if (it_pos!=init_pos.end()){
				part_count = 1;	// part_count initialization
			}
			else{	// last particles on the last position
				it_pos--;
			}
		}
		part_vect[i].M = it_pos->second;
		part_vect[i].mesh_index = it_pos->first;
		part_vect[i].t = 0;
		part_vect[i].no = i;
	}

	return part_vect;
}

// particle displacement considering the surrounding fractures
void Transport::finite_matrix_displacement(Particle & pa){

	// 0. Variables
	// 0.1. Inputs
	FractureMesh current_mesh = net_mesh.return_mesh(pa.mesh_index);
	pointcpp<double> M_init = pa.M, M_init_proj, M_out;
	double porosity = phys_param.porosity, Dm = phys_param.Dm, transfer_time;
	int transfer_mesh=-1;
	double particle_time;

	// 1. Particle displacement until an intersection
	while (!M_init.identic(M_out, EPSILON)){

		// 1.1. Determination of the characteristic times (depending on fracture discretization)
		// 1.1.1. Determination of the fracture extremity (in the sense of the flow)
		current_mesh.distance_from_M_to_extremity(M_init, M_out);
		// 1.1.2. Determination of the neighboring fractures
		M_init_proj = current_mesh.Mesh_Scale_Advection(M_init, EPSILON);	//small displacement to get the "good" meshes of projection
		Projection_Map projections = Orthogonal_Projection_M(M_init_proj, M_out);
		// 1.1.3. Determination of the corresponding advection and diffusion times (depending on projected distance and transfer_proba_max)
		double t_max = current_mesh.Mesh_Time_Advection(M_init,M_out);
		double t_advec = current_mesh.Advection_Time_Computation(M_init,projections,porosity,Dm,num_param.proba_transfer,t_max);
		double t_diff = current_mesh.Get_Total_Time_From_Advec_Time(t_advec,Dm,porosity,rng_tracker)-t_advec;

		// 1.2. Corresponding displacements
		// 1.2.1. Displacement by advection (and update of the orthogonal projections)
		M_init = current_mesh.Mesh_Scale_Advection(M_init, t_advec);
		if (M_init.identic(M_out, EPSILON)){
			M_init_proj = current_mesh.Mesh_Scale_Advection(M_init, t_advec-EPSILON, true);
		}
		else{M_init_proj = M_init;}
		projections = Orthogonal_Projection_M(M_init_proj, M_out);

		// 1.2.2. Displacement by diffusion
		// if the particle transfers to a neighbor fracture by diffusion
		if (Transfer_Probability_and_Transfer_Time(projections,t_diff,transfer_mesh,transfer_time)){
			particle_time = t_advec+transfer_time;	// particle time updating
			current_mesh = net_mesh.return_mesh(transfer_mesh);	// mesh current updating
			pa.mesh_index = transfer_mesh;
			find_projected_value(projections, transfer_mesh, M_init);	// particle position updating
		}
		// if the particle does not transfer (go back to the initial fracture after diffusion in the matrix)
		else{
			particle_time = t_advec+t_diff;	// particle time updating
		}

		// 1.2.3. Updating of particle and positions characteristics
		pa.t += particle_time;
		pa.M = M_init;
		current_mesh.distance_from_M_to_extremity(M_init, M_out);

	}
}

// particle displacement considering an infinite surrounding matrix
void Transport::infinite_matrix_displacement(Particle & pa){

	// 0. Variables
	FractureMesh current_mesh = net_mesh.return_mesh(pa.mesh_index);
	pointcpp<double> M_init = pa.M, M_out;

	// 1. Displacement from the first to the second extremity of the mesh
	// definition of M_out (extremity opposite to M_init)
	current_mesh.distance_from_M_to_extremity(M_init, M_out);
	double advection_time = current_mesh.Mesh_Time_Advection(M_init, M_out);
	double particle_time = current_mesh.Get_Total_Time_From_Advec_Time(advection_time,this->phys_param.Dm,this->phys_param.porosity,this->rng_tracker);

	// 2. Updating of particle characteristics
	pa.M = M_out;
	pa.t += particle_time;

}

// ADDITION OF NEIGHBOURS EVEN IF THEIR FLUX IS LOWER THAN EPS_POINT_OUT_LIMIT_NO_FLOW TO OBSERVE DIFFUSION IN FRACTURED MEDIA (DR)
// affect the probabilities depending on the possible path with respect to the streamline routine method
bool Transport::select_neighbour_slr(vector<int>& possible_neighbours,vector<double>& possible_probabilities,FractureMesh Mesh_current,CgalPoint2D orip){
	size_t nb_pos = 0;				// Number of possible pathways from "ori" (only segments in which flux leaves "ori")
	double total_flow = 0.;			// Total flow leaving "ori"
	// the possible probabilities are filled with flow values, so compute the total flow and count the number of effectively possible neighbours (flow > 0)
	for(size_t i=0 ; i<possible_probabilities.size() ; i++){
		if(possible_probabilities[i]>0){		//#JRTRACKER+ROMAIN:PRECISION en train de retirer les bras morts (correction prŽalable des flux pour avoir un flux nul sur les artes "bras morts")
			nb_pos++;
			total_flow += possible_probabilities[i];
		}
	}
	// SLR as a function of the number of possible paths
	switch(nb_pos){
		case(0) : { // limit point, the particle must go outside the system
			if(!this->domain.on_output_limit(orip)){
				cout << "WARNINIG in select_neighbour_slr (Transport.cpp): a point with no fluxes leaving the segment edge ori" << endl;
				return false;
			}
			for(size_t i=0 ; i<possible_probabilities.size() ; i++){
				possible_probabilities[i] = 1;
			}
			return true;
				  } // case 0
		case(1) : { // only one possible path, its probability is set to one
			for(size_t i=0 ; i<possible_probabilities.size() ; i++){
				if(possible_probabilities[i]>0)
					possible_probabilities[i] = 1;
				else
					possible_probabilities[i] = 0;
			}
				return true;
				  } // case 1
		case(2) : { // two possibilities : continue in the same fracture or turn (case 1) or turn left or right (case 2)
			// check if we are in case 1 or 2
			int q1 = -1; int q4 = -1;
			for(size_t i=0 ; i<possible_probabilities.size() ; i++){
				CgalVector2D vect1(Mesh_current.p_ori.p, Mesh_current.p_tar.p);
				int index_neigh = possible_neighbours[i];
				FractureMesh mesh_neigh = net_mesh[index_neigh];
				CgalVector2D vect2(orip,mesh_neigh.p_tar.p);
				if((fabs((double)(vect1.x()/vect2.x()-vect1.y()/vect2.y()))<1e-10) && (possible_probabilities[i]>0))
					q1 = i; // index of the continuation in the same fracture
				else if(possible_probabilities[i]>0)
					q4 = i; // index of the bifurcation
			}
			if(q1!=-1 && q4 != -1){ // we find a continuation and a bifurcation, so apply SLR
				if(possible_probabilities[q1]<possible_probabilities[q4]){	// bifurcation dominates, so goes into it
					for(size_t i=0 ; i<possible_probabilities.size() ; i++)
						possible_probabilities[i] = 0.;
					possible_probabilities[q4] = 1;
				}else{														// straight line dominate, mixing probabilities
					double pq4 = possible_probabilities[q4]/possible_probabilities[q1];
					double pq1 = 1.-pq4;
					for(size_t i=0 ; i<possible_probabilities.size() ; i++)
						possible_probabilities[i] = 0.;
					possible_probabilities[q1] = pq1;
					possible_probabilities[q4] = pq4;
				}
			}else{ // we did not find a straight line, so there is two bifurcations, flow-equal probabilities
				for(size_t i=0 ; i<possible_probabilities.size() ; i++){
					if(possible_probabilities[i]>0)
						possible_probabilities[i] /= total_flow;
					else
						possible_probabilities[i] = 0;
				}
			}
				return true;
				  } // case 2
		case(3) : { // all path are possible excepting going back: probabilities = prorata of fluxes
			for(size_t i=0 ; i<possible_probabilities.size() ; i++){
				if(possible_probabilities[i]>0)
					possible_probabilities[i] /= total_flow;
				else
					possible_probabilities[i] = 0;
			}
				return true;
				  } // case 3
		default : { // more than 3 possibilities, wrong case
			cout << "WARNINIG in select_neighbour_slr (Transport.cpp): find more than 3 possible paths, error in the flow computation" << endl;
			return false;
				  }
	}
	return false;
}

// return the index of the next mesh from the current position and the current mesh of the particle
int Transport::Mesh_Neighbour_Choice_And_Test(pointcpp<double> particle_pos, FractureMesh Mesh_current, RngStream_a& rng){

	// 1. Variables
	CgalPoint2D current_pos(particle_pos.i,particle_pos.j);

	// 2. Determination of the corresponding mesh extremity
	FluxPoint2D ori;
	if(CGAL::squared_distance(Mesh_current.p_ori.p,current_pos)<CGAL::squared_distance(Mesh_current.p_tar.p,current_pos))
		ori = Mesh_current.p_ori;
	else
		ori = Mesh_current.p_tar;	// Transforms Flux_point_2D_neighbour in Flux_point_2D

	// 3. Determination of neighboring meshes and the associated probabilities
	vector<int> possible_neighbours;
	vector<double> possible_probabilities;
	vector<int> neighbours = Mesh_current.neigh_meshes;
	// Sorts among the neighbors of "Mesh_current->p_tar" and "Mesh_current->p_ori" the neighbors of "ori"
	size_t neighbours_size=neighbours.size() ;
	for(size_t i=0 ; i<neighbours_size ; i++){
		FractureMesh mesh_neigh = net_mesh[neighbours[i]];
		FluxPoint2D origin = mesh_neigh.p_ori;
		FluxPoint2D target = mesh_neigh.p_tar;
		if (origin.index == ori.index){
			possible_neighbours.push_back(neighbours[i]);
			double flux = net_mesh[neighbours[i]].velocity;		// Gets the flow of the mesh that will be used do determine probabilities
			possible_probabilities.push_back(flux);
		}
		else if (target.index == ori.index){
			possible_neighbours.push_back(neighbours[i]);
			double flux = -net_mesh[neighbours[i]].velocity;		// Gets the flow of the mesh that will be used do determine probabilities
			possible_probabilities.push_back(flux);
		}

	}
	// If no neighbors for the point
	if (possible_neighbours.empty()){
		std::cout << "WARNINIG in Mesh_Neighbour_Choice_And_Test (Transport.cpp): a point is connected to nothing" << endl;
	}

	// 4. Computation of probabilities for the different possible outlet from the inlet point "ori"
	select_neighbour_slr(possible_neighbours,possible_probabilities,Mesh_current,ori.p);

	// 5. Choice of the outlet mesh from the probability distribution
	// Faire une fonction dans RNGstream_a.h
	// cumulative for possible probabilities
	int index = rng.drawPdf(possible_probabilities);

	return possible_neighbours[index];
}


// displacement until the right border of the domain
void Transport::particle_displacement(Particle & pa){

	// 0. Variables
	int simu_option = this->num_param.simu_option;
	int new_mesh_index = -1;

	// 1. Particle displacement until the right boundary of the domain
	while (!domain.on_output_limit(pa.M)){

		// 1.1. Local displacement until a segment extremity
		// 1.1.1. Infinite matrix case
		if (simu_option==INFINITE_MATRIX){
			infinite_matrix_displacement(pa);
		}
		// 1.1.2. Finite matrix case
		else if (simu_option==FINITE_MATRIX){
			finite_matrix_displacement(pa);
		}
		else{cout << "WARNING in particle_displacement (Trasnport.cpp): simulation option not defined" << endl;}

		// 1.2. Choice of the next fracture
		if (!domain.on_output_limit(pa.M)){	// particle did not reach the end of the domain -> determination of the next mesh
			// determination of the next fracture from proportional flow
			FractureMesh mesh_current = net_mesh[pa.mesh_index];
			new_mesh_index = Mesh_Neighbour_Choice_And_Test(pa.M,mesh_current,rng_tracker);
			// particle updating
			pa.mesh_index = new_mesh_index;
		}

	}

}

map<int,double> Transport::Particles_Transport(){

	// 0. Variables
	int nb_part = this->num_param.nb_part;
	map<int,double> arrival_times;	// <particle identifier,particle arrival time>

	// 1. Particles initialization and injection on the left border
	vector<Particle> part_vect = Particles_Injection();

	// 2. Particle displacement
	cout << "particle number = " ;
	for (int i=0;i<nb_part;i++){
		double a = (i+1)/1000.;
		int b = (i+1)/1000;
		if (a==(double)b){cout << "\n" << i+1 << "/" << nb_part;}
		//cout << i << endl;
		Particle pa = part_vect[i];
		particle_displacement(pa);
		arrival_times[pa.no] = pa.t;
	}
	cout << endl;

	return arrival_times;
}
