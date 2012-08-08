/*
 * Transport.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "../Input_Output/Parameters.h"
#include "../Domain_Definition/Domain.h"
#include "../Domain_Definition/NetworkMeshes.h"
#include "../Transport/Particle.h"
#include "../Utilitaries/RandomNumber.h"
#include "../Utilitaries/Constantes.h"
#include "../Utilitaries/Segment.h"
#include "../Utilitaries/Point_Cgal.h"
#include "../Utilitaries/Structures.h"

#include <boost/math/special_functions/erf.hpp>



class PhysicsParam{
public:
	double Dm;	// matrix diffusion coefficient
	double porosity;	// matrix porosity
public:
	PhysicsParam(){};
	virtual ~PhysicsParam(){};
	PhysicsParam(double Dm_,double porosity_){Dm = Dm_; porosity = porosity_;};
};

class NumericalParam{
public:
	int nb_part;	// number of particles for transport simulation
	double proba_transfer;	// determine segment discretization for particle transport
	int simu_option;
public:
	NumericalParam(){};
	virtual ~NumericalParam(){};
	NumericalParam(int nb_part_,double proba_transfer_,int simu_option_){nb_part = nb_part_; proba_transfer = proba_transfer_; simu_option = simu_option_;};
};

class Transport{
public:
	RngStream_a rng_tracker;		///< Random Number Generator for diffusion and dispersion
	NetworkMeshes net_mesh;	// description of the fracture network
	Domain domain;	// geometric description of the domain
	PhysicsParam phys_param;
	NumericalParam num_param;
	Distribution Transfer_Time_Distribution;	///< Transfer time distributions
public:
	Transport(){};
	virtual ~Transport(){};
	Transport(RngStream_a,NetworkMeshes,Domain,Parameters);
	std::vector<Particle> Particles_Injection();
	void infinite_matrix_displacement(Particle&);
	void finite_matrix_displacement(Particle&);
	void particle_displacement(Particle&);
	std::map<int,double> Particles_Transport();
	double Get_Total_Time_From_Advec_Time(double,double,double,double);
	int Mesh_Neighbour_Choice_And_Test(pointcpp<double>,FractureMesh,RngStream_a&);
	bool select_neighbour_slr(std::vector<int>&,std::vector<double>&,FractureMesh,CgalPoint2D);

	// functions for orthogonal projection between fractures
	bool Projection_Periodic_Boundaries(CgalPoint2D,Segment2D,Segment2D,CgalPoint2D,CgalPoint2D,CgalPoint2D&,int&,double&);
	bool Closest_Orthogonal_Projection(CgalPoint2D,Segment2D,CgalPoint2D&,int&,double&);
	Projection_Map Orthogonal_Projection_M(pointcpp<double>,pointcpp<double>);
	void M_project_in_system(CgalPoint2D&,const double&);
	void M_project_in_system(pointcpp<double> &,const double &);

	// functions for transferring to neighboring fractures
	bool Transfer_Probability_and_Transfer_Time(const Projection_Map&,const double,int&,double&);
	bool Transfer_Probability_Computation(const std::pair<double,double> &,const double,double &);
	bool Transfer_Time_Computation(const std::pair<double,double>&,const double,double&);
};

// functions for transferring probability
double Transfer_Probability_FPTD(const double,const double,const double);
bool Transfer_Probability_Feller(const std::pair<double,double> &,const double,const double,double&);
bool Feller_Analytical_Solution(double&,const double,const std::pair<double,double>&,const double);
bool Feller_Analytical_Solution_Single(double&,const double,const double,const double,const double);
double Transfer_Time_FPTD(const double,const double,const double);
bool Transfer_Time_Feller(const std::pair<double,double>&,const double,const double,RngStream_a&,double&,Distribution&);
bool Transfer_Time_Distribution_Computation(Distribution &,const std::pair<double,double> &,const double,const double,RngStream_a&);
bool Feller_Maximal_Time_Determination(double &,const double,const std::pair<double,double> &,const double);
bool Get_Time_Feller(std::pair<std::vector<double>,std::vector<double> >&,const double,double&);
bool Barrier_Choice(const std::map<int,double>&,const double,int&);




#endif /* TRANSPORT_H_ */
