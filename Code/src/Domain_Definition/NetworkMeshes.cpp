/*
 * NetworkMeshes.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#include "NetworkMeshes.h"
#include <fstream>

using namespace std;


NetworkMeshes::NetworkMeshes(){}
NetworkMeshes::~NetworkMeshes(){}
NetworkMeshes::NetworkMeshes(string code_path, string file_name_domain){

	// local variables
	string str_line, x_str, y_str;
	int pos1, pos2, mesh_index, nb_segment, index;
	double velocity, aperture;

	// Domain definition from files
	string file_name = code_path + "/Input/Domain_files/";
	file_name += file_name_domain;
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		// number of segments to describe the fracture network
		fichier >> nb_segment;
		meshes = vector<FractureMesh>(nb_segment);
		// each line describes a segment (or mesh) composing the fracture network
		fichier >> str_line;
		for (int i=0;i<nb_segment;i++){
			// definition of the segement's first extremity (coordinate and index)
			pos1 = 0;
			pos1 = str_line.find("(");
			pos2 = str_line.find(",",pos1);
			x_str = str_line.substr(pos1+1,pos2-1);
			pos1 = pos2;
			pos2 = str_line.find(")",pos1);
			y_str = str_line.substr(pos1+1,pos2-pos1-1);
			fichier >> index;
			FluxPoint2D p_origin(atof(x_str.c_str()),atof(y_str.c_str()),index);

			// definition of the segement's second extremity
			fichier >> str_line;
			pos1 = 0;
			pos1 = str_line.find("(");
			pos2 = str_line.find(",",pos1);
			x_str = str_line.substr(pos1+1,pos2-1);
			pos1 = pos2;
			pos2 = str_line.find(")",pos1);
			y_str = str_line.substr(pos1+1,pos2-pos1-1);
			fichier >> index;
			FluxPoint2D p_target(atof(x_str.c_str()),atof(y_str.c_str()),index);

			// flow velocity of the studied segment
			fichier >> velocity;
			// segment aperture
			fichier >> aperture;
			// studied segment numbering
			fichier >> mesh_index;
			// numbering of the neighboring meshes
			vector<int> neigh_meshes;
			fichier >> str_line;
			while ((str_line.find("(")==str_line.npos)&&(!fichier.eof())){
				neigh_meshes.push_back(atoi(str_line.c_str()));
				fichier >> str_line;
			}

			// mesh affectation
			FractureMesh mesh(p_origin,p_target,velocity,aperture,mesh_index,neigh_meshes);
			this->meshes[mesh_index] = mesh;
		}
	}

	else{
		cout << "WARNING in Parameters (Parameters.cpp) : unknown file" << endl;
	}

	fichier.close();
}

FractureMesh NetworkMeshes::return_mesh(int mesh_ind){
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		if (it->mesh_index == mesh_ind){
			return *it;
		}
	}
	cout << "WARNING in return_mesh (NetworkMeshes.cpp): mesh number " << mesh_ind << "not found" << endl;
	return FractureMesh();
}


