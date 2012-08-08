/*
 * Parameters.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#include "Parameters.h"
#include <fstream>
#include <iostream>

using namespace std;

Parameters::Parameters(){
	cout << "Enter the path for Input and Output folders" << endl;
	cin >> code_path;
	string file_name = code_path+"/Input/File_names.txt";
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		fichier >> file_name_param;
		fichier >> file_name_domain;
	}
	else{cout << "Parameters files not found" << endl;}
	read_param();
}


// Parameter reading and affectation
void Parameters::read_param(){
	string file_name = code_path+"/Input/Param_files/";
	file_name += file_name_param;
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		fichier >> this->Lx;
		fichier >> this->Ly;
		fichier >> this->Dm;
		fichier >> this->porosity;
		fichier >> this->nb_part;
		fichier >> this->proba_transfer;
		fichier >> this->simu_option;
		fichier >> this->Nt;
	}
	else{
		cout << "WARNING in Parameters (Parameters.cpp) : unknown file" << endl;
	}
	fichier.close();
}
