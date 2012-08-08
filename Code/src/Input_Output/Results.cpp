/*
 * Results.cpp
 *
 *  Created on: Jun 29, 2012
 *      Author: delphine
 */


#include "Results.h"
#include "../Utilitaries/Structures.h"
#include <fstream>
#include <iostream>
#include <math.h>

using namespace std;



// determine the cumulative distribution of arrival times
Results::Results(map<int,double> arrival_times_,Parameters param){
	arrival_times = arrival_times_;
	Nt = param.Nt;

	// output files name and localization
	string file_name_domain = param.file_name_domain,file_name_param = param.file_name_param;
	int pos = file_name_param.find(".");
	output_file = param.code_path+"/Output/";
	output_file += file_name_param.substr(0,pos);
	output_file += "_";
	pos = file_name_domain.find(".");
	output_file += file_name_domain.substr(0,pos);
	output_file += "_meth";
	stringstream out;
	out << param.simu_option;
	output_file += out.str();
	output_file += ".txt";
}


void Results::post_processing(){
	// time scale of arrival times
	double min_time = get_min_value(arrival_times);
	double max_time = get_max_value(arrival_times);
	vector<double> time = define_scale(min_time,max_time,Nt,SCALE_LOG);
	// cumulative distribution computation
	double current_time; int nb_part;
	for (int i=0;i<Nt;i++){
		current_time = time[i];
		nb_part = 0;
		for (map<int,double>::iterator it=arrival_times.begin();it!=arrival_times.end();it++){
			if (it->second<=current_time) nb_part++;
		}
		cum_dist_times[current_time] = nb_part;
	}
}

void Results::writing(){
	int nb_part = arrival_times.size();
	std::ofstream output(output_file.c_str(), std::ofstream::out);
	for (std::map<double,int>::iterator it=cum_dist_times.begin(); it!=cum_dist_times.end(); it++)
		output << it->first << "	"<< (double)it->second/nb_part << endl;
	output.close();
}



