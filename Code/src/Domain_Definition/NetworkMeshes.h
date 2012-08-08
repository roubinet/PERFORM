/*
 * NetworkMeshes.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef NETWORKMESHES_H_
#define NETWORKMESHES_H_

#include "FractureMesh.h"

class NetworkMeshes {
public:
	std::vector<FractureMesh> meshes;	// list of meshes composing the fracture network
public:
	NetworkMeshes();
	virtual ~NetworkMeshes();
	NetworkMeshes(std::string,std::string);
	FractureMesh return_mesh(int);
	FractureMesh operator[](int index){return meshes[index];}
};

#endif /* NETWORKMESHES_H_ */
