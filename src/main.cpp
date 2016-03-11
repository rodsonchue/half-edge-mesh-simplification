/*	CS3242 3D Modeling And Animation
*	Programming Assignment I
*	School of Computing
*	National University of Singapore
*/
#include "stdafx.h"

/*
Author: Rodson Chue Le Sheng
Matric No. A0110787A
Last updated: 24/02/2016

Note: Code is based on a skeleton template provided as part of a module assignment

Mesh reduction by Edge Collapse based on Melax's Selection criteria

A publically accessible document regarding this approach is accessible here:
http://dev.gameres.com/program/Visual/3D/PolygonReduction.pdf

Majority of implementation is done in the Mesh class. (Mesh.h + Mesh.cpp)
*/

int main(int argc, char* argv[]){
	if(argc!=4) std::cout<<"Usage: ./exe input output faceCnt\n";
	else{
		Mesh mesh;
		mesh.simplifyMesh(argv[1],argv[2],atoi(argv[3]));
	}
	return 0;
}