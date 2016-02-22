#include "stdafx.h"

/*
Author: Rodson Chue Le Sheng
Matric No. A0110787A

Note: Code is based on a skeleton template provided as part of a module assignment

Mesh reduction by edge collapse based on Melax's Edge Collapse algorithm @see mesh.h
*/

int main(int argc, char* argv[]){
	if(argc!=4) std::cout<<"Usage: ./exe input output faceCnt\n";
	else{
		Mesh mesh;
		mesh.simplifyMesh(argv[1],argv[2],atoi(argv[3]));
	}
	return 0;
}