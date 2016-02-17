#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

typedef unsigned int uint;

//data structure of indexed face set
typedef struct {
	//3d coordinates
	float x;
	float y;
	float z;
} Vertex;

typedef struct{
	//three vertex ids
	uint a,b,c;
} Face;

//data structure of Halfedge
typedef struct {
	//A vertex can be expressed in terms of its 3d coordinates
	float x;
	float y;
	float z;

	//It also stores an outgoing half-edge
	HEEdge* outEdge;
} HEVertex;

typedef struct {
	//An edge stores a vertex that it points to
	HEVertex* endVertex;

	//The adjacent face
	HEFace* adjFace;

	//The next and previous half edges on the same adj face
	HEEdge* prevEdge;
	HEEdge* nextEdge;

	//Its twin half edge (points in opposite direction)
	HEEdge* twinEdge;
} HEEdge;

typedef struct {
	//A face stores one of its adjacent Half Edge
	HEEdge* edge;
} HEFace;

class Mesh{
private:
	std::vector<Vertex> V;
	std::vector<Face> F;

	std::vector<HEVertex> HEV;
	std::vector<HEEdge> HEE;
	std::vector<HEFace> HEF;
public:
	Mesh() {};
	Mesh(const char*);
	//load a Mesh from .mesh file
	void loadMF(const char*);
	//write a Mesh to .mesh file (no header)
	void writeMF(const char*);
	//simplify a mesh
	void simplifyMesh(const char* input, const char* output, int faceCnt);
	//turn indexed face set to halfedge
	void convertMesh();
	//turn halfedge to indexed face set
	void revertMesh();
	//helper methods
	std::vector<HEVertex> neighborVertices(HEVertex v);
	std::vector<HEFace> neighborFaces(HEVertex v);
	std::vector<HEVertex> adjacentVertices(HEFace f);
	//return vertex count
	int Vcnt();
	//return face count
	int Fcnt();
};
#endif