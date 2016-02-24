/*
Author: Rodson Chue Le Sheng
Matric No. A0110787A
Last updated: 24/02/2016

Note: Code is based on a skeleton template provided as part of a module assignment

Mesh reduction by Edge Collapse based on Melax's Selection criteria

A publically accessible document regarding this approach is accessible here: 
http://dev.gameres.com/program/Visual/3D/PolygonReduction.pdf

Differences: Using half-edge data structure instead of the suggested "triangle" structure
*/
#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

typedef unsigned int uint;

//data structure of indexed face set
typedef struct Vertex {
	//3d coordinates
	float x;
	float y;
	float z;
} Vertex;

typedef struct Face {
	//three vertex ids
	uint a,b,c;
} Face;

//data structure of Halfedge
//Vertex of Half-Edge Data Structure
typedef struct HEVertex {
	//A vertex can be expressed in terms of its 3d coordinates
	float x;
	float y;
	float z;

	//It also stores an outgoing half-edge
	struct HEEdge* outEdge = nullptr;

	//Validity of vertex, vertex without an outgoing half-edge may or may not be valid
	bool isValid = false;

	bool operator==(HEVertex other) const {
		return x == other.x
			&& y == other.y
			&& z == other.z;
	}

	//Used when using struct as key of map, is the default comparator
	bool operator<(HEVertex other) const {
		return x < other.x
			|| (x == other.x && y < other.y)
			|| (x == other.x && y == other.y && z < other.z);
	}

} HEVertex;

//Edge of Half-Edge Data Structure
typedef struct HEEdge {
	//States validity of half edge
	bool isValid = false;

	//An edge stores a vertex that it points to
	struct HEVertex* endVertex = nullptr;

	//The adjacent face
	struct HEFace* adjFace = nullptr;

	//The next and previous half edges on the same adj face
	struct HEEdge* prevEdge = nullptr;
	struct HEEdge* nextEdge = nullptr;

	//Its twin half edge (points in opposite direction)
	struct HEEdge* twinEdge = nullptr;

	//The cost of collapsing edge
	float cost = FLT_MAX;
} HEEdge;

//Face of Half-Edge Data Structure
typedef struct HEFace {
	//A face stores one of its adjacent Half Edge
	struct HEEdge* edge = nullptr;

} HEFace;

class Mesh{
private:
	std::vector<Vertex> V;
	std::vector<Face> F;

	std::vector<HEVertex> HEV;
	std::vector<HEEdge> HEE;
	std::vector<HEFace> HEF;

	//Internal functions
	void replaceVertex(HEFace* f, HEVertex* u, HEVertex* v);
	std::vector<HEFace*> filterFaceWithVertices(std::vector<HEFace*> vector, HEVertex* u, HEVertex* v);
	HEEdge* getHalfEdge(HEFace* face, HEVertex* u, HEVertex* v);
	void invalidateEdgePairs(HEVertex* u, HEVertex* v);
	std::vector<HEFace*> getFacesWithVertices(HEVertex* u, HEVertex* v);
	void makeTwins(HEEdge* edge, HEEdge* otherEdge);
	void invalidateFace(HEFace* f);
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
	//Collapses a vertex to another vertex (or nothing)
	void collapseVertex(HEVertex* u, HEVertex* v);
	//helper methods
	std::vector<HEVertex*> neighborVertices(HEVertex* v);
	std::vector<HEFace*> neighborFaces(HEVertex* v);
	std::vector<HEVertex*> adjacentVertices(HEFace* f);
	std::vector<HEEdge*> adjacentEdges(HEFace* f);
	void removeInvalids();
	//return vertex count
	int Vcnt();
	//return face count
	int Fcnt();
};
#endif