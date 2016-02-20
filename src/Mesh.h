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
	struct HEEdge* outEdge;

	bool operator==(HEVertex other) const
	{
		return x == other.x
			&& y == other.y
			&& z == other.z;
	}

	//Used when using struct as key of map
	bool operator<(HEVertex other) const
	{
		return x < other.x
			|| (x == other.x && y < other.y)
			|| (x == other.x && y == other.y && z < other.z);
	}

} HEVertex;

//Edge of Half-Edge Data Structure
typedef struct HEEdge {
	//An edge stores a vertex that it points to
	struct HEVertex* endVertex;

	//The adjacent face
	struct HEFace* adjFace;

	//The next and previous half edges on the same adj face
	struct HEEdge* prevEdge;
	struct HEEdge* nextEdge;

	//Its twin half edge (points in opposite direction)
	struct HEEdge* twinEdge;

	bool operator!=(HEEdge other) const
	{
		return !(endVertex == other.endVertex);
	}

	bool operator==(HEEdge other) const
	{
		return endVertex == other.endVertex;
	}

} HEEdge;

//Face of Half-Edge Data Structure
typedef struct HEFace {
	//A face stores one of its adjacent Half Edge
	struct HEEdge* edge;
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