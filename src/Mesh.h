/*	CS3242 3D Modeling And Animation
*	Programming Assignment I
*	School of Computing
*	National University of Singapore
*/
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
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h> //for rng

typedef unsigned int uint;

//data structure of indexed face set
typedef struct Vertex {
	//3d coordinates
	float x;
	float y;
	float z;
} Vertex;

//A vertex can also be used as 3 dimensional vector
typedef Vertex Vector3;

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

	//Whether edge collapse is dangerous (non-manifold causing)
	bool isDanger = false;

	//Index to identify each unique half edge
	int index;

	bool operator==(HEEdge other) const {
		return index == other.index;
	}

	bool operator<(HEEdge other) const {
		//If both are invalid or valid
		if (isValid == other.isValid) {
			if (isDanger == other.isDanger) {
				return cost < other.cost;
			}
			else {
				//this is not dangerous but the other is
				return !isDanger;
			}
		} else {
			//this is valid and the other is not.
			return isValid;
		}
	}

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

	std::vector<HEVertex*> HEV;
	std::vector<HEEdge*> HEE;
	std::vector<HEFace*> HEF;

	//Quick way to check number of valid faces
	int numValidFaces = 0;

	//Internal functions
	void replaceVertex(HEFace* f, HEVertex* u, HEVertex* v);
	std::vector<HEFace*> filterFaceWithVertex(std::vector<HEFace*> vector, HEVertex* v);
	HEEdge* getHalfEdge(HEFace* face, HEVertex* u, HEVertex* v);
	void invalidateEdgePairs(HEVertex* u, HEVertex* v);
	std::vector<HEFace*> getFacesWithVertices(HEVertex* u, HEVertex* v);
	void makeTwins(HEEdge* edge, HEEdge* otherEdge);
	void invalidateFace(HEFace* f);
	bool isValidFace(HEFace* f);
	bool faceHasVertex(HEFace* f, HEVertex* v);
	void debugVertex(HEVertex* v); //Remove after TODO
	int getNumVerticesAdjacentTo(HEVertex* u, HEVertex* v);
	std::vector<HEEdge*> getIncomingEdges(HEVertex* v);
	std::vector<HEEdge*> getOutgoingEdges(HEVertex* v);
	std::vector<HEEdge*> getAllNeighbourhoodEdges(HEVertex* v);
	std::set<HEEdge*> getTwoRingNeighbourhoodEdges(HEVertex* v);
	void computeCollapseCost(HEEdge* e);
	int selectLeastCostEdge();
	float magnitude(HEVertex* u, HEVertex* v);
	float computeDotProduct(HEFace* faceA, HEFace* faceB);
	float computeDotProduct(Vector3 vA, Vector3 vB);
	Vector3 computeNormalizedNormal(HEFace* face);
	Vector3 computeNormalizedNormal(Vector3 a, Vector3 b, Vector3 c);
	Vector3 assignOtherIfNotEqual(HEVertex* current, HEVertex* comparedV, HEVertex* otherV);
	bool checkGeometry(HEEdge* edge);
	bool faceFlips(HEFace* f, HEVertex* u, HEVertex* v);
	bool canCollapse(HEEdge* __edge);
public:
	Mesh() {};
	Mesh(const char*);
	//load a Mesh from .mesh file
	void loadMF(const char*);
	//write a Mesh to .mesh file (no header)
	void writeMF(const char*);
	//simplify a mesh
	void simplifyMesh(const char* input, const char* output, int faceCnt);
	//Performs edge collapses depending on required face count
	void performCollapses(int faceCnt);
	//turn indexed face set to halfedge
	void convertMesh();
	//turn halfedge to indexed face set
	void revertMesh();
	//Collapses an edge (if it is valid)
	void collapseEdge(HEEdge* e);
	//helper methods
	std::vector<HEVertex*> neighborVertices(HEVertex* v);
	std::vector<HEEdge*> outgoingEdges(HEVertex* v);
	std::vector<HEFace*> neighborFaces(HEVertex* v);
	std::vector<HEVertex*> adjacentVertices(HEFace* f);
	std::vector<HEEdge*> adjacentEdges(HEFace* f);
	std::pair<HEVertex*, HEVertex*> getVertices(HEEdge* e);
	void removeInvalids();
	//return vertex count
	int Vcnt();
	//return face count
	int Fcnt();
};
#endif