#include "Mesh.h"

Mesh::Mesh(const char* filename){	
	loadMF(filename);
}

void Mesh::loadMF(const char* filename){
	if(V.size()>0) V.clear();
	if(F.size()>0) F.clear();
	std::ifstream infile;
	infile.open(filename, std::ios::in);
	std::string strbuff;
	while(std::getline(infile,strbuff)){
		std::stringstream ss;
		ss<<strbuff;
		char type;
		ss>>type;
		if(type=='v'){
			Vertex v;
			ss>>v.x>>v.y>>v.z;
			V.push_back(v);
		}
		else if(type=='f'){
			Face f;
			ss>>f.a>>f.b>>f.c;
			F.push_back(f);
		}
	}
	infile.close();
}

void Mesh::writeMF(const char* filename){
	std::ofstream outfile;
	outfile.open(filename, std::ios::out);
	std::string strbuff;
	for(uint i=0;i<V.size();i++){
		outfile<<"v "<<V[i].x<<" "<<V[i].y<<" "<<V[i].z<<std::endl;
	}
	for(uint i=0;i<F.size();i++){
		outfile<<"f "<<F[i].a<<" "<<F[i].b<<" "<<F[i].c<<std::endl;
	}
	outfile.close();
}

void Mesh::simplifyMesh(const char* input, const char* output, int faceCnt){
	// you may assume inputs are always valid
	loadMF(input);
	convertMesh();
	std::cout<<"Original face count: "<<HEF.size()<<std::endl;
	// do mesh simplification here

	revertMesh();
	writeMF();
}

void Mesh::convertMesh()
{
	//TODO
	for (Face face : F)
	{
		//Each 3-verticed face has 3 Half-edges
		HEEdge halfEdge;
		halfEdge.endVertex;

	}
}

void Mesh::revertMesh()
{
	//TODO
}

bool operator==(const HEVertex& lhs, const HEVertex& rhs)
{
	return lhs.x == rhs.x
		&& lhs.y == rhs.y
		&& lhs.z == rhs.z;
}

bool operator!=(const HEEdge& lhs, const HEEdge& rhs)
{
	return !(lhs.endVertex == rhs.endVertex);
}

bool operator==(const HEEdge& lhs, const HEEdge& rhs)
{
	return lhs.endVertex == rhs.endVertex;
}

/*
*	General approach: get all neighbour vertices by getting all outgoing half-edges from v
*	See also: HEEdge equality using != operator (minimalistic)
*/
std::vector<HEVertex> Mesh::neighborVertices(HEVertex v)
{
	std::vector<HEVertex> neighbors;
	
	HEEdge* firstOutEdge = v.outEdge;
	HEEdge* currOutEdge = v.outEdge;

	do {
		//The end vertex of this half-edge is a neighbor vertex
		neighbors.push_back(*currOutEdge->endVertex);

		//Get the next outgoing half-edge
		//The next edge of its twin is also an outgoing half-edge from v
		currOutEdge = currOutEdge->twinEdge->nextEdge;

	} while (currOutEdge != firstOutEdge);

	return neighbors;
}

/*
*	General approach: get all neighbour faces by getting all outgoing half-edges from v
*	See also: HEEdge equality using != operator (minimalistic)
*/
std::vector<HEFace> Mesh::neighborFaces(HEVertex v)
{
	std::vector<HEFace> neighbors;

	HEEdge* firstOutEdge = v.outEdge;
	HEEdge* currOutEdge = v.outEdge;

	do {
		//The end vertex of this half-edge is a neighbor vertex
		neighbors.push_back(*currOutEdge->adjFace);

		//Get the next outgoing half-edge
		//The next edge of its twin is also an outgoing half-edge from v
		currOutEdge = currOutEdge->twinEdge->nextEdge;

	} while (currOutEdge != firstOutEdge);

	return neighbors;
}

/*
*	General approach: get all adjacent vertices by getting all half-edges adjacent to face
*	See also: HEEdge equality using != operator (minimalistic)
*/
std::vector<HEVertex> Mesh::adjacentVertices(HEFace f)
{
	std::vector<HEVertex> adjVertices;

	HEEdge* firstOutEdge = f.edge;
	HEEdge* currOutEdge = f.edge;

	do {
		//The end vertex of this half-edge is a neighbor vertex
		adjVertices.push_back(*currOutEdge->adjFace);

		//Get the next outgoing half-edge
		//The next edge is also adjacent to the same face
		currOutEdge = currOutEdge->nextEdge;

	} while (currOutEdge != firstOutEdge);

	return adjVertices;
}

int Mesh::Vcnt(){
	return HEV.size();
}

int Mesh::Fcnt(){
	return HEF.size();
}