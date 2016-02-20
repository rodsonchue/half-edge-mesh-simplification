#include "stdafx.h"
#include "Mesh.h"
#include <map>

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
	std::cout << "Original face count: " << F.size() << std::endl;
	std::cout<<"HalfEdge face count: "<<HEF.size()<<std::endl;
	// do mesh simplification here

	revertMesh();
	writeMF(output);
}

void Mesh::convertMesh()
{
	if (HEV.size()>0) HEV.clear();
	if (HEE.size()>0) HEE.clear();
	if (HEF.size()>0) HEF.clear();

	//Prepare the HEVertex
	for (Vertex vertex : V)
	{
		std::cout << "Working on Vertex (" << vertex.x << ", " << vertex.y << ", " << vertex.z << ")" << std::endl;
		HEVertex heVertex;
		heVertex.x = vertex.x;
		heVertex.y = vertex.y;
		heVertex.z = vertex.z;

		HEV.push_back(heVertex);
	}

	//Prepare the HEEdge and HEFace
	//And update HEVertex's outgoing edge
	std::map<std::pair<uint, uint>, HEEdge*> edgeMap;
	std::map<std::pair<uint, uint>, HEEdge*>::iterator edgeMapIt;

	for (Face face : F)
	{
		std::cout << "Working on Face with Vertex Indexes (" << face.a << ", " << face.b << ", " << face.c << ")" << std::endl;
		HEFace heFace;

		/*
		                 a
		                / |\ halfEdge
		halfEdgeNext  |/    \
		              b------>c
		                halfEdgePrev
		*/

		//Each half edge takes their respective vertex index
		HEEdge halfEdge;
		halfEdge.endVertex = &HEV[face.a];
		halfEdge.adjFace = &heFace;
		HEE.push_back(halfEdge);
		edgeMap[std::pair<uint,uint>(face.a, face.c)] = &halfEdge;

		edgeMapIt = edgeMap.find(std::pair<uint, uint>(face.c, face.a));
		if (edgeMapIt != edgeMap.end()) {
			edgeMap[std::pair<uint, uint>(face.a, face.c)]->twinEdge = edgeMapIt->second;
			edgeMap[std::pair<uint, uint>(face.c, face.a)]->twinEdge = &halfEdge;
		}

		HEEdge halfEdgeNext;
		halfEdgeNext.endVertex = &HEV[face.b];
		halfEdgeNext.adjFace = &heFace;
		HEE.push_back(halfEdgeNext);
		edgeMap[std::pair<uint, uint>(face.b, face.a)] = &halfEdgeNext;

		edgeMapIt = edgeMap.find(std::pair<uint, uint>(face.a, face.b));
		if (edgeMapIt != edgeMap.end()) {
			edgeMap[std::pair<uint, uint>(face.b, face.a)]->twinEdge = edgeMapIt->second;
			edgeMap[std::pair<uint, uint>(face.a, face.b)]->twinEdge = &halfEdgeNext;
		}

		HEEdge halfEdgePrev;
		halfEdgePrev.endVertex = &HEV[face.c];
		halfEdgePrev.adjFace = &heFace;
		HEE.push_back(halfEdgePrev);
		edgeMap[std::pair<uint, uint>(face.c, face.b)] = &halfEdgePrev;

		edgeMapIt = edgeMap.find(std::pair<uint, uint>(face.b, face.c));
		if (edgeMapIt != edgeMap.end()) {
			edgeMap[std::pair<uint, uint>(face.b, face.c)]->twinEdge = edgeMapIt->second;
			edgeMap[std::pair<uint, uint>(face.c, face.b)]->twinEdge = &halfEdgePrev;
		}

		//Update HEVertex outgoing edges
		HEV[face.a].outEdge = &halfEdgeNext;
		HEV[face.b].outEdge = &halfEdgePrev;
		HEV[face.c].outEdge = &halfEdge;

		//Face stores one of its half edge
		heFace.edge = &halfEdge;
		HEF.push_back(heFace);

		//Link the 3 half edges
		halfEdge.nextEdge = &halfEdgeNext;
		halfEdge.prevEdge = &halfEdgePrev;

		halfEdgeNext.nextEdge = &halfEdgePrev;
		halfEdgeNext.prevEdge = &halfEdge;

		halfEdgePrev.nextEdge = &halfEdge;
		halfEdgePrev.prevEdge = &halfEdgeNext;
	}
}

void Mesh::revertMesh()
{
	//Empty the current Face and Vertices
	F.clear();
	V.clear();

	//Lookup table to get index given a HEVertex key
	std::map<HEVertex, uint> vertexIndex;

	for (HEVertex eaVertex : HEV) {
		std::cout << "Working on HEVertex (" << eaVertex.x << ", "
											<< eaVertex.y << ", "
											<< eaVertex.z << ", "
											<< ")" << std::endl;

		Vertex vertex;
		vertex.x = eaVertex.x;
		vertex.y = eaVertex.y;
		vertex.z = eaVertex.z;
		vertexIndex[eaVertex] = (uint) V.size();
		V.push_back(vertex);
	}

	for (HEFace eaFace : HEF) {
		
		std::cout << "Working on HEFace containing endVertex (" << eaFace.edge->endVertex->x << ", "
																<< eaFace.edge->endVertex->y << ", "
																<< eaFace.edge->endVertex->z << ", "
																<< ")" << std::endl;
		Face face;
		face.a = vertexIndex[*eaFace.edge->endVertex];
		face.b = vertexIndex[*eaFace.edge->nextEdge->endVertex];
		face.c = vertexIndex[*eaFace.edge->prevEdge->endVertex];
		F.push_back(face);
	}
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
		adjVertices.push_back(*currOutEdge->endVertex);

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