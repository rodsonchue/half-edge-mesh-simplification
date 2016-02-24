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
		heVertex.isValid = true;

		HEV.push_back(heVertex);
	}

	//Prepare the HEEdge and HEFace
	//And update HEVertex's outgoing edge

	//Stores HEEdges that MAY* have a twin edge, not all edges may have twins (boundary edges)
	std::map<std::pair<uint, uint>, HEEdge*> edgeMap;
	std::map<std::pair<uint, uint>, HEEdge*>::iterator edgeMapIt;

	for (Face face : F)
	{
		//std::cout << "Working on Face with Vertex Indexes (" << face.a << ", " << face.b << ", " << face.c << ")" << std::endl;
		HEFace* heFace = (HEFace*) malloc(sizeof(HEFace));

		/*
		                 a
		                / |\ halfEdge (point of reference)
		halfEdgeNext  |/    \
		              b------>c
		                halfEdgePrev
		*/

		//Each half edge takes their respective vertex index

		//Half edge (c->a)
		HEEdge* halfEdge = (HEEdge*)malloc(sizeof(HEEdge));
		halfEdge->endVertex = &HEV[face.a];
		halfEdge->adjFace = heFace;
		HEE.push_back(*halfEdge);
		edgeMap[std::pair<uint,uint>(face.a, face.c)] = halfEdge;

		//Half edge (a->b)
		HEEdge* halfEdgeNext = (HEEdge*)malloc(sizeof(HEEdge));
		halfEdgeNext->endVertex = &HEV[face.b];
		halfEdgeNext->adjFace = heFace;
		HEE.push_back(*halfEdgeNext);
		edgeMap[std::pair<uint, uint>(face.b, face.a)] = halfEdgeNext;

		//Half edge (b->c)
		HEEdge* halfEdgePrev = (HEEdge*)malloc(sizeof(HEEdge));
		halfEdgePrev->endVertex = &HEV[face.c];
		halfEdgePrev->adjFace = heFace;
		HEE.push_back(*halfEdgePrev);

		//Update HEVertex outgoing edges
		HEV[face.a].outEdge = halfEdgeNext;
		HEV[face.b].outEdge = halfEdgePrev;
		HEV[face.c].outEdge = halfEdge;

		//Link the 3 half edges
		halfEdge->nextEdge = halfEdgeNext;
		halfEdge->prevEdge = halfEdgePrev;

		halfEdgeNext->nextEdge = halfEdgePrev;
		halfEdgeNext->prevEdge = halfEdge;

		halfEdgePrev->nextEdge = halfEdge;
		halfEdgePrev->prevEdge = halfEdgeNext;

		//Link twin edges, remove from map once linked (no longer needed to be accessed)
		//Increases access time to the remaining K,V pairs
		edgeMapIt = edgeMap.find(std::pair<uint, uint>(face.c, face.a));
		if (edgeMapIt != edgeMap.end()) {
			edgeMap[std::pair<uint, uint>(face.a, face.c)]->twinEdge = edgeMapIt->second;
			edgeMap[std::pair<uint, uint>(face.c, face.a)]->twinEdge = halfEdge;
			edgeMap.erase(std::pair<uint, uint>(face.a, face.c));
			edgeMap.erase(std::pair<uint, uint>(face.c, face.a));
		}
		edgeMapIt = edgeMap.find(std::pair<uint, uint>(face.a, face.b));
		if (edgeMapIt != edgeMap.end()) {
			edgeMap[std::pair<uint, uint>(face.b, face.a)]->twinEdge = edgeMapIt->second;
			edgeMap[std::pair<uint, uint>(face.a, face.b)]->twinEdge = halfEdgeNext;
			edgeMap.erase(std::pair<uint, uint>(face.b, face.a));
			edgeMap.erase(std::pair<uint, uint>(face.a, face.b));
		}
		edgeMap[std::pair<uint, uint>(face.c, face.b)] = halfEdgePrev;
		edgeMapIt = edgeMap.find(std::pair<uint, uint>(face.b, face.c));
		if (edgeMapIt != edgeMap.end()) {
			edgeMap[std::pair<uint, uint>(face.b, face.c)]->twinEdge = edgeMapIt->second;
			edgeMap[std::pair<uint, uint>(face.c, face.b)]->twinEdge = halfEdgePrev;
			edgeMap.erase(std::pair<uint, uint>(face.b, face.c));
			edgeMap.erase(std::pair<uint, uint>(face.c, face.b));
		}

		//These Half Edges can now be considered valid as they constitute to form a face
		halfEdge->isValid = true;
		halfEdgeNext->isValid = true;
		halfEdgePrev->isValid = true;
		
		//Face stores one of its half edge
		heFace->edge = halfEdge;
		HEF.push_back(*heFace);
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
		/*
		std::cout << "Working on HEVertex (" << eaVertex.x << ", "
											<< eaVertex.y << ", "
											<< eaVertex.z << ", "
											<< ")" << std::endl;*/

		Vertex vertex;
		vertex.x = eaVertex.x;
		vertex.y = eaVertex.y;
		vertex.z = eaVertex.z;
		vertexIndex[eaVertex] = (uint) V.size();
		V.push_back(vertex);
	}

	for (HEFace eaFace : HEF) {
		/*
		std::cout << "Working on HEFace containing endVertex (" << eaFace.edge->endVertex->x << ", "
																<< eaFace.edge->endVertex->y << ", "
																<< eaFace.edge->endVertex->z << ", "
																<< ")" << std::endl;*/
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
std::vector<HEVertex*> Mesh::neighborVertices(HEVertex* v)
{
	std::vector<HEVertex*> neighbors;
	
	HEEdge* firstOutEdge = v->outEdge;
	HEEdge* currOutEdge = v->outEdge;

	do {
		//The end vertex of this half-edge is a neighbor vertex
		neighbors.push_back(currOutEdge->endVertex);

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
std::vector<HEFace*> Mesh::neighborFaces(HEVertex* v)
{
	std::vector<HEFace*> neighbors;

	HEEdge* firstOutEdge = v->outEdge;
	HEEdge* currOutEdge = v->outEdge;

	do {
		//The end vertex of this half-edge is a neighbor vertex
		neighbors.push_back(currOutEdge->adjFace);

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
std::vector<HEVertex*> Mesh::adjacentVertices(HEFace* f)
{
	std::vector<HEVertex*> adjVertices;

	HEEdge* firstOutEdge = f->edge;
	HEEdge* currOutEdge = f->edge;

	do {
		//The end vertex of this half-edge is a neighbor vertex
		adjVertices.push_back(currOutEdge->endVertex);

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

void Mesh::collapseVertex(HEVertex* u, HEVertex* v) {
	if (!v) {
		//u is a vertex by itself, we just delete it
		delete u;
		return;
	}

	//Otherwise v is related to other half edges and faces

	//Update references upon vertex collapse
	std::vector<HEFace*> adjFaces = neighborFaces(u);
	std::vector<HEFace*> filteredFaces = filterFaceWithVertices(adjFaces, u, v);

	if (filteredFaces.size == 2) {
			//We begin by looking at the u->v half edge
			HEEdge* uvEdge = getHalfEdge(filteredFaces[0], u, v);
			//w is the 3rd vertex on that face
			HEEdge* wvEdge = uvEdge->nextEdge->twinEdge;
			HEEdge* uwEdge = uvEdge->prevEdge->twinEdge;
			//On the twin side,
			HEEdge* vuEdge = uvEdge->twinEdge;
			//x is the 3rd vertex on the other face also containing u,v
			HEEdge* vxEdge = vuEdge->prevEdge->twinEdge;
			HEEdge* xuEdge = vuEdge->nextEdge->twinEdge;

			//Update twin edge pairings
			makeTwins(uwEdge, wvEdge);
			makeTwins(vxEdge, xuEdge);

			//Invalidate Faces VWU, VUX
			invalidateFace(filteredFaces[0]);
			invalidateFace(filteredFaces[1]);
	}
	else {
		//only 1 face, is corner case
		//We begin by looking at the u->v half edge
		HEEdge* uvEdge = getHalfEdge(filteredFaces[0], u, v);

		//Non null means is on the correct side
		if (uvEdge) {
			//w is the 3rd vertex on that face
			HEEdge* wvEdge = uvEdge->nextEdge->twinEdge;
			HEEdge* uwEdge = uvEdge->prevEdge->twinEdge;

			//Update twin edge pairings
			makeTwins(uwEdge, wvEdge);

			//Invalidate Face
			invalidateFace(filteredFaces[0]);
		}
		else { //need to look at the other half
			HEEdge* vuEdge = getHalfEdge(filteredFaces[0], v, u);

			//x is the 3rd vertex on the other face also containing u,v
			HEEdge* vxEdge = vuEdge->prevEdge->twinEdge;
			HEEdge* xuEdge = vuEdge->nextEdge->twinEdge;

			//Update twin edge pairings
			makeTwins(vxEdge, xuEdge);

			//Invalidate Face
			invalidateFace(filteredFaces[0]);
		}
	}

	//Update remaining valid faces to use v instead of u
	for (HEFace* eaFace : adjFaces)
	{
		if (!eaFace->edge) {
			replaceVertex(eaFace, u, v);
		}
	}

	//Invalidate vertex u
	u->isValid = false;

	//TODO
	//Recompute edge collapse cost in neighbourhood
}

void Mesh::replaceVertex(HEFace* f, HEVertex* u, HEVertex* v) {
	bool notFound = true;
	HEEdge* initialEdge = f->edge;
	HEEdge* currEdge = f->edge;
	do {
		if (currEdge->endVertex == u) {
			currEdge->endVertex = v;
			return;
		}
		currEdge = currEdge->nextEdge;
	} while (currEdge != initialEdge);
	
	//Only reaches here if the vertex u is not part of the face f
	return;
}

std::vector<HEEdge*> Mesh::adjacentEdges(HEFace* f) {
	std::vector<HEEdge*> adjEdges;
	if (!f->edge) {
		return adjEdges;
	}

	HEEdge* initialEdge = f->edge;
	HEEdge* currEdge = f->edge;
	do {
		adjEdges.push_back(currEdge);
		currEdge = currEdge->nextEdge;
	} while (currEdge != initialEdge);

	return adjEdges;
}

void Mesh::removeInvalids() {
	std::vector<HEVertex> newHEV;
	for (HEVertex v : HEV) {
		if (v.isValid) {
			newHEV.push_back(v);
		}
	}
	HEV = newHEV;

	std::vector<HEFace> newHEF;
	for (HEFace f : HEF) {
		if (f.edge != nullptr) {
			newHEF.push_back(f);
		}
	}
	HEF = newHEF;

	//*Half edges not neccessarily need to be removed
}

std::vector<HEFace*> Mesh::filterFaceWithVertices(std::vector<HEFace*> vector,
	HEVertex* u, HEVertex* v) {
	std::vector<HEFace*> filteredFaces;
	for (HEFace* face : vector) {
		bool hasU = false, hasV = false;
		for (HEVertex* vertex : adjacentVertices(face)) {
			if (vertex == u) {
				hasU = true;
			}
			if (vertex == v) {
				hasV = true;
			}
		}

		if (hasU && hasV) {
			filteredFaces.push_back(face);
		}
	}
	return filteredFaces;
}

HEEdge* Mesh::getHalfEdge(HEFace* face, HEVertex* u, HEVertex* v) {
	HEEdge* edge = face->edge;

	if (edge->endVertex == v && edge->prevEdge->endVertex == u) {
		return edge;
	}
	else if (edge->prevEdge->endVertex == v && edge->nextEdge->endVertex == u) {
		return edge->prevEdge;
	}
	else if (edge->nextEdge->endVertex == v && edge->endVertex == u) {
		return edge->nextEdge;
	}

	return nullptr;
}

void Mesh::invalidateEdgePairs(HEVertex* u, HEVertex* v) {
	//Get both half edges
	std::vector<HEFace*> faces = getFacesWithVertices(u, v);
	//There are either two faces, one, or none (vertices are not directly connected)
	HEEdge* edge;
	if (faces.size == 2) {
		edge = getHalfEdge(faces[0], u, v);
		edge->twinEdge->isValid = false;
		edge->isValid = false;
	}
	else if (faces.size == 1) {
		if (!(edge = getHalfEdge(faces[0], u, v))) {
			edge = getHalfEdge(faces[0], v, u);
		}
		edge->isValid = false;
	}

	//otherwise no edge to invalidate
}

std::vector<HEFace*> Mesh::getFacesWithVertices(HEVertex* u, HEVertex* v) {
	return filterFaceWithVertices(neighborFaces(u), u, v);
}

void Mesh::makeTwins(HEEdge* edge, HEEdge* otherEdge) {
	edge->twinEdge = otherEdge;
	otherEdge->twinEdge = edge;
}

void Mesh::invalidateFace(HEFace* f) {
	//When the face is invalid, so are its adjacent edges
	for (HEEdge* edge : adjacentEdges(f)) {
		edge->isValid = false;
	}

	f->edge = nullptr;
}