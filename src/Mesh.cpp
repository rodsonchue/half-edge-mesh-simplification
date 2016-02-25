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
	std::cout<<"HalfEdge face count: "<< HEF.size()<< std::endl;
	
	std::cout << "hello0" << std::endl; //TODO
	if (numValidFaces > faceCnt) {
		std::cout << "performed collapse" << std::endl; //TODO
		performCollapses(faceCnt);
	}
	std::cout << "hello1" << std::endl;//TODO
	removeInvalids();
	std::cout << "hello2" << std::endl;//TODO

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
		//std::cout << "Working on Vertex (" << vertex.x << ", " << vertex.y << ", " << vertex.z << ")" << std::endl;
		HEVertex* heVertex = (HEVertex*)malloc(sizeof(HEVertex));
		heVertex->x = vertex.x;
		heVertex->y = vertex.y;
		heVertex->z = vertex.z;
		heVertex->isValid = true;

		HEV.push_back(heVertex);
	}

	//Prepare the HEEdge and HEFace
	//And update HEVertex's outgoing edge

	//Stores HEEdges that MAY* have a twin edge, not all edges may have twins (boundary edges)
	std::map<std::pair<uint, uint>, HEEdge*> edgeMap;
	std::map<std::pair<uint, uint>, HEEdge*>::iterator edgeMapIt;
	int nextIndex = 0;

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

		//Shift index by -1 since mesh file starts from 1
		uint va = face.a - 1, vb = face.b - 1, vc = face.c - 1;
		
		//Half edge (c->a)
		HEEdge* halfEdge = (HEEdge*)malloc(sizeof(HEEdge));
		halfEdge->endVertex = HEV[va];
		halfEdge->adjFace = heFace;
		halfEdge->isValid = true;
		halfEdge->index = nextIndex;
		nextIndex++;
		HEE.push_back(halfEdge);
		edgeMap[std::pair<uint,uint>(vc, va)] = halfEdge;

		//Half edge (a->b)
		HEEdge* halfEdgeNext = (HEEdge*)malloc(sizeof(HEEdge));
		halfEdgeNext->endVertex = HEV[vb];
		halfEdgeNext->adjFace = heFace;
		halfEdgeNext->isValid = true;
		halfEdgeNext->index = nextIndex;
		nextIndex++;
		HEE.push_back(halfEdgeNext);
		edgeMap[std::pair<uint, uint>(va, vb)] = halfEdgeNext;

		//Half edge (b->c)
		HEEdge* halfEdgePrev = (HEEdge*)malloc(sizeof(HEEdge));
		halfEdgePrev->endVertex = HEV[vc];
		halfEdgePrev->adjFace = heFace;
		halfEdgePrev->isValid = true;
		halfEdge->index = nextIndex;
		nextIndex++;
		HEE.push_back(halfEdgePrev);
		edgeMap[std::pair<uint, uint>(vb, vc)] = halfEdgePrev;

		//Update HEVertex outgoing edges
		HEV[va]->outEdge = halfEdgeNext;
		HEV[vb]->outEdge = halfEdgePrev;
		HEV[vc]->outEdge = halfEdge;

		//Link the 3 half edges
		halfEdge->nextEdge = halfEdgeNext;
		halfEdge->prevEdge = halfEdgePrev;

		halfEdgeNext->nextEdge = halfEdgePrev;
		halfEdgeNext->prevEdge = halfEdge;

		halfEdgePrev->nextEdge = halfEdge;
		halfEdgePrev->prevEdge = halfEdgeNext;

		//Link twin edges, remove from map once linked (no longer needed to be accessed)
		//Increases access time to the remaining K,V pairs
		edgeMapIt = edgeMap.find(std::pair<uint, uint>(va, vc));
		if (edgeMapIt != edgeMap.end()) {
			edgeMap[std::pair<uint, uint>(va, vc)]->twinEdge = edgeMap[std::pair<uint, uint>(vc, va)];
			edgeMap[std::pair<uint, uint>(vc, va)]->twinEdge = edgeMap[std::pair<uint, uint>(va, vc)];
			edgeMap.erase(std::pair<uint, uint>(va, vc));
			edgeMap.erase(std::pair<uint, uint>(vc, va));
		}
		edgeMapIt = edgeMap.find(std::pair<uint, uint>(vb, va));
		if (edgeMapIt != edgeMap.end()) {
			edgeMap[std::pair<uint, uint>(vb, va)]->twinEdge = edgeMap[std::pair<uint, uint>(va, vb)];
			edgeMap[std::pair<uint, uint>(va, vb)]->twinEdge = edgeMap[std::pair<uint, uint>(vb, va)];
			edgeMap.erase(std::pair<uint, uint>(vb, va));
			edgeMap.erase(std::pair<uint, uint>(va, vb));
		}
		edgeMapIt = edgeMap.find(std::pair<uint, uint>(vc, vb));
		if (edgeMapIt != edgeMap.end()) {
			edgeMap[std::pair<uint, uint>(vc, vb)]->twinEdge = edgeMap[std::pair<uint, uint>(vb, vc)];
			edgeMap[std::pair<uint, uint>(vb, vc)]->twinEdge = edgeMap[std::pair<uint, uint>(vc, vb)];
			edgeMap.erase(std::pair<uint, uint>(vc, vb));
			edgeMap.erase(std::pair<uint, uint>(vb, vc));
		}
		
		//Face stores one of its half edge
		heFace->edge = halfEdge;
		HEF.push_back(heFace);
	}

	//We know there are this many faces
	numValidFaces = HEF.size();

	//Compute collapse cost for each edge
	for (HEEdge* eaEdge : HEE) {
		//TODO
	}
}

void Mesh::revertMesh()
{
	//Empty the current Face and Vertices
	F.clear();
	V.clear();

	//Lookup table to get index given a HEVertex key
	std::map<HEVertex*, uint> vertexIndex;

	for (HEVertex* eaVertex : HEV) {
		/*
		std::cout << "Working on HEVertex (" << eaVertex.x << ", "
											<< eaVertex.y << ", "
											<< eaVertex.z << ", "
											<< ")" << std::endl;*/

		Vertex vertex;
		vertex.x = eaVertex->x;
		vertex.y = eaVertex->y;
		vertex.z = eaVertex->z;
		vertexIndex[eaVertex] = (uint) V.size();
		V.push_back(vertex);
	}

	for (HEFace* eaFace : HEF) {
		/*
		std::cout << "Working on HEFace containing endVertex (" << eaFace.edge->endVertex->x << ", "
																<< eaFace.edge->endVertex->y << ", "
																<< eaFace.edge->endVertex->z << ", "
																<< ")" << std::endl;*/
		Face face;
		face.a = vertexIndex[eaFace->edge->endVertex] + 1;
		face.b = vertexIndex[eaFace->edge->nextEdge->endVertex] + 1;
		face.c = vertexIndex[eaFace->edge->prevEdge->endVertex] + 1;
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
	if (!v->isValid) {
		std::cout << "neighborFaces: invalid vertex" << std::endl; //TODO
		return neighbors;
	}

	HEEdge* firstOutEdge = v->outEdge;
	HEEdge* currOutEdge = v->outEdge;
	bool forceExit = false;

	do {
		std::cout << "neighborFaces: in first half" << std::endl; //TODO
		//The face of this half-edge is a neighbor face
		neighbors.push_back(currOutEdge->adjFace);

		//Get the next outgoing half-edge
		//The next edge of its twin is also an outgoing half-edge from v
		if (currOutEdge->twinEdge != nullptr) {
			currOutEdge = currOutEdge->twinEdge->nextEdge;
		}
		else {
			forceExit = true;
		}

	} while (currOutEdge != firstOutEdge && !forceExit);

	if (forceExit) {
		std::cout << "neighborFaces: in second half" << std::endl; //TODO
		currOutEdge = v->outEdge;
		forceExit = false;
		//Cover the other rotational direction of neighbours
		while (!forceExit){
			//If there is still an adjacent face, continue to traverse
			if (currOutEdge->prevEdge->twinEdge != nullptr) {
				currOutEdge = currOutEdge->prevEdge->twinEdge;

				//That face is also a neighbour face
				neighbors.push_back(currOutEdge->adjFace);
			}
			else {
				//No more neighbours in this direction, exit
				forceExit = true;
			}
		}
	}

	return neighbors;
}

/*
*	A face is bound by 3 vertices, forming a triangle
*/
std::vector<HEVertex*> Mesh::adjacentVertices(HEFace* f)
{
	std::vector<HEVertex*> adjVertices;

	adjVertices.push_back(f->edge->endVertex);
	adjVertices.push_back(f->edge->nextEdge->endVertex);
	adjVertices.push_back(f->edge->prevEdge->endVertex);

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
		std::cout << "u is a lone vertex" << std::endl; //TODO
		return;
	}

	//Otherwise v is related to other half edges and faces
	std::cout << "helloworld" << std::endl; //TODO
	//Update references upon vertex collapse
	std::vector<HEFace*> adjFaces = neighborFaces(u);
	std::cout << "utest1" << std::endl; //TODO
	std::cout << "Number of Faces(prefilter): " << adjFaces.size() << std::endl; //TODO
	std::vector<HEFace*> filteredFaces = filterFaceWithVertex(adjFaces, v);

	std::cout << "Number of Faces: " << filteredFaces.size() << std::endl; //TODO
	if (filteredFaces.size() == 2) {
		std::cout << "Two Filtered Faces" << std::endl; //TODO
		//We begin by looking at the u->v half edge
		HEEdge* uvEdge = getHalfEdge(filteredFaces[0], u, v);
		//std::cout << "test1" << std::endl; //TODO
		if (uvEdge == nullptr) {
			uvEdge = getHalfEdge(filteredFaces[1], u, v);
		}
		//w is the 3rd vertex on that face
		HEEdge* wvEdge = uvEdge->nextEdge->twinEdge;
		HEEdge* uwEdge = uvEdge->prevEdge->twinEdge;
		//std::cout << "test2" << std::endl; //TODO

		//On the twin side,
		HEEdge* vuEdge = uvEdge->twinEdge;
		//x is the 3rd vertex on the other face also containing u,v
		HEEdge* vxEdge = vuEdge->prevEdge->twinEdge;
		HEEdge* xuEdge = vuEdge->nextEdge->twinEdge;
		//std::cout << "test3" << std::endl; //TODO

		//Update twin edge pairings
		makeTwins(uwEdge, wvEdge);
		makeTwins(vxEdge, xuEdge);

		//Invalidate Faces VWU, VUX
		invalidateFace(filteredFaces[0]);
		invalidateFace(filteredFaces[1]);
	}
	else {
		std::cout << "One Filtered Face" << std::endl; //TODO
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
	//std::cout << "Begin replacing u with v" << std::endl; //TODO

	//Update remaining valid faces to use v instead of u
	for (HEFace* eaFace : adjFaces)
	{
		if (isValidFace(eaFace)) {/*
			std::cout << "Working on adj Face with vertices "; //TODO
			for (HEVertex* vert : adjacentVertices(eaFace)) {
				std::cout << debugVertex(vert);
			}
			std::cout << std::endl; //TODO
			*/
			replaceVertex(eaFace, u, v);
		}
	}
	std::cout << "Invalidating vertex "; debugVertex(u); std::cout << std::endl; //TODO

	//Invalidate vertex u
	invalidV.push_back(u);
	u->isValid = false;
	u->outEdge = nullptr;

	//TODO
	//Recompute edge collapse cost in neighbourhood
}

bool Mesh::isValidFace(HEFace* f) {
	return f->edge != nullptr;
}

void Mesh::replaceVertex(HEFace* f, HEVertex* u, HEVertex* v) {
	bool notFound = true;
	HEEdge* initialEdge = f->edge;
	HEEdge* currEdge = f->edge;
	do {/*
		std::cout << "Current Vertex is (" << currEdge->endVertex->x << ", "
			<< currEdge->endVertex->y << ", " 
			<< currEdge->endVertex->z << ")" << std::endl; //TODO*/
		if (currEdge->endVertex == u) {
			//std::cout << "Replace one" << std::endl; //TODO
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
	std::vector<HEVertex*> newHEV;
	for (HEVertex* v : HEV) {
		if (v->isValid) {
			newHEV.push_back(v);
		}
	}
	HEV = newHEV;

	std::vector<HEFace*> newHEF;
	for (HEFace* f : HEF) {
		if (f->edge != nullptr) {
			newHEF.push_back(f);
		}
	}
	HEF = newHEF;

	//*Half edges not neccessarily need to be removed
}

std::vector<HEFace*> Mesh::filterFaceWithVertex(std::vector<HEFace*> vector, HEVertex* v) {
	std::vector<HEFace*> filteredFaces;
	for (HEFace* face : vector) {
		if (faceHasVertex(face, v)) {
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
	if (faces.size() == 2) {
		edge = getHalfEdge(faces[0], u, v);
		edge->twinEdge->isValid = false;
		edge->isValid = false;
	}
	else if (faces.size() == 1) {
		if (!(edge = getHalfEdge(faces[0], u, v))) {
			edge = getHalfEdge(faces[0], v, u);
		}
		edge->isValid = false;
	}

	//otherwise no edge to invalidate
}

std::vector<HEFace*> Mesh::getFacesWithVertices(HEVertex* u, HEVertex* v) {
	return filterFaceWithVertex(neighborFaces(u), v);
}

void Mesh::makeTwins(HEEdge* edge, HEEdge* otherEdge) {
	edge->twinEdge = otherEdge;
	otherEdge->twinEdge = edge;
}

void Mesh::invalidateFace(HEFace* f) {
	//When the face is invalid, so are its adjacent edges
	for (HEEdge* edge : adjacentEdges(f)) {
		std::cout << "Invalidating edge pointing to "; debugVertex(edge->endVertex); std::cout << std::endl; //TODO
		edge->isValid = false;
	}
	std::cout << "Invalidating face " << std::endl; //TODO
	f->edge = nullptr;
	numValidFaces--;
}

void Mesh::performCollapses(int faceCnt) {
	//Perform the required number of edge collapses
	//May over-perform as some collapses remove 1 face, some remove 2
	while(numValidFaces > faceCnt){
		//Version 1:
		//Randomly select one
		int select;
		bool isValidEdge = false;
		while (!isValidEdge){
			select = rand() % HEE.size();
			//std::cout << "Selected value: " << select; //TODO
			if (HEE[select]->isValid) {
				isValidEdge = true;
				//std::cout << " is valid"  << std::endl; //TODO
			}
		}

		HEVertex* vertexU = HEE[select]->prevEdge->endVertex;
		std::cout << "Selected Edge's vertex u is: (" << vertexU->x  << ", "
											<< vertexU->y << ", "
											<< vertexU->z << ")" << std::endl;
		HEVertex* vertexV = HEE[select]->endVertex;
		std::cout << "Selected Edge's vertex v is: (" << vertexV->x << ", "
			<< vertexV->y << ", "
			<< vertexV->z << ")" << std::endl;

		collapseVertex(vertexU, vertexV);
		std::cout << "Current number of faces: " << numValidFaces << std::endl;
	}
}

bool Mesh::faceHasVertex(HEFace* f, HEVertex* v) {
	for (HEEdge* edge : adjacentEdges(f)) {
		if (edge->endVertex == v) {
			return true;
		}
	}
	return false;
}

void Mesh::debugVertex(HEVertex* v) {
	std::cout << "(" << v->x << ", " << v->y << ", " << v->z << ")";;
}