#include "stdafx.h"
#include "Mesh.h"
#include <map>
#include <algorithm>
#include <math.h>
#include <assert.h>

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
	std::cout << "after collapses" << std::endl;//TODO
	removeInvalids();
	std::cout << "after removing invalids" << std::endl;//TODO

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

	//Check whether mesh is "closed"
	if (!edgeMap.empty()) {
		std::cout << "there are edges with no twins!" << std::endl;
	}

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
	std::cout << "in neighbor vertices" << std::endl; //TODO
	std::vector<HEVertex*> neighbors;
	HEEdge* firstOutEdge = v->outEdge;
	HEEdge* currOutEdge = v->outEdge;

	if (!currOutEdge->isValid) {
		std::cout << "invalid edge" << std::endl; //TODO
		return neighbors;
	}
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
	if (v==nullptr || !v->isValid) {
		std::cout << "neighborFaces: invalid vertex" << std::endl; //TODO
		return neighbors;
	}

	HEEdge* firstOutEdge = v->outEdge;
	HEEdge* currOutEdge = v->outEdge;
	bool forceExit = false;

	do {
		//std::cout << "neighborFaces: in first half" << std::endl; //TODO
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

void Mesh::collapseEdge(HEEdge* _edge) {

	//Edges on one side of the face
	HEEdge* e = _edge;
	HEEdge* en = _edge->nextEdge;
	HEEdge* ep = _edge->prevEdge;

	//Edges on the opposite face
	HEEdge* o = _edge->twinEdge;
	HEEdge* on = _edge->twinEdge->nextEdge;
	HEEdge* op = _edge->twinEdge->prevEdge;

	//The faces
	HEFace* fe = e->adjFace;
	HEFace* fo = o->adjFace;

	//The vertices
	HEVertex* ve = e->endVertex;
	HEVertex* vo = o->endVertex;

	//All edges that once ended at vo will now end at ve
	//(vo is being collapsed to ve)
	std::cout << "before incoming edges" << std::endl; //TODO
	for (HEEdge* eaIncomingEdge : getIncomingEdges(vo)) {
		eaIncomingEdge->endVertex = ve;
	}std::cout << "after incoming edges" << std::endl; //TODO
	
	//Reassign ve with any valid outEdge if needed
	if (ve->outEdge == o || ve->outEdge == en) {
		ve->outEdge = op->twinEdge;
		std::cout << "Reassigning ve outedge" << std::endl; //TODO
	}
	assert(ve->outEdge->isValid);

	//Since we are dealing with triangles, 
	//The opposite halfedges of the now-invalid half edges now need to be "twinned"
	makeTwins(en->twinEdge, ep->twinEdge);
	makeTwins(on->twinEdge, op->twinEdge);
	std::cout << "after twins" << std::endl; //TODO

	//Faces fe,fo also collapse (with their adjacent edges)
	invalidateFace(fe);
	invalidateFace(fo);
	std::cout << "after invalidate face" << std::endl; //TODO

	//Invalidate the collapsed vertex
	vo->isValid = false;

	std::cout << "ve's out edge is " << ve->outEdge->isValid << std::endl; //TODO

	//Recompute edge collapse cost in neighbourhood of ve
	for (HEVertex* v : neighborVertices(ve)) {
		for (HEEdge* eaEdge : outgoingEdges(v)) {
			computeCollapseCost(eaEdge);
		}
	}
	//And for itself as well
	for (HEEdge* eaEdge : outgoingEdges(ve)) {
		computeCollapseCost(eaEdge);
	}
	
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
	if (edge != nullptr) {
		edge->twinEdge = otherEdge;
	}
	
	if (otherEdge != nullptr) {
		otherEdge->twinEdge = edge;
	}
}

void Mesh::invalidateFace(HEFace* f) {
	//When the face is invalid, so are its adjacent edges
	for (HEEdge* edge : adjacentEdges(f)) {
		//std::cout << "Invalidating edge pointing to "; debugVertex(edge->endVertex); std::cout << std::endl; //TODO
		edge->isValid = false;
	}
	//std::cout << "Invalidating face " << std::endl; //TODO
	f->edge = nullptr;
	numValidFaces--;
}

void Mesh::performCollapses(int faceCnt) {
	//First compute all collapse costs
	//Not needed if doing random-select collapse
	for (HEEdge* eaEdge : HEE) {
		computeCollapseCost(eaEdge);
	}

	//Perform the required number of edge collapses
	while(numValidFaces > faceCnt){
		int select;
		int prevSelect = -1;
		/*
		//Version 1:
		//Randomly select one
		int select;
		bool isValidEdge = false;
		while (!isValidEdge){
			select = rand() % HEE.size();
			std::cout << "Selected value: " << select; //TODO
			if (HEE[select]->isValid) {
				isValidEdge = true;
				std::cout << " is valid"  << std::endl; //TODO
			}
		}
		collapseEdge(HEE[select]);
		*/

		//Version 2:
		//Select based on collapse cost (do least cost first)
		bool isGood = false;
		do {
			select = selectLeastCostEdge();
			std::cout << "Selected edge index: " << select << std::endl;
			std::cout << "With cost: " << HEE[select]->cost << std::endl;
			std::cout << "is valid: " << HEE[select]->isValid << std::endl;
			std::cout << "is dangerous: " << HEE[select]->isDanger << std::endl;

			if (!canCollapse(HEE[select])) {
				if (prevSelect == select) {
					std::cout << "Cannot collapse even the worst case" << std::endl; //TODO
					return;
				}
				else {
					//Redo the loop with updated edge conditions
					prevSelect = select;
				}
			}
			else {
				//Can be collapsed
				isGood = true;
			}

		} while (!isGood);
		std::cout << "Can collapse" << std::endl; //TODO

		collapseEdge(HEE[select]);

		//Keep track of previous to detect uncollapsable cycles
		prevSelect = select;

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

int Mesh::getNumVerticesAdjacentTo(HEVertex* u, HEVertex* v) {
	int counter = 0;
	std::cout << "u: " << u->isValid; debugVertex(u); std::cout << std::endl; //TODO
	std::cout << "v: " << v->isValid; debugVertex(v); std::cout << std::endl; //TODO
	for (HEVertex* eaNeighborOfU : neighborVertices(u)) {
		for (HEVertex* eaNeighborOfV : neighborVertices(v)) {
			if (eaNeighborOfU == eaNeighborOfV) counter++;
		}
	}
	return counter;
}

std::vector<HEEdge*> Mesh::getIncomingEdges(HEVertex* v) {
	std::vector<HEEdge*> incomingEdges;
	HEEdge* firstEdge = v->outEdge->twinEdge;
	HEEdge* currentEdge = v->outEdge->twinEdge;

	do {
		incomingEdges.push_back(currentEdge);
		currentEdge = currentEdge->nextEdge->twinEdge;
	} while (firstEdge != currentEdge);

	for (auto eaEdge : incomingEdges) {
		assert(eaEdge->isValid && eaEdge->endVertex == v);
	}
	
	return incomingEdges;
}

bool Mesh::canCollapse(HEEdge* __edge) {
	std::cout << "test can collapse" << std::endl; //TODO
	//An invalid edge can't be collapsed
	if (!__edge->isValid) {
		std::cout << "is invalid edge" << std::endl; //TODO
		return false;
	}
	std::cout << "not invalid edge" << std::endl; //TODO
	//See: http://stackoverflow.com/questions/27049163/mesh-simplification-edge-collapse-conditions

	//Connectivity checks
	HEVertex* ve = __edge->endVertex;
	HEVertex* vo = __edge->twinEdge->endVertex;

	std::cout << "ve: " << ve->isValid; debugVertex(ve); std::cout << std::endl; //TODO
	std::cout << "vo: " << vo->isValid; debugVertex(vo); std::cout << std::endl; //TODO

	//The vertices must have exactly two common neighbor vertices
	//Otherwise the collapse would cause non-manifold
	//Therefore the edge
	std::cout << "can get ve,vo" << std::endl; //TODO
	if (getNumVerticesAdjacentTo(ve, vo) != 2) {
		std::cout << "Not exactly 2 common vertex neighbours" << std::endl; //TODO
		//So that we don't consider this edge anymore, we set it dangerous
		__edge->isDanger = true;
		return false;
	}
	std::cout << "sharing 2 adj vertices" << std::endl; //TODO

	//Geometry Checks
	//TBC, might not be needed for 3D meshes

	return true;
}

std::vector<HEEdge*> Mesh::getAllNeighbourhoodEdges(HEVertex* v) {
	assert(v->isValid);
	
	std::vector<HEEdge*> allEdges;
	HEEdge* firstOutgoingEdge = v->outEdge;
	HEEdge* currentOutgoingEdge = v->outEdge;

	do {
		allEdges.push_back(currentOutgoingEdge);
		allEdges.push_back(currentOutgoingEdge->twinEdge);
		currentOutgoingEdge = currentOutgoingEdge->twinEdge->nextEdge;
	} while (firstOutgoingEdge != currentOutgoingEdge);

	for (auto eaEdge : allEdges) {
		assert(eaEdge->isValid);
	}

	return allEdges;
}

int Mesh::selectLeastCostEdge() {
	int select = 0;
	for (int index = 1; index < HEE.size(); index++) {
		if (*HEE[index] < *HEE[select]) {
			select = index;
		}
	}
	assert(HEE[select]->isValid);
	return select;
}

std::vector<HEEdge*> Mesh::outgoingEdges(HEVertex* v) {
	std::vector<HEEdge*> outgoing;
	HEEdge* firstOutEdge = v->outEdge;
	HEEdge* currOutEdge = v->outEdge;

	if (!currOutEdge->isValid) {
		std::cout << "outgoingEdges: invalid edge" << std::endl; //TODO
		return outgoing;
	}
	do {
		//The end vertex of this half-edge is a neighbor vertex
		outgoing.push_back(currOutEdge);
		//Get the next outgoing half-edge
		//The next edge of its twin is also an outgoing half-edge from v
		currOutEdge = currOutEdge->twinEdge->nextEdge;

	} while (currOutEdge != firstOutEdge);

	for (HEEdge* e : outgoing) {
		assert(e->isValid);
	}

	return outgoing;
}

void Mesh::computeCollapseCost(HEEdge* e) {
	/*
	The collapse cost is the maximum among either side face,
	which is the minimum dot product of norms amongst all
	neighbouring faces with that side face
	*/
	HEVertex* u = e->prevEdge->endVertex;
	HEVertex* v = e->endVertex;
	float edgeLength = magnitude(u, v);
	float curvature = 0;

	//Find the two faces that are adjacent/opposite e
	std::vector<HEFace*> sideFaces;
	for (HEFace* eaF : neighborFaces(u)) {
		if (faceHasVertex(eaF, v)) {
			sideFaces.push_back(eaF);
		}
	}

	//Use the faces to determine curvature term
	for (HEFace* eaF : neighborFaces(u)) {
		float minCurvature = 1;
		for (HEFace* eaSideFace : sideFaces) {
			//Compute dot product
			float dotProduct = computeDotProduct(eaF, eaSideFace);
			minCurvature = std::min(minCurvature, (1 - dotProduct) / 2.0f);
		}

		curvature = std::max(curvature, minCurvature);
	}

	//The is the cost of collapsing this edge
	e->cost = edgeLength * curvature;
}

float Mesh::magnitude(HEVertex* u, HEVertex* v) {
	//By pythagoras theorem, it is the square root of
	//the sum of the difference in each dimension squared.
	return std::sqrtf(
		std::powf(u->x - v->x, 2.0) +
		std::powf(u->y - v->y, 2.0) +
		std::powf(u->z - v->z, 2.0)
		);
}

float Mesh::computeDotProduct(HEFace* faceA, HEFace* faceB) {
	//Compute the normals first (normalized)
	Vector3 faceANormal = computeNormalizedNormal(faceA);
	Vector3 faceBNormal = computeNormalizedNormal(faceB);

	//The dot product of a 3d vector is the product of sums on each dimension
	return faceANormal.x * faceBNormal.x +
		faceANormal.y * faceBNormal.y +
		faceANormal.z * faceBNormal.z;
}

Vector3 Mesh::computeNormalizedNormal(HEFace* face) {
	//Produce two vectors, to compute the normal
	HEVertex* current = face->edge->endVertex;
	HEVertex* next = face->edge->nextEdge->endVertex;
	HEVertex* prev = face->edge->prevEdge->endVertex;

	//current -> next
	Vector3 a;
	a.x = next->x - current->x;
	a.y = next->y - current->y;
	a.z = next->z - current->z;
	//current -> prev
	Vector3 b;
	b.x = next->x - prev->x;
	b.y = next->y - prev->y;
	b.z = next->z - prev->z;

	//Using coordinate notation formula:
	Vector3 normal;
	normal.x = a.y*b.z - a.z*b.y;
	normal.y = a.z*b.x - a.x*b.z;
	normal.z = a.x*b.y - a.y*b.x;

	return normal;
}

std::pair<HEVertex*, HEVertex*> Mesh::getVertices(HEEdge* e) {
	assert(e && e->isValid);
	return std::pair<HEVertex*, HEVertex*>(
		e->prevEdge->endVertex,
		e->endVertex);
}

void Mesh::debugVertex(HEVertex* v) {
	std::cout << "(" << v->x << ", " << v->y << ", " << v->z << ")";;
}