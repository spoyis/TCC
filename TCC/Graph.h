#pragma once
#include <utility>
#include <iostream>
#include "TypeValues.h"
#include "UnionFind.h"

template <typename edge, typename vertex>
class Graph { // begin class Graph
private:
	class AdjacencyMatrix; // see AdjacencyMatrix.h

public:
	using edge_t = edge;
	using vertex_t = vertex;
	

	Graph(int n) {
		adjMatrix = new(n) AdjacencyMatrix;
		_vertexVal = new UnionFind<vertex>(n);
		_vertexCount = n;
	}

	// Move constructor
	Graph(Graph&& other) noexcept {
		this->adjMatrix = other.adjMatrix;
		this->_vertexVal = other._vertexVal;
		this->_vertexCount = other._vertexCount;

		other.adjMatrix = nullptr;
		other._vertexVal = nullptr;
	}

	// Move assignment operator
	Graph& operator=(const Graph&& other) {
		delete this->adjMatrix;
		delete this->_vertexVal;

		this->_vertexVal = other._vertexVal;
		this->adjMatrix = other.adjMatrix;
		this->_vertexCount = other._vertexCount;

		adjMatrix->graphObj = this;
		other.adjMatrix = nullptr;
		other._vertexVal = nullptr;

		return *this;
	}

	// destructor
	~Graph() {
		delete adjMatrix;
		delete _vertexVal;
	}

	bool operator==(Graph& other) {
		if (other._vertexCount != _vertexCount) return false;

		for (long i = 0; i < _vertexCount; i++)
			for (long j = 0; j < _vertexCount; j++)
				if (this->operator[](i)[j] != other[i][j]) return false;

		return true;
	}

	// returns adjacencymatrix ROW
	AdjacencyMatrix& operator[](int n) { return *(adjMatrix + (n * _vertexCount)); }

	// checks if the given graph is the complete graph or not
	// undirected graph
	bool isUndirectedComplete() {
		const edge* matrix = (edge*)adjMatrix;

		for (int i = 0; i < _vertexCount; i++)
			for (int j = i + 1; j < _vertexCount; j++)
			{
				auto i_root = _vertexVal->findRoot(i);
				auto j_root = _vertexVal->findRoot(j);
				if (i_root == j_root) continue;
				if (!matrix[i_root * _vertexCount + j_root]) return false;
			}
		return true;
	}
	// checks if the given graph is the complete graph or not
	// directed graph
	bool isDirectedComplete() {
		// TODO: MAKE THE LOGIC HERE BETTER, FIX FOR JOINED VERTICES
		const edge* matrix = (edge*)adjMatrix;
		for (int i = 0; i < _vertexCount; i++) {
			for (int j = 0; j < _vertexCount; j++) {
				if (i == j) continue;
				if (!matrix[i * _vertexCount + j])
					return false;
			}
		}
		return true;
	}

	Graph clone() {
		return Graph(*this);
	}

	Graph* cloneHeap() {
		return new Graph(*this);
	}

	// returns copy of current vertex data
	vertex getVertexData(int index) const {
		return _vertexVal->operator[](index);
	}

	void setVertexData(int index, vertex data) {
		_vertexVal->operator[](index) = data;
	}

	void joinVertices(int index1, int index2) {
		int root1 = _vertexVal->findRoot(index1);
		int root2 = _vertexVal->findRoot(index2);

		_vertexVal->unionOp(root1, root2);
		joinEdgeData(std::min(index1, index2), std::max(index1, index2));
	}

	void joinEdgeData(int index1, int index2) {
		for (long i = 0; i < _vertexCount; i++) {
			this->operator[](index1)[i] = this->operator[](index2)[i];
		}

		for (long i = 0; i < _vertexCount; i++) {
			this->operator[](i)[index1] = this->operator[](i)[index2];
		}
	}

	int getRoot(int index) {
		return _vertexVal->findRoot(index);
	}

private:
	AdjacencyMatrix* adjMatrix;
	UnionFind<vertex>* _vertexVal;
	int _vertexCount;

protected:
	// copy constructor privated, helps to not call it by mistake in the code, should you WANT to use it, refer to clone method instead.
	Graph(const Graph& other) {
		_vertexCount = other._vertexCount;
		adjMatrix = other.adjMatrix->clone(_vertexCount);
		_vertexVal = other._vertexVal->clone();
		//adjMatrix->graphObj = this;
	}
}; // end class Graph

#include "AdjacencyMatrix.h"
// G
// escolher dois vertices não adjacentes
// ---> contrai
// ---> adiciona aresta

// DIMACS