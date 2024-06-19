#pragma once
#include <utility>
#include <iostream>
#include "TypeValues.h"
#include "UnionFind.h"
#include <vector>

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

	auto getVertexCount() { return _vertexCount; };

	void setVertexData(int index, vertex data) {
		_vertexVal->operator[](index) = data;
	}

	void joinVertices(int index1, int index2) {
		int root1 = _vertexVal->findRoot(index1);
		int root2 = _vertexVal->findRoot(index2);

		_vertexVal->unionOp(root1, root2);
		joinEdgeData(std::min(index1, index2), std::max(index1, index2));
		joinVertexData(std::min(index1, index2), std::max(index1, index2));
	}


	bool areJoinable(int index1, int index2){
		return areJoinableHelper<vertex>(index1, index2);
	}

	std::vector<long> neighboringVertices(int vertex) {
		std::vector<long> output;
		const edge* matrix = (edge*)adjMatrix;
		for (long i = 0; i < _vertexCount; i++) {
			if (i == vertex) continue;
			if (matrix[vertex * _vertexCount + i]) output.push_back(i);
		}

		return output;
	}


	template <class Vector>
	typename std::enable_if<is_specialization<Vector, std::vector>::value, bool>::type
	areJoinableHelper(int index1, int index2) {
		UnionFind<vertex>& unionfind = *_vertexVal;
		long minSize = 1;
		long size = 0;

		for (auto color_1: unionfind[index1]) {
			for (auto color_2 : unionfind[index2]) {
				if (color_1 == color_2) { size++; break; }
			}

			if (size == minSize) return true;
		}

		return false;
	}

	template<typename vertex_t>
	typename std::enable_if<!is_specialization<vertex_t, std::vector>::value, bool>::type
	areJoinableHelper(int index1, int index2) {
			return true;
	}

	// BEGIN JOINVERTEXDATA METHODS
	void joinVertexData(int index1, int index2) {
		UnionFind<vertex>& unionfind = *_vertexVal;
		vertex& v1 = unionfind[index1];
		vertex& v2 = unionfind[index1];

		joinVertexDataHelper<vertex>(v1, v2);
	}
	// -- Template method for any std::vector, calculates intersection.
	template <class Vector>
	typename std::enable_if<is_specialization<Vector, std::vector>::value>::type
	joinVertexDataHelper(Vector v1, const Vector v2) {
		auto right = v1.size();
		long left = 0;

		while(left < right) {
			auto currentElement = v1[left];
			bool found = false;
			for (auto otherElement : v2) {
				if (currentElement == otherElement) {
					found = true;
					break;
				}
			}
			if (!found) {
				std::swap(v1[left], v1[--right]);
			}
			else left++;
		}
		v1.resize(right);
	}
	// -- OVERLOADED FOR NON STD::VECTORS
	template<typename vertex_t>
	typename std::enable_if<!is_specialization<vertex_t, std::vector>::value>::type
	joinVertexDataHelper(const vertex_t& v1, const vertex_t&v2) {
		std::cout << "Given type is not an std::vector... method not implemented lol" << std::endl;
	}
	// END JOINVERTEXDATA METHODS 

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
// DIMACS