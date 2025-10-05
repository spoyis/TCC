#pragma once
#include <utility>
#include <iostream>
#include "TypeValues.h"
#include "UnionFind.h"
#include <vector>
#include <random>


std::mt19937 gen(5489); // RNG

template <typename edge, typename vertex>
class Graph { // begin class Graph
private:
	class AdjacencyMatrix; // see AdjacencyMatrix.h

public:
	using edge_t = edge;
	using vertex_t = vertex;
	using clique_t = std::vector<long>;
	

	Graph(int n) {
		adjMatrix = new(n) AdjacencyMatrix;
		_vertexVal = new UnionFind<vertex>(n);
		_vertexLabels = std::make_shared<std::vector<std::string>>(n);
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
	Graph& operator=(Graph&& other) {
		delete this->adjMatrix;
		delete this->_vertexVal;

		this->_vertexVal = other._vertexVal;
		this->adjMatrix = other.adjMatrix;
		this->_vertexCount = other._vertexCount;
		this->_vertexLabels = std::move(other._vertexLabels);

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

	// returns REFERENCE of current vertex data
	vertex& getVertexData(int index) const {
		return _vertexVal->operator[](index);
	}

	auto getVertexCount() { return _vertexCount; };

	std::string getVertexLabel(int index) { return (*_vertexLabels)[index]; }

	long getEdgeCount(){
		long output = 0;
		for (long i = 0; i < _vertexCount; i++) {
			auto rootI = this->getRoot(i);
			if (rootI != i) continue;
			for (long j = 0; j < _vertexCount; j++) {
				auto rootJ = this->getRoot(j);
				if (rootJ != j || rootJ == rootI) continue;
				if (this->operator[](rootI)[rootJ]) output++;
			}
		}

		return output;
	}

	void setVertexLabel(std::string label, int index) {
		(*_vertexLabels)[index] = label;
	}

	void setVertexData(int index, vertex data) {
		_vertexVal->operator[](index) = data;
	}

	void joinVertices(int index1, int index2) {
		int root1 = _vertexVal->findRoot(index1);
		int root2 = _vertexVal->findRoot(index2);
		//if (root1 == root2) std::cout << "you're trying to join a vertex thats already been joined, you sure?\n";

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

	long findCliqueGreedy() {

	}

	clique_t findCliqueGrasp(const clique_t& initialClique = {}, int maxIterations = 5) {
	clique_t bestClique;

	for (int iter = 0; iter < maxIterations; iter++) {
		clique_t clique;

		if (initialClique.empty()) {
			std::uniform_int_distribution<> distrib(0, _vertexCount - 1);
			int current = getRoot(distrib(gen));
			clique.push_back(current);
		} else {
			clique = initialClique;
		}

		// Candidate set = all other vertices
		std::vector<int> candidates;
		for (int i = 0; i < _vertexCount; i++) {
			int root = getRoot(i);
			if (std::find(clique.begin(), clique.end(), i) == clique.end() && root == i) {
				candidates.push_back(i);
			}
		}

		while (!candidates.empty()) {
			// Build restricted candidate list (RCL)
			struct CandidateScore {
				int v;
				int score;
			};
			std::vector<CandidateScore> scored;

			for (int v : candidates) {
				// Check if v is adjacent to all current clique members
				bool isNeighbor = true;
				for (int c : clique) {
					if (!this->operator[](v)[c]) {
						isNeighbor = false;
						break;
					}
				}
				if (!isNeighbor) continue;

				// Heuristic: degree of v
				int degree = 0;
				for (int k = 0; k < _vertexCount; k++) {
					int rootK = getRoot(k);
					if (rootK != k)continue;
					if (rootK != v && this->operator[](v)[rootK]) {
						degree++;
					}
				}

				scored.push_back({ v, degree });
			}

			if (scored.empty()) break;

			// Find max score
			int maxScore = -1;
			for (auto& s : scored) {
				if (s.score > maxScore) maxScore = s.score;
			}

			// RCL: (top scorers within 80% of max)
			std::vector<int> rcl;
			for (auto& s : scored) {
				if (s.score >= (maxScore * 0.8)) {
					rcl.push_back(s.v);
				}
			}

			// Pick random from RCL
			std::uniform_int_distribution<> d(0, (int)rcl.size() - 1);
			int chosen = rcl[d(gen)];
			clique.push_back(chosen);

			// Update candidates: keep only vertices adjacent to all in clique
			std::vector<int> newCandidates;
			for (int v : candidates) {
				if (v == chosen) continue;
				bool ok = true;
				for (int c : clique) {
					if (!this->operator[](v)[c]) {
						ok = false;
						break;
					}
				}
				if (ok) newCandidates.push_back(v);
			}
			candidates.swap(newCandidates);
		}

		if (clique.size() > bestClique.size()) {
			bestClique = clique;
		}
	}

	return bestClique;
}


	clique_t findCliqueRandom(const clique_t& initialClique = {}) {
		clique_t clique;

		if (initialClique.empty()) {
			// Create initial clique with random vertex if none provided
			std::uniform_int_distribution<> distrib(0, _vertexCount - 1);
			int randomNumber = this->getRoot(distrib(gen));
			clique = { randomNumber };
		}
		else {
			// Use provided initial clique
			clique = initialClique;
		}

    
    // Create candidates array with all vertices
    std::vector<long> candidates(_vertexCount);
    for (long i = 0; i < _vertexCount; i++) {
        candidates[i] = i;
    }
    
    long n = _vertexCount;
    while (n > 0) {
        std::uniform_int_distribution<long> dist(0, n - 1);
        long idx = dist(gen);
        auto rootI = this->getRoot(candidates[idx]);
        
        // Skip if not a root vertex
        if (rootI != candidates[idx]) {
            std::swap(candidates[idx], candidates[n - 1]);
            --n;
            continue;
        }
        
        // Check if this vertex can be added to the clique
        bool flag = true;
        for (long j = 0; j < clique.size(); j++) {
            int cliqueVertex = clique[j];
            if (cliqueVertex == rootI || !this->operator[](rootI)[cliqueVertex]) {
                flag = false;
                break;
            }
        }
        
        if (flag) {
            clique.push_back(rootI);
        }
        
        // Remove this candidate from consideration
        std::swap(candidates[idx], candidates[n - 1]);
        --n;
    }
    
    return clique;
	}


	template <class Vector>
	typename std::enable_if<is_specialization<Vector, std::vector>::value, bool>::type
	areJoinableHelper(int index1, int index2) {
		UnionFind<vertex>& unionfind = *_vertexVal;

		for (const auto& color_1 : unionfind[index1]) {
			for (const auto& color_2 : unionfind[index2]) {
				if (color_1 == color_2) {
					return true;
				}
			}
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
		vertex& v2 = unionfind[index2];

		joinVertexDataHelper<vertex>(v1, v2);
	}
	// -- Template method for any std::vector, calculates intersection.
	template <class Vector>
	typename std::enable_if<is_specialization<Vector, std::vector>::value>::type
	joinVertexDataHelper(Vector& v1, const Vector& v2) {
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
		std::cout << "Given type is not an std::vector... method not implemented lol (ignore this line)" << std::endl;
	}
	// END JOINVERTEXDATA METHODS 

	void joinEdgeData(int index1, int index2) {
		for (long i = 0; i < _vertexCount; i++) {
			this->operator[](index1)[i] |= this->operator[](index2)[i];
		}

		for (long i = 0; i < _vertexCount; i++) {
			this->operator[](i)[index1] |= this->operator[](i)[index2];
		}
	}

	int getRoot(int index) {
		return _vertexVal->findRoot(index);
	}

	long getVertexDegree(int index) {
		long output = 0;
		for (long i = 0; i < _vertexCount; i++) {
			if (this->operator[](index)[i] && i != index) output++;
		}
		return output;
	}

private:
	AdjacencyMatrix* adjMatrix;
	UnionFind<vertex>* _vertexVal;
	int _vertexCount;
	std::shared_ptr<std::vector<std::string>> _vertexLabels;


protected:
	// copy constructor privated, helps to not call it by mistake in the code, should you WANT to use it, refer to clone method instead.
	Graph(const Graph& other) {
		_vertexCount = other._vertexCount;
		adjMatrix = other.adjMatrix->clone(_vertexCount);
		_vertexVal = other._vertexVal->clone();
		_vertexLabels = other._vertexLabels;
		//adjMatrix->graphObj = this;
	}
}; // end class Graph

#include "AdjacencyMatrix.h"
// DIMACS