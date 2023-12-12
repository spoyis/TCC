#include <iostream>

template <typename edge, typename vertex>
class Graph<edge, vertex>::AdjacencyMatrix { // begin class AdjacencyMatrix
public:
	edge data;

	static void* operator new(size_t size, int vertices) {
		auto p = malloc(size * (vertices * vertices));
		return p;
	}

	AdjacencyMatrix* clone(const int vertexCount) {
		AdjacencyMatrix* copy = new(vertexCount) AdjacencyMatrix;
		const auto adjMatrixSize = vertexCount * vertexCount;
		for (long i = 0; i < adjMatrixSize; i++) {
			*(&copy->data + i) = *(&data + i);
		}

		return copy;
	}

	edge& operator[](const int n) {
		return *(&data + n); 
	}

	edge  operator[](const int n) const {
		return  (edge) *(&data + n); 
	}
}; // end class AdjacencyMatrix