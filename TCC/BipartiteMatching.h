#pragma once
#include "Graph.h"
constexpr long invalid_flow_edge = -1;

class BipartiteMatching { // begin class BipartiteMatching
public:
  BipartiteMatching(long firstPartitionSize, long secondPartitionSize) // begin constructor
    : flowGraph(firstPartitionSize + secondPartitionSize + 2),
      firstPartitionSize(firstPartitionSize),
      secondPartitionSize(secondPartitionSize),
      target(firstPartitionSize + secondPartitionSize + 1),
      visited(firstPartitionSize + secondPartitionSize + 2, 0)
  {
    // Explicitly default-initialize all edges
    long vertexCount = firstPartitionSize + secondPartitionSize + 2;
    for (long u = 0; u < vertexCount; u++) {
      for (long v = 0; v < vertexCount; v++) {
        flowGraph[u][v] = FordFulkersonEdge(); // Calls default ctor, sets capacity & flow to invalid_flow_edge
      }
    }

    // Then set edges from second partition vertices to target
    for (long i = firstPartitionSize + 1; i < target; i++) {
      flowGraph[i][target] = FordFulkersonEdge(1,0);
      flowGraph[target][i] = FordFulkersonEdge(0,0);
    }
  } // end constructor

  void addEdge(long firstPartitionVertex, long secondPartitionVertex) {
    flowGraph[firstPartitionVertex + 1][secondPartitionVertex + firstPartitionSize + 1] = FordFulkersonEdge(1, 0);
    flowGraph[secondPartitionVertex + 1 + firstPartitionSize][firstPartitionVertex + 1] = FordFulkersonEdge(0, 0);
  }

  void defineInitialCapacity(long firstPartitionVertex, long capacity) {
    flowGraph[source][firstPartitionVertex + 1] = FordFulkersonEdge(capacity, 0);
    flowGraph[firstPartitionVertex + 1][source] = FordFulkersonEdge(0, 0);
  }

  long solve() // Ford Fulkerson algorithm
  {
    long output = 0;
    while (dfs(source)) {
      visitIndex++;
      output++;
    }
    return output;
  }
  
private:
  
  struct FordFulkersonEdge {
    long capacity{invalid_flow_edge };
    long flow{ invalid_flow_edge };

    FordFulkersonEdge(long cap = invalid_flow_edge, long flw = invalid_flow_edge)
      : capacity(cap), flow(flw) {}

    long remainingCapacity() {
      return capacity - flow;
    }
  };
  Graph<FordFulkersonEdge, long>flowGraph;
  std::vector<long> visited;
  long firstPartitionSize;
  long secondPartitionSize;
  long source{0};
  long target;
  long visitIndex{1};
  
  long dfs(long vertex) {
    if (vertex == target) return 1;
    visited[vertex] = visitIndex;

    for (long i = 1; i <= target; i++) {
      if (i == vertex || visited[i] == visitIndex) continue;
      FordFulkersonEdge& outgoingEdge = flowGraph[vertex][i];
      if (outgoingEdge.capacity != invalid_flow_edge && outgoingEdge.remainingCapacity()) {
        if (dfs(i)) {
          augment(vertex, i);
          return 1;
        }
      }
    }
    return 0;
  }

  void augment(long v1, long v2) {
    FordFulkersonEdge& outgoingEdge = flowGraph[v1][v2];
    FordFulkersonEdge& residualEdge = flowGraph[v2][v1];

    outgoingEdge.flow += 1;
    residualEdge.flow -= 1;
  }

}; // end class BipartiteMatching