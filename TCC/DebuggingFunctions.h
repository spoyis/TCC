#pragma once
#include <vector>
#include "Graph.h"

template <typename edge, typename vertex>
bool isThisAClique(std::vector<long> vertices, Graph<edge, vertex>* graph) {
  if (vertices.size() <= 1)
    return true;


  for (long i = 0; i < vertices.size(); ++i) {
    long v1 = vertices[i];
    for (long j = i + 1; j < vertices.size(); ++j) {
      long v2 = vertices[j];
      if (!(*graph)[v1][v2]) {
        return false;
      }
    }
  }

  return true; 
}