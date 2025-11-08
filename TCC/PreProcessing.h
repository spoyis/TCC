#pragma once
#include "Graph.h"

namespace Coloring {
  template <class graph_t>
  class Checker; // forward declaration
}

namespace Coloring::Optimization {

  template <typename edge, typename vertex>
  void sameTimeDifferentRoomsOptimization(Graph<edge, vertex> *_originalGraph, std::vector<std::vector<long>>& roomData) {
    std::vector<std::pair<long, long>> nonNeighbors;
    Graph<edge, vertex>& graph = *_originalGraph;
    std::vector<long> roots = graph.getRoots();
    for (long i = 0; i < roots.size(); i++)
    {
      long rootI = roots[i];
      for (long j = 0; j < roots.size(); j++) {
        long rootJ = roots[j];
        if (i == j) continue;
        if (roomData[rootI].size() > 1 || roomData[rootJ].size() > 1) continue;

        if (!graph[rootI][rootJ]) nonNeighbors.push_back({ i,j });
      }
    }

    long counter = 0;
    for (auto& vertexPair : nonNeighbors) {
      auto v1 = vertexPair.first;
      auto v2 = vertexPair.second;

      if (roomData[v1][0] == roomData[v2][0]) {
        graph[v1][v2] = { 1 };
        graph[v2][v1] = { 1 };

        counter++;
        //std::cout << "[OPTIMIZING] CLASS " << graph.getVertexLabel(v1) << " AND " << graph.getVertexLabel(v2) << " SHARE THE SAME ONE ROOM\n";
        //std::cout << "[DEBUG] Added edge between " << v1 << " and " << v2 << std::endl;

      }
    }

    std::cout << "[PREPROCESS][DIFFERENTROOMSSAMETIME] created " << counter << " edges.\n";

  }

  template <typename edge, typename vertex>
  void preprocessNonJoinableEdges(Graph<edge, vertex>*_originalGraph) {
    Graph<edge, vertex>& graph = *_originalGraph;
    long vertexCount = graph.getVertexCount();
    std::vector<std::pair<long, long>> nonNeighbors;

    // collect all non-adjacent vertex pairs
    for (long i = 0; i < vertexCount; i++) {
      long rootI = graph.getRoot(i);
      if (rootI < i) continue;

      for (long j = 0; j < vertexCount; j++) {
        long rootJ = graph.getRoot(j);
        if (rootJ < j) continue;
        if (i != j && !graph[rootI][rootJ]) {
          nonNeighbors.emplace_back(i, j);
        }
      }
    }

    long optimizationCounter = 0;
    for (auto& [v1, v2] : nonNeighbors) {
      if (!graph.areJoinable(v1, v2)) {
        graph[v1][v2] = { 1 };
        graph[v2][v1] = { 1 };
        optimizationCounter++;
      }
    }

    if (optimizationCounter != 0) {
      std::cout << "[PREPROCESS][NON_JOINABLE] Added "
        << optimizationCounter << " edges.\n";
    }
  }

  template <typename edge, typename vertex>
  void propagateSingleColorConstraints(Graph<edge, vertex>* _originalGraph) {
    Graph<edge, vertex>& graph = *_originalGraph;
    long vertexCount = graph.getVertexCount();
    bool changed = true;
    long propagationCounter = 0;

    while (changed) {
      changed = false;

      for (long i = 0; i < vertexCount; i++) {
        auto& colors = graph.getVertexData(i);
        if (colors.size() == 1) {
          auto fixedColor = colors.front();

          for (long j = 0; j < vertexCount; j++) {
            if (i == j || !graph[i][j]) continue;

            if (graph.removeColor(j, fixedColor)) {
              propagationCounter++;
              changed = true;
            }
          }
        }
      }
    }

    if (propagationCounter != 0) {
      std::cout << "[PREPROCESS][SINGLECOLOR] Propagated "
        << propagationCounter
        << " color removals.\n";
    }
  }


  template <typename graph_t>
  long contractUnitIntersectionSingles(typename Coloring::Checker<graph_t>::Zykov* zykovRoot) {
    auto& g = *zykovRoot->getGraph();
    bool changed = true;
    long counter = 0;

    while (changed) {
      changed = false;
      const long n = g.getVertexCount();

      for (long i = 0; i < n; ++i) {
        long ri = g.getRoot(i);
        if (ri != i) continue;

        auto& Li = g.getVertexData(ri);
        if (Li.size() != 1) continue;

        for (long j = i + 1; j < n; ++j) {
          long rj = g.getRoot(j);
          if (rj != j || ri == rj || g[ri][rj]) continue;

          auto& Lj = g.getVertexData(rj);
          if (Lj.size() != 1) continue;

          if (Li[0] == Lj[0]) {
            g.joinVertices(ri, rj);
            counter++;
            changed = true;
            break;
          }
        }
        if (changed) break;
      }
    }

    zykovRoot->setVertexCount(g.getVertexCount() - counter);

    std::cout << "[PREPROCESS][SINGLEINTERSECTION] Propagated "
      << counter << " vertex joinings.\n";
    return counter;
  }
}