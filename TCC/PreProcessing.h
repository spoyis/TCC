#pragma once
#include "Graph.h"


namespace Coloring::Validator {
  template <typename edge, typename vertex>
  bool canJoinBasedOnRoomConstraint(
    long v1,
    long v2,
    Graph<edge, vertex>& g,
    const std::vector<std::vector<long>>& roomData
  );
}

namespace Coloring {
  template <class graph_t>
  class Checker; // forward declaration
}

namespace Coloring::Optimization {

  template <typename edge, typename vertex>
  void sameTimeDifferentRoomsOptimization(Graph<edge, vertex> *_originalGraph, std::vector<std::vector<long>>& roomData, bool shouldPrint = true) {
    Graph<edge, vertex>& graph = *_originalGraph;
    std::vector<long> roots = graph.getRoots();
    long counter = 0;

    for (long i = 0; i < roots.size(); i++)
    {
      long rootI = roots[i];
      for (long j = i + 1; j < roots.size(); j++) {
        long rootJ = roots[j];
        if (graph[rootI][rootJ]) continue;  

        if (Validator::canJoinBasedOnRoomConstraint(rootI, rootJ, graph, roomData) == false) {
          graph[rootI][rootJ] = { 1 };
          graph[rootJ][rootI] = { 1 };

          counter++;
        }
      }
    }
    if(shouldPrint)
    std::cout << "[PREPROCESS][DIFFERENTROOMSSAMETIME] created " << counter << " edges.\n";

  }

  template <typename edge, typename vertex>
  void preprocessNonJoinableEdges(Graph<edge, vertex>*_originalGraph, bool shouldPrint = true) {
    Graph<edge, vertex>& graph = *_originalGraph;
    long vertexCount = graph.getVertexCount();
    std::vector<std::pair<long, long>> nonNeighbors;

    // collect all non-adjacent vertex pairs
    for (long i = 0; i < vertexCount; i++) {
      long rootI = graph.getRoot(i);
      if (rootI < i) continue;

      for (long j = i + 1; j < vertexCount; j++) {
        long rootJ = graph.getRoot(j);
        if (rootJ < j) continue;
        if (i != j && !graph[rootI][rootJ]) {
          nonNeighbors.emplace_back(rootI, rootJ);
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

    if (optimizationCounter != 0 && shouldPrint) {
      std::cout << "[PREPROCESS][NON_JOINABLE] Added "
        << optimizationCounter << " edges.\n";
    }
  }

  template <typename edge, typename vertex>
  void propagateSingleColorConstraints(Graph<edge, vertex>* _originalGraph, bool shouldPrint = true) {
    Graph<edge, vertex>& graph = *_originalGraph;
    
    bool changed = true;
    long propagationCounter = 0;
    auto roots = graph.getRoots();
    long vertexCount = roots.size();
    while (changed) {
      changed = false;

      for (long i = 0; i < vertexCount; i++) {
        long ri = roots[i];
        auto& colors = graph.getVertexData(ri);
        if (colors.size() == 1) {
          auto fixedColor = colors.front();

          for (long j = 0; j < vertexCount; j++) {
            long rj = roots[j];
            if (ri == rj || !graph[ri][rj]) continue;

            if (graph.removeColor(rj, fixedColor)) {
              propagationCounter++;
              changed = true;
            }
          }
        }
      }
    }

    if (propagationCounter != 0 && shouldPrint) {
      std::cout << "[PREPROCESS][SINGLECOLOR] Propagated "
        << propagationCounter
        << " color removals.\n";
    }
  }


  template <typename graph_t>
  long contractUnitIntersectionSingles(typename Coloring::Checker<graph_t>::Zykov* zykovRoot, std::vector<std::vector<long>>& roomData, bool shouldPrint = true) {
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

    zykovRoot->setVertexCount(g.getRoots().size());
    if(counter > 0 && shouldPrint)
    std::cout << "[PREPROCESS][SINGLEINTERSECTION] Propagated "
      << counter << " vertex joinings.\n";
    return counter;
  }
}