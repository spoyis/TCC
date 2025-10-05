#pragma once
#include <functional>
#include "Graph.h"
#include <limits.h>

namespace Coloring::Heuristic
{// begin namespace Coloring::Heuristic

  enum Enums{
    STRATEGY_FIRST_VALID,
    STRATEGY_LOWEST_DEGREE,
    STRATEGY_RANDOM_FROM_CLIQUE,
    STRATEGY_SHARED_NEIGHTBORS,
    
    EVAL_NAIVE,

    EVAL_DEFAULT = EVAL_NAIVE,
    STRATEGY_DEFAULT = STRATEGY_FIRST_VALID,
  };

  template<typename edge>
  using eval = long(*)(edge);

  template<typename edge, typename vertex>
  using strategy = std::pair<long, long>(*)(Graph<edge,vertex>&, eval<edge> evalfunc, std::vector<long>, long, long);

  template<typename edge, typename vertex>
  std::pair<long, long> firstValid(Graph<edge, vertex>& g, eval<edge> evalfunc, std::vector<long> clique, long lastVertex, long vertexCount) {
    for (int i = 0; i < clique.size(); i++)
    {
      long root1 = g.getRoot(clique[i]);
      for (int j = 0; j <= lastVertex; j++) {
        long root2 = g.getRoot(j);
        //std::cout << "CHECKING " << i << " AND " << j << "  == [ " << (root1 == root2) << "]\n";
        if (root1 == root2) continue;
        //std::cout << "DIDN'T CONTINUE; EVAL RESULT == " << evalfunc(g[root1][root2]) << '\n';
        if (evalfunc(g[root1][root2])) return { root1,root2 };
      }
    }

    return {-1,-1};
  }

  // Chooses v1 as min degree vertex (sorted by degree)
 // v2 is first non neighbor
  template<typename edge, typename vertex>
  std::pair<long, long> lowestDegree(Graph<edge, vertex>& g, eval<edge> evalfunc, std::vector<long> clique, long lastVertex, long vertexCount) {
    // Create vector of pairs (degree, vertex) for sorting
    std::vector<std::pair<long, long>> degreeVertexPairs;

    for (int i = 0; i < clique.size(); i++) {
      long root = g.getRoot(clique[i]);
      auto degree = g.getVertexDegree(root);
      degreeVertexPairs.push_back({ degree, root });
    }

    if (degreeVertexPairs.empty()) return { -1, -1 };

    // Sort by degree (ascending order - lowest degrees first)
    std::sort(degreeVertexPairs.begin(), degreeVertexPairs.end());

    // Prepare candidates list once
    std::vector<long> candidates(lastVertex + 1);
    for (long i = 0; i <= lastVertex; ++i) {
      candidates[i] = i;
    }

    // Try each vertex in order of increasing degree
    for (const auto& degreeVertexPair : degreeVertexPairs) {
      long bestVertex = degreeVertexPair.second;
      long bestDeg = degreeVertexPair.first;

      //std::cout << "TRYING VERTEX == " << g.getVertexLabel(bestVertex) << " WITH " << bestDeg << " NEIGHBORS\n";

      // Reset candidates for each vertex attempt
      long n = lastVertex + 1;

      while (n > 0) {
        std::uniform_int_distribution<long> dist(0, n - 1);
        long idx = dist(gen);
        long root2 = g.getRoot(candidates[idx]);
        if (root2 != bestVertex && evalfunc(g[bestVertex][root2])) {
          return { bestVertex, root2 };
        }
        std::swap(candidates[idx], candidates[n - 1]);
        --n;
      }
    }

    return { -1, -1 };
  }

  // v1 is a random vertex from clique
  // v2 is a random vertex not in clique
  template<typename edge, typename vertex>
  std::pair<long, long> randomFromClique(Graph<edge, vertex>& g, eval<edge> evalfunc, std::vector<long> clique, long lastVertex, long vertexCount) {
    std::uniform_int_distribution<long> distClique(0, clique.size() - 1);

    // Try each vertex in clique as root1 (in random order)
    std::vector<long> cliqueOrder = clique;
    std::shuffle(cliqueOrder.begin(), cliqueOrder.end(), gen);

    for (long cliqueVertex : cliqueOrder) {
      long root1 = g.getRoot(cliqueVertex);

      // Try to find a non-neighbor for this root1
      std::vector<long> candidates(lastVertex + 1);
      for (long i = 0; i <= lastVertex; ++i) {
        candidates[i] = i;
      }
      long n = lastVertex + 1;

      while (n > 0) {
        std::uniform_int_distribution<long> dist(0, n - 1);
        long idx = dist(gen);
        long root2 = g.getRoot(candidates[idx]);
        if (root2 != root1 && evalfunc(g[root1][root2])) {
          return { root1, root2 };
        }
        std::swap(candidates[idx], candidates[n - 1]);
        --n;
      }
    }

    return { -1, -1 }; // Should never happen if graph is not complete
  }

  // Returns the two vertices that share the most common neighbors
  template<typename edge, typename vertex>
  std::pair<long, long> mostSharedNeighbors(Graph<edge, vertex>& g, eval<edge> evalfunc, std::vector<long> clique, long lastVertex, long vertexCount) {
    long bestV1 = -1;
    long bestV2 = -1;
    long maxSharedNeighbors = -1;

    for (long i = 0; i <= lastVertex; i++) {
      long root1 = g.getRoot(i);

      for (long j = i + 1; j <= lastVertex; j++) {
        long root2 = g.getRoot(j);

        if (root1 == root2) continue;
        if (!evalfunc(g[root1][root2])) continue;

        long sharedCount = 0;
        for (long k = 0; k <= lastVertex; k++) {
          long root3 = g.getRoot(k);
          if (root3 == root1 || root3 == root2) continue;

          // Check if k is a neighbor of both root1 and root2
          bool isNeighborOf1 = !evalfunc(g[root1][root3]); 
          bool isNeighborOf2 = !evalfunc(g[root2][root3]); 

          if (isNeighborOf1 && isNeighborOf2) {
            sharedCount++;
          }
        }

        if (sharedCount > maxSharedNeighbors) {
          maxSharedNeighbors = sharedCount;
          bestV1 = root1;
          bestV2 = root2;
        }
      }
    }

    return { bestV1, bestV2 };
  }

  template<typename edge>
  long naiveEval(edge e) {
    return !e ? 1 : 0;
  }

}// end namespace Coloring::Heuristic