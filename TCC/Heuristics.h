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
    long root1 = g.getRoot(clique[distClique(gen)]);


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

    return { -1, -1 };
  }

  template<typename edge>
  long naiveEval(edge e) {
    return !e ? 1 : 0;
  }

}// end namespace Coloring::Heuristic