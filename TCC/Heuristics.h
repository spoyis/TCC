#pragma once
#include <functional>
#include "Graph.h"
#include <limits.h>

namespace Coloring::Heuristic
{// begin namespace Coloring::Heuristic

  enum Enums{
    STRATEGY_FIRST_VALID,
    STRATEGY_LOWEST_DEGREE,
    
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
      for (int j = i + 1; j <= lastVertex; j++) {
        long root2 = g.getRoot(j);
        //std::cout << "CHECKING " << i << " AND " << j << "  == [ " << (root1 == root2) << "]\n";
        if (root1 == root2) continue;
        //std::cout << "DIDN'T CONTINUE; EVAL RESULT == " << evalfunc(g[root1][root2]) << '\n';
        if (evalfunc(g[root1][root2])) return { root1,root2 };
      }
    }

    return {-1,-1};
  }

  // Chooses v1 as min degree vertex
  // v2 is first non neighbor
  template<typename edge, typename vertex>
  std::pair<long, long> lowestDegree(Graph<edge, vertex>& g, eval<edge> evalfunc, std::vector<long> clique, long lastVertex, long vertexCount) {
    long bestVertex = -1;
    long bestDeg = LONG_MAX;
    
    for (int i = 0; i < clique.size(); i++) {
      long root = g.getRoot(clique[i]);
      auto degree = g.getVertexDegree(root);
      if (degree < vertexCount - 1 && degree < bestDeg) {
        bestDeg = degree;
        bestVertex = root;
      }
    }

    
    if (bestVertex == -1) return { -1, -1 };
    //std::cout << "BEST DEGREE VERTEX == " << g.getVertexLabel(bestVertex) << " WITH " << bestDeg << "NEIGHBORS\n";
    std::vector<long> candidates(lastVertex + 1);
    for (long i = 0; i <= lastVertex; ++i) {
      candidates[i] = i;
    }

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

    return { -1, -1 };
  }

  template<typename edge>
  long naiveEval(edge e) {
    return !e ? 1 : 0;
  }

}// end namespace Coloring::Heuristic