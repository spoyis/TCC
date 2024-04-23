#pragma once
#include <functional>
#include "Graph.h"

namespace Coloring::Heuristic
{// begin namespace Coloring::Heuristic

  enum Enums{
    STRATEGY_FIRST_VALID, EVAL_NAIVE,

    EVAL_DEFAULT = EVAL_NAIVE,
    STRATEGY_DEFAULT = STRATEGY_FIRST_VALID,
  };

  template<typename edge>
  using eval = long(*)(edge);

  template<typename edge, typename vertex>
  using strategy = std::pair<long, long>(*)(Graph<edge,vertex>&, eval<edge>, long);

  template<typename edge, typename vertex>
  std::pair<long, long> firstValid(Graph<edge, vertex>& g, eval<edge> evalfunc, long lastVertex) {
    for (int i = 0; i <= lastVertex; i++)
      for (int j = i + 1; j <= lastVertex; j++) {
        long root1 = g.getRoot(i);
        long root2 = g.getRoot(j);
        if (root1 == root2) continue;
        if (evalfunc(g[root1][root2]) and g.areJoinable(root1, root2)) return { root1,root2 };
      }
        

    return {-1,-1};
  }

  template<typename edge>
  long naiveEval(edge e) {
    return !e ? 1 : 0;
  }

}// end namespace Coloring::Heuristic