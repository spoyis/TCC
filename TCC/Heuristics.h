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
    STRATEGY_RANDOM_COLOR_INTERSECTION,

    EVAL_NAIVE,

    EVAL_DEFAULT = EVAL_NAIVE,
    STRATEGY_DEFAULT = STRATEGY_FIRST_VALID,
  };

  template<typename edge>
  using eval = long(*)(edge);

  template<typename edge, typename vertex>
  using strategy = std::pair<long, long>(*)(
    Graph<edge, vertex>&,
    eval<edge> evalfunc,
    const std::vector<std::vector<long>>& roomData,
    std::vector<long> clique,
    long lastVertex,
    long vertexCount
    );

  template<typename edge, typename vertex>
  std::pair<long, long> firstValid(
    Graph<edge, vertex>& g,
    eval<edge> evalfunc,
    const std::vector<std::vector<long>>&,
    std::vector<long> clique,
    long lastVertex,
    long vertexCount) {
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
  std::pair<long, long> lowestDegree(
    Graph<edge, vertex>& g,
    eval<edge> evalfunc,
    const std::vector<std::vector<long>>& ,
    std::vector<long> clique,
    long lastVertex,
    long vertexCount) {
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
  std::pair<long, long> randomFromClique(
    Graph<edge, vertex>& g,
    eval<edge> evalfunc,
    const std::vector<std::vector<long>>&,
    std::vector<long> clique,
    long lastVertex,
    long vertexCount) {
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
  std::pair<long, long> mostSharedNeighbors(
    Graph<edge, vertex>& g,
    eval<edge> evalfunc,
    const std::vector<std::vector<long>>& roomData,
    std::vector<long> clique,
    long lastVertex,
    long vertexCount
  ) {
    struct PairScore {
      long v1;
      long v2;
      long sharedNeighbors;
      long roomIntersection;
    };

    std::vector<PairScore> scores;

    for (long i = 0; i <= lastVertex; i++) {
      long root1 = g.getRoot(i);

      for (long j = i + 1; j <= lastVertex; j++) {
        long root2 = g.getRoot(j);
        if (root1 == root2) continue;
        if (!evalfunc(g[root1][root2])) continue; 

        // Shared neighbor score 
        long sharedCount = 0;
        for (long k = 0; k <= lastVertex; k++) {
          long root3 = g.getRoot(k);
          if (root3 == root1 || root3 == root2) continue;

          bool isNeighborOf1 = !evalfunc(g[root1][root3]);
          bool isNeighborOf2 = !evalfunc(g[root2][root3]);
          if (isNeighborOf1 && isNeighborOf2) sharedCount++;
        }

        // Room intersection score
        const auto& R1 = roomData[root1];
        const auto& R2 = roomData[root2];
        long roomShared = 0;
        for (auto r1 : R1)
          for (auto r2 : R2)
            if (r1 == r2) roomShared++;

        scores.push_back({ root1, root2, sharedCount, roomShared });
      }
    }

    if (scores.empty()) return { -1, -1 };

    // Sort primarily by shared neighbors, then room intersection
    std::sort(scores.begin(), scores.end(), [](const PairScore& a, const PairScore& b) {
      if (a.sharedNeighbors == b.sharedNeighbors)
        return a.roomIntersection > b.roomIntersection;
      return a.sharedNeighbors > b.sharedNeighbors;
      });

    // Top 1% random pick
    size_t topCount = std::max<size_t>(1, scores.size() / 100);
    std::uniform_int_distribution<size_t> dist(0, topCount - 1);

    auto& chosen = scores[dist(gen)];
    return { chosen.v1, chosen.v2 };
  }



  // Randomly selects a non-adjacent vertex pair with the largest color-list intersection
  template<typename edge, typename vertex>
  std::pair<long, long> randomLargestColorIntersection(
    Graph<edge, vertex>& g,
    eval<edge> evalfunc,
    const std::vector<std::vector<long>>& ,
    std::vector<long> clique,
    long lastVertex,
    long vertexCount)
  {
    long maxIntersection = -1;
    std::vector<std::pair<long, long>> bestPairs;

    for (long i = 0; i <= lastVertex; i++) {
      long root1 = g.getRoot(i);

      for (long j = i + 1; j <= lastVertex; j++) {
        long root2 = g.getRoot(j);
        if (root1 == root2) continue;
        if (!evalfunc(g[root1][root2])) continue; // skip adjacent


        int intersection = g.colorIntersectionSize(root1, root2);
        if (intersection == 0) continue;

        if (intersection > maxIntersection) {
          maxIntersection = intersection;
          bestPairs.clear();
          bestPairs.emplace_back(root1, root2);
        }
        else if (intersection == maxIntersection) {
          bestPairs.emplace_back(root1, root2);
        }
      }
    }

    if (bestPairs.empty()) return { -1, -1 };

    std::uniform_int_distribution<size_t> dist(0, bestPairs.size() - 1);
    return bestPairs[dist(gen)];
  }

  template<typename edge>
  long naiveEval(edge e) {
    return !e ? 1 : 0;
  }

}// end namespace Coloring::Heuristic