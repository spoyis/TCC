#pragma once
#include "Heuristics.h"
#include "Graph.h"
#include "BipartiteMatching.h"

namespace Coloring {// begin namespace Coloring


  template<class graph_t, typename color>
  class Checker {// begin class Coloring::Checker
  public:
    Checker(graph_t* graph, std::vector<std::vector<long>> roomData) 
      : roomData(graph->getVertexCount()),
        optimizationStrategies( 1, false )
    {
      _originalGraph = graph;
      _initialVertexCount = _bestAnswer = _originalGraph->getVertexCount();

      for (long vertex = 0; vertex < _initialVertexCount; vertex++) {
        this->roomData[vertex] = std::move(roomData[vertex]);
      }
    }

    enum OptimizationStrategy {
      ROOT_NODE_SAME_TIME_DIFFERENT_ROOMS = 0,
    };

    Checker(graph_t* graph, std::vector<std::vector<long>> roomData, std::vector<OptimizationStrategy> strategies)
      : roomData(graph->getVertexCount()),
        optimizationStrategies(1, false) 
    {
      _originalGraph = graph;
      _initialVertexCount = _bestAnswer = _originalGraph->getVertexCount();

      for (long vertex = 0; vertex < _initialVertexCount; vertex++) {
        this->roomData[vertex] = std::move(roomData[vertex]);
      }

      for (auto& strategy : strategies) {
        this->optimizationStrategies[strategy] = true;
      }
    }

    int run() {
      if (isEnabled(ROOT_NODE_SAME_TIME_DIFFERENT_ROOMS)) {
        std::cout << optimizationStrategies[0] << " THE FUCK\n";
        sameTimeDifferentRoomsOptimization();
      }

      Zykov* root = new Zykov(_originalGraph->cloneHeap(), _initialVertexCount, this);
      int result = root->run(getStrategy(), getEval());
      delete root;
      return result;
    }

  private:
    class Zykov;
    using edge = typename graph_t::edge_t;
    using vertex = typename graph_t::vertex_t;

    long _bestAnswer;
    long _initialVertexCount;
    Graph<edge, vertex>* _originalGraph;
    Heuristic::Enums strategy{ Heuristic::STRATEGY_DEFAULT };
    Heuristic::Enums eval{ Heuristic::EVAL_DEFAULT };
    std::vector<bool> optimizationStrategies;
    std::vector<std::vector<long>> roomData;

    Heuristic::strategy<edge,vertex> getStrategy() {
      switch (strategy) {
      case Heuristic::STRATEGY_FIRST_VALID:
        return &Heuristic::firstValid<edge, vertex>;
      }
    }

    Heuristic::eval<edge> getEval() {
      switch (eval) {
      case Heuristic::EVAL_NAIVE:
        return &Coloring::Heuristic::naiveEval<edge>;
      }
    }

    bool isEnabled(OptimizationStrategy givenStrategy) {
      return optimizationStrategies[givenStrategy];
    }

    void sameTimeDifferentRoomsOptimization() {
      std::vector<std::pair<long, long>> nonNeighbors;
      Graph<edge, vertex>& graph = *_originalGraph;

      for (long i = 0; i < _initialVertexCount; i++)
        for (long j = 0; j < _initialVertexCount; j++) {
          if (i == j) continue;
          if (roomData[i].size() > 1 || roomData[j].size() > 1) continue;

          if (!graph[i][j]) nonNeighbors.push_back({ i,j });
        }
      for (auto& vertexPair : nonNeighbors) {
        auto v1 = vertexPair.first;
        auto v2 = vertexPair.second;

        if (roomData[v1][0] == roomData[v2][0]) {
          _originalGraph[v1][v2] = {1};
          _originalGraph[v2][v1] = {1};
        }
      }
    }

  }; // end class Coloring::Checker


  template<class graph_t, typename color>
  class Checker<graph_t, color>::Zykov  
  { // begin class Coloring::Checker::Zykov
    enum PruneMotive {
      not_pruned, bad_clique_contraction, no_valid_intermediate_coloring, no_valid_final_coloring
    };

  private:
    long _currentVertexCount;
    long _rightmostVertexIndex;
    long _bestAnswer;
    Graph<edge, vertex>* _graph;
    Checker* _checkerPtr;
    bool _isValidNode{true};
    PruneMotive pruneMotive{ not_pruned };
    
    
    // consider saving information to make/unmake changes to graph when searching through zykov tree.

  public:
    Zykov(Graph<edge,vertex>* graph, int size, Checker* ptr) {
      _graph = graph;
      _currentVertexCount = _bestAnswer = size;
      _rightmostVertexIndex = size - 1;
      _checkerPtr = ptr;
    }

    ~Zykov() { delete _graph; }

    long run(Heuristic::strategy<edge,vertex> heuristic, Heuristic::eval<edge> eval) {
      auto nonNeighboringVertices = heuristic(*_graph, eval, _rightmostVertexIndex); // std::pair<long,long>

      // heuristic returns -1 when no valid vertex pair is found.
      if(nonNeighboringVertices.first != -1){ 
        Zykov contractVertices(*_graph, nonNeighboringVertices, this);
        Zykov addEdge(*_graph, nonNeighboringVertices, edge{1}, this);

        _bestAnswer = std::min(_bestAnswer, std::min(
          contractVertices.run(heuristic, eval),
          addEdge.run(heuristic, eval)
        ));
      }

      return _graph->isUndirectedComplete() ? std::min(_currentVertexCount, _bestAnswer) : _bestAnswer;
    }

  private:

    // add edge constructor
    Zykov(Graph<edge,vertex>& g, std::pair<long,long> vertices, edge edgeValue, Zykov* parent) {
      _graph = g.cloneHeap();

      Graph<edge, vertex>& graph = *_graph;
      graph[vertices.first][vertices.second] = edgeValue;
      graph[vertices.second][vertices.first] = edgeValue;
      _checkerPtr = parent->_checkerPtr;
      _currentVertexCount = parent->_currentVertexCount;
      _bestAnswer = parent->_bestAnswer;
      _rightmostVertexIndex = parent->_rightmostVertexIndex;
    }

    // contract vertices constructor
    Zykov(Graph<edge, vertex>& g, std::pair<long, long> vertices, Zykov* parent) {
      _graph = g.cloneHeap();

      Graph<edge, vertex>& graph = *_graph;
      
      graph.joinVertices(vertices.first, vertices.second);
      _checkerPtr = parent->_checkerPtr;
      _currentVertexCount = parent->_currentVertexCount - 1;
      _bestAnswer = parent->_bestAnswer;
      _rightmostVertexIndex = parent->_rightmostVertexIndex;

    }

  }; // end class Coloring::Checker::Zykov
}// end namespace Coloring

