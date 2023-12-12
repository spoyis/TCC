#pragma once
#include "Heuristics.h"
#include "Graph.h"

namespace Coloring {// begin namespace Coloring


  template<class graph_t, typename color>
  class Checker {// begin class Coloring::Checker
  public:
    Checker(graph_t* graph, int size) {
      _originalGraph = graph;
      _initialVertexCount = _bestAnswer = size;
    }

    int run() {
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

  }; // end class Coloring::Checker


  template<class graph_t, typename color>
  class Checker<graph_t, color>::Zykov  
  { // begin class Coloring::Checker::Zykov
  private:
    long _currentVertexCount;
    long _rightmostVertexIndex;
    long _bestAnswer;
    Graph<edge, vertex>* _graph;
    Checker* _checkerPtr;
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
