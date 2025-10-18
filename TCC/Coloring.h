#pragma once
#include "Heuristics.h"
#include "Graph.h"
#include "BipartiteMatching.h"
#include <queue>
#include "validator.h"
#include <stdexcept>

long globalCounter = 0;

namespace Coloring {// begin namespace Coloring

  enum OptimizationStrategy {
    ROOT_NODE_SAME_TIME_DIFFERENT_ROOMS = 0,
    VERTICES_WITH_NO_COLOR_INTERSECTION = 1
  };

  template<class graph_t>
  class Checker {// begin class Coloring::Checker
    class Zykov;
    struct ZykovPtrComparator {
      bool operator()(Zykov* a, Zykov* b) const {
        return *a > *b;
      }
    };

  private:
    enum PruneMotive {
      not_pruned, bad_clique_contraction, no_valid_intermediate_coloring, no_valid_final_coloring
    };
  public:
    Checker(graph_t* graph, std::vector<std::vector<long>>& roomData) 
      : roomData(graph->getVertexCount()),
        optimizationStrategies( 2, false )
    {
      _originalGraph = graph;
      _initialVertexCount = _bestAnswer = _originalGraph->getVertexCount();

      for (long vertex = 0; vertex < _initialVertexCount; vertex++) {
        this->roomData[vertex] = std::move(roomData[vertex]);
      }
    }

    Checker(graph_t* graph, std::vector<std::vector<long>>& roomData, std::vector<OptimizationStrategy> strategies)
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

    // Move constructor
    Checker(Checker&& other) noexcept
      : _bestAnswer(other._bestAnswer),
      _initialVertexCount(other._initialVertexCount),
      _originalGraph(other._originalGraph),
      optimizationStrategies(std::move(other.optimizationStrategies)),
      roomData(std::move(other.roomData))
    {
      other._originalGraph = nullptr;
      other._bestAnswer = 0;
      other._initialVertexCount = 0;
    }

    int run() {
      if (Validator::naiveCheck<edge, vertex>(*_originalGraph, roomData) == Validator::MISSING_COLORS) return -1;

      if (isEnabled(ROOT_NODE_SAME_TIME_DIFFERENT_ROOMS)) {
        sameTimeDifferentRoomsOptimization();
      }
      std::priority_queue<Zykov*, std::vector<Zykov*>, ZykovPtrComparator> searchQueue;
      // Create root Zykov
      Zykov* root = new Zykov(_originalGraph->cloneHeap(), _initialVertexCount, this);
      searchQueue.push(root);
      long long step = 0;
      const int BATCH_SIZE = 10;

      while (!searchQueue.empty()) {
        if (hasSolution) {
          break;
        }
        // Extract batch of nodes (up to BATCH_SIZE)
        std::vector<Zykov*> batch;
        int batchCount = std::min(BATCH_SIZE, static_cast<int>(searchQueue.size()));

        for (int i = 0; i < batchCount; i++) {
          batch.push_back(searchQueue.top());
          searchQueue.pop();
        }

        // Process each node in the batch
        for (Zykov* current : batch) {
          if (step % 50 == 0 && _isRootChecker) {
            std::cout << "[CHECKER] STEP " << step++ << " LOWERBOUND == " << current->getLowerBound() << " UPPERBOUND == " << current->getUpperBound() << '\n';
            std::cout << "[CHECKER] DENSITY " << current->getGraphDensity() << '\n';
            std::cout << "[CHECKER] DEPTH " << current->getNodeDepth() << '\n';
            std::cout << "[CHECKER] BESTANSWER: " << _bestAnswer << "\n";
            std::cout << "[CHECKER] PRUNED NODE STATS:\n";
            std::cout << "[CHECKER][VALID NODES]: " << processedNodeStatistics[Validator::VALID] << '\n';
            std::cout << "[CHECKER][TIMESLOT][INVALID_JOINING]:" << processedNodeStatistics[Validator::INVALID_JOINING] << '\n';
            std::cout << "[CHECKER][TIMESLOT][PIGEONHOLE]: " << processedNodeStatistics[Validator::INVALID_TIMESLOT_PIGEONHOLE] << '\n';
            std::cout << "[CHECKER][TIMESLOT][NO_COLORING]: " << processedNodeStatistics[Validator::INVALID_TIMESLOT_COLORING] << '\n';
            std::cout << "[CHECKER][ROOM][PIGEONHOLE]: " << processedNodeStatistics[Validator::INVALID_ROOM_PIGEONHOLE] << '\n';
            std::cout << "[CHECKER][ROOM][NO_COLORING]: " << processedNodeStatistics[Validator::INVALID_ROOM_COLORING] << std::endl;
          }
          else step++;

          if (!_isRootChecker && step % 100 == 0) {
            std::cout << "[AUXILIARY CHECK] STEP: " << step << '\n';
          }
          

          extraValidationStepCounter++;
          nodesExplored++;
          // Pruning: if this node's bound is > current best, skip it
          if (current->getLowerBound() >= _bestAnswer){
            nodesPruned++;
            delete current;
            continue;
          }
          // Process current node
          auto result = current->processNode(getStrategy(), getEval(), searchQueue);
          processedNodeStatistics[result]++;

          if (shouldRunRecursiveCheck(result) && false) {
            // create new problem instance
            // the instance will be a subgraph of the original graph with the failing clique.
            Checker<graph_t> recursiveChecker = createRecursiveInstance(current);
            recursiveChecker.setStrategy(this->strategy);
            for (long i = 0; i < optimizationStrategies.size(); i++) {
              recursiveChecker.setOptimizationStrategy(static_cast<OptimizationStrategy>(i), optimizationStrategies[i]);
            }
            recursiveChecker.setRecursiveInstanceFlag();
            auto recursiveResult = recursiveChecker.run();
            if (!recursiveChecker.hasSolution ) {
              std::cout << "HOLY SHIT!!!!\n";
              std::cout << "HOLY SHIT!!!!\n";
              std::cout << "HOLY SHIT!!!!\n";
              std::cout << "HOLY SHIT!!!!\n";
              std::cout << "HOLY SHIT!!!!\n";
              std::cout << "HOLY SHIT!!!!\n";
              std::cout << "HOLY SHIT!!!!\n";
              std::cout << "HOLY SHIT!!!!\n";
              std::cout << "HOLY SHIT!!!!\n";
              break;
            }
            std::cout << "RECURSION RESULTED IN " << recursiveResult << '\n';
            extraValidationStepCounter = 0;
          }

          if (shouldRunBadCliqueCleanup(result)) {
            runBadCliqueCleanup(current, searchQueue);
            extraValidationStepCounter = 0;
          }

          delete current;
        }
      }
      // Clean up any remaining nodes in queue
      while (!searchQueue.empty()) {
        delete searchQueue.top();
        searchQueue.pop();
      }
      return _bestAnswer;
    }

    void setStrategy(Heuristic::Enums value) {
      strategy = value;
    }

    void setOptimizationStrategy(OptimizationStrategy strategy, bool value) {
      optimizationStrategies[strategy] = value;
    }

  private:
   
    using edge = typename graph_t::edge_t;
    using vertex = typename graph_t::vertex_t;
    bool _isRootChecker{ true };
    bool _areBadCliqueCleanupsEnabled{ true };
    bool _areRecursiveChecksEnabled{ true };
    long _bestAnswer;
    long _initialVertexCount;
    long extraValidationStepCounter = { 0 };
    bool hasSolution{ false };
    Graph<edge, vertex>* _originalGraph;
    Heuristic::Enums strategy{ Heuristic::STRATEGY_DEFAULT };
    Heuristic::Enums eval{ Heuristic::EVAL_DEFAULT };
    std::vector<bool> optimizationStrategies{false, false};
    std::vector<std::vector<long>> roomData;
    std::vector<long> processedNodeStatistics{0,0,0,0,0,0};

    // Statistics
    long nodesExplored = 0;
    long nodesPruned = 0;

    Heuristic::strategy<edge,vertex> getStrategy() {
      switch (strategy) {
      case Heuristic::STRATEGY_FIRST_VALID:
        return &Heuristic::firstValid<edge, vertex>;
      
      case Heuristic::STRATEGY_LOWEST_DEGREE:
        return &Heuristic::lowestDegree<edge, vertex>;

      case Heuristic::STRATEGY_RANDOM_FROM_CLIQUE:
        return &Heuristic::randomFromClique<edge, vertex>;

      case Heuristic::STRATEGY_SHARED_NEIGHTBORS:
        return &Heuristic::mostSharedNeighbors<edge, vertex>;
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

    bool shouldRunRecursiveCheck(long result) {
      return result != Validator::VALID && extraValidationStepCounter >= 100 && _isRootChecker;
    }

    bool shouldRunBadCliqueCleanup(long result)  {
      return (shouldRunRecursiveCheck(result)); // same logic, for now - change if necessary in the future.
    }

    bool shouldZykovNodesSaveBadCliques() {
      return (_areRecursiveChecksEnabled || _areBadCliqueCleanupsEnabled) && extraValidationStepCounter >= 100 && _isRootChecker;
    }

    void sameTimeDifferentRoomsOptimization() {
      std::vector<std::pair<long, long>> nonNeighbors;
      Graph<edge, vertex>& graph = *_originalGraph;

      for (long i = 0; i < _initialVertexCount; i++)
      {
        for (long j = 0; j < _initialVertexCount; j++) {
          if (i == j) continue;
          if (roomData[i].size() > 1 || roomData[j].size() > 1) continue;

          if (!graph[i][j]) nonNeighbors.push_back({ i,j });
        }
      }

      for (auto& vertexPair : nonNeighbors) {
        auto v1 = vertexPair.first;
        auto v2 = vertexPair.second;

        if (roomData[v1][0] == roomData[v2][0]) {
          graph[v1][v2] = {1};
          graph[v2][v1] = {1};

          std::cout << "[OPTIMIZING] CLASS " << graph.getVertexLabel(v1) << " AND " << graph.getVertexLabel(v2) << " SHARE THE SAME ONE ROOM\n";
          //std::cout << "[DEBUG] Added edge between " << v1 << " and " << v2 << std::endl;
        }
      }
    }

    Checker<graph_t> createRecursiveInstance(Zykov* node) {
      std::vector<long> vertices = node->getJoinedVerticesFromClique(node->getBadClique());

      std::vector<std::vector<long>>* _newRoomData = new std::vector<std::vector<long>>(vertices.size());
      Graph<edge, vertex>* _newGraph = new Graph<edge,vertex>(vertices.size());
      Graph<edge, vertex>& newGraph = *_newGraph;
      Graph<edge, vertex>& originalGraph = *_originalGraph;

      for (long i = 0; i < vertices.size(); ++i) {
        long oldI = vertices[i];
        //set vertex data (timeslots)
        auto vertexData = originalGraph.getVertexData(oldI);
        newGraph.setVertexData(i, vertexData);
        // set room data
        (*_newRoomData)[i] = roomData[oldI];
        // set labels
        newGraph.setVertexLabel(node->getVertexLabel(oldI), i);

        // copy adjacency
        for (long j = 0; j < vertices.size(); ++j) {
          long oldJ = vertices[j];
          newGraph[i][j] = originalGraph.at(oldI, oldJ);
        }
      }

      return Checker<graph_t>(_newGraph, *_newRoomData);
    }

    void setRecursiveInstanceFlag() {
      _isRootChecker = false;
    }

    void runBadCliqueCleanup(Zykov* badNode, std::priority_queue<Zykov*, std::vector<Zykov*>, ZykovPtrComparator>& searchQueue) {
      std::vector<long> badClique = badNode->getBadClique();
      Graph<edge, vertex>& badGraph = *badNode->getGraph();

      bool reduced = true;
      long reductionCounter = 0;
      while (reduced && badClique.size() > 1) {
        reduced = false;
        for (long i = 0; i < badClique.size(); ++i) {
          std::vector<long> subClique;
          for (long j = 0; j < badClique.size(); ++j) {
            if (j != i) subClique.push_back(badClique[j]);
          }
          auto validation = Validator::clique(subClique, badGraph, roomData);
          if (validation != Validator::VALID) {
            badClique = subClique;
            reduced = true;
            reductionCounter++;
            break; // restart from the smaller clique
          }
        }
      }
      std::cout << "[CLEANUP] REDUCED BAD CLIQUE BY  " << reductionCounter << '\n' << "[CLEANUP] BAD CLIQUE SIZE IS " << badClique.size() << '\n';

      std::vector<long> badVertices = badNode->getJoinedVerticesFromClique(badClique);

      std::vector<Zykov*> extractQueue;
      while (!searchQueue.empty()) {
        extractQueue.push_back(searchQueue.top());
        searchQueue.pop();
      }

      std::vector<Zykov*>duplicateQueue;

      for (Zykov* currentNode : extractQueue) {

        Graph<edge, vertex>& currentGraph = *currentNode->getGraph();
        bool sameClique = true;

        // 1. Check that all roots in badClique are still roots in this node
        for (auto rootVertex : badClique) {
          if (currentGraph.getRoot(rootVertex) != rootVertex) {
            sameClique = false;
            break;
          }
        }

        if (!sameClique) {
          duplicateQueue.push_back(currentNode);
          continue;
        }

        // check 2: all non-root vertices in badVertices have same root in both
        for (long v : badVertices) {
          
          long rBad = badGraph.getRoot(v);
          long rCur = currentGraph.getRoot(v);
          if (rBad != rCur) {
            sameClique = false;
            break;
          }
        }

        if (!sameClique) {
          duplicateQueue.push_back(currentNode);
          continue;
        }

        // check 3: the roots in badClique form the same clique in current
        for (long i = 0; i < badClique.size() && sameClique; ++i) {
          for (long j = i + 1; j < badClique.size(); ++j) {
            long a = badClique[i];
            long b = badClique[j];
            if (!currentGraph[a][b]) {
              sameClique = false;
              break;
            }
          }
        }

        if (!sameClique) {
          duplicateQueue.push_back(currentNode); // keep it — it’s different
        }
        else // if sameClique == true, skip reinsertion, clique is invalid.
          delete currentNode; 
      }

      for (Zykov* node : duplicateQueue) {
        searchQueue.push(node);
      }


      std::cout << "[CLEANUP] Removed "
        << (extractQueue.size() - duplicateQueue.size())
        << " nodes with duplicate bad clique\n";
    }

  }; // end class Coloring::Checker


  template<class graph_t>
  class Checker<graph_t>::Zykov  
  { // begin class Coloring::Checker::Zykov

  private:
    long _currentVertexCount;
    long _rightmostVertexIndex;
    long _bestAnswer;
    std::vector<long> badClique;
    Graph<edge, vertex>* _graph;
    Checker* _checkerPtr;
    bool _isValidNode{true};
    PruneMotive pruneMotive{ not_pruned };
    double _upperBound;
    double _lowerBound{ 0 };
    double _graphDensity;
    long _depth{ 0 };
    std::string _lastOperation{ "none" }; // debug helper!
  public:

    // for branch and bound comparison.
    bool operator>(const Zykov& other) const {
      //return this->_currentVertexCount > other._currentVertexCount; // vertex count
      //return this->getGraphDensity() < other.getGraphDensity();   // density
      return (this->getUpperBound() - this->getLowerBound()) > (other.getUpperBound() - other.getLowerBound()); // bound differential.
    }

    Zykov(Graph<edge,vertex>* graph, int size, Checker* ptr) {
      _graph = graph;
      _currentVertexCount = _bestAnswer = size;
      _rightmostVertexIndex = size - 1;
      _checkerPtr = ptr;
      computeBounds(ptr->isEnabled(VERTICES_WITH_NO_COLOR_INTERSECTION));
    }

    long getLowerBound() const { return _lowerBound;}
    long getUpperBound() const { return _upperBound;}
    std::vector<long> getBadClique() { if (badClique.size() == 0) std::cout << "[THIS SHOULDNT HAPPEN] asked for an empty \"bad\" clique ?!\n"; return badClique; }
    double getGraphDensity() const { return _graphDensity; }
    long getNodeDepth() const { return _depth; }
    Graph<edge, vertex>* getGraph() const { return _graph; }
    std::string getVertexLabel(long index) const { return _graph->getVertexLabel(index); }

    std::vector<long> getJoinedVerticesFromClique(const  std::vector<long>& clique) {
      std::vector<long> output(clique);
      // all vertices in clique are in the original graph, we need to find the ones we "missed" when we joined them.
      Graph<edge, vertex>& graph = *_graph;

      for (auto& v : clique) {
        if (graph.getRoot(v) != v) {
          throw std::logic_error("Clique contains a non-root vertex");
        }
        for (long i = 0; i < _graph->getVertexCount(); i++) {
          auto root = graph.getRoot(i);
          if (root != i && root == v) {
            output.push_back(i);
          }
        }
      }

      return output;
    }

    ~Zykov() { delete _graph; }

    void edgeAdditionOptimization(const std::vector<long>& clique) {
      Graph<edge, vertex>& graph = *_graph;
      long n = graph.getVertexCount();
      int edgesAdded = 0;

      // Iterate over all pairs of vertices outside the clique
      for (long x = 0; x < n; ++x) {
        auto rootX = graph.getRoot(x);
        if (std::find(clique.begin(), clique.end(), rootX) != clique.end()) continue; // skip if in clique
        if (rootX != x) continue; // skip if not a root

        for (long y = x + 1; y < n; ++y) {
          auto rootY = graph.getRoot(y);
          if (std::find(clique.begin(), clique.end(), rootY) != clique.end()) continue;
          if (rootY != y) continue; // skip if not a root
          if (graph[x][y]) continue; // skip existing edges

          // Check if every vertex in clique is adjacent to at least one of x or y
          bool allCovered = true;
          for (auto z : clique) {
            if (!graph[x][z] && !graph[y][z]) {
              allCovered = false;
              break;
            }
          }

          if (allCovered) {
            graph[x][y] = edge{ 1 };
            graph[y][x] = edge{ 1 };
            edgesAdded++;
          }
        }
      }
    }

    
    long processNode(
      Heuristic::strategy<edge, vertex> heuristic,
      Heuristic::eval<edge> eval,
      std::priority_queue<Zykov*, std::vector<Zykov*>, ZykovPtrComparator>& searchQueue) {
      
      // Validate multiple random cliques for more aggressive pruning
      const int NUM_CLIQUES_TO_VALIDATE = 5;
      std::vector<long> largestClique;
      bool saveBadCliques = _checkerPtr->shouldZykovNodesSaveBadCliques();

      for (int i = 0; i < NUM_CLIQUES_TO_VALIDATE; ++i) {
        auto clique = (*_graph).findCliqueGrasp();
        auto cliqueValidation = Validator::clique<edge, vertex>(clique, *_graph, _checkerPtr->roomData);

        if (cliqueValidation != Validator::VALID) {
          // If any clique is invalid, we can prune this node
          if (saveBadCliques) {
            badClique = clique;
          }
          return cliqueValidation;
        }
        
        if (clique.size() > largestClique.size()) {
          largestClique = clique;
        }
      }
      _lowerBound = std::max(_lowerBound, (double)largestClique.size());
      edgeAdditionOptimization(largestClique);
      
      
      // this returns {-1, -1} if no valid vertex pair is found
      auto nonNeighboringVertices = heuristic(*_graph, eval, largestClique, _rightmostVertexIndex, _currentVertexCount);
      //std::cout << "SELECTED " << nonNeighboringVertices.first << " AND " << nonNeighboringVertices.second << '\n';
      // If no valid vertex pair found, check if complete
      if (nonNeighboringVertices.first == -1) {
        if (_graph->isUndirectedComplete()) {
          _checkerPtr->hasSolution = true;
          std::cout << "GRAFO COMPLETO\n";
          if (largestClique.size() < _checkerPtr->_bestAnswer)
            _checkerPtr->_bestAnswer = largestClique.size();
        }
        else {
          std::cout << "[ERROR] THIS ISN'T SUPPOSED TO HAPPEN EVER!!!\n";
        }
        return Validator::VALID;
      }

      //std::cout << "FOUND " << nonNeighboringVertices.first << " AND " << nonNeighboringVertices.second << '\n';

      // Child 1: Add edge
      Zykov* addEdgeChild = new Zykov(*_graph, nonNeighboringVertices, edge{ 1 }, this);
      auto addEdgeClique = addEdgeChild->_graph->findCliqueRandom({ nonNeighboringVertices.first, nonNeighboringVertices.second });
      auto addEdgeValidation = Validator::clique<edge, vertex>(addEdgeClique, *addEdgeChild->_graph, _checkerPtr->roomData);
      if (addEdgeValidation == Validator::VALID) {
        searchQueue.push(addEdgeChild);
      }
      else {
        _checkerPtr->processedNodeStatistics[addEdgeValidation]++;
        delete addEdgeChild;
      }

      if ((*_graph).areJoinable(nonNeighboringVertices.first, nonNeighboringVertices.second)) {
        // Child 2: Contract vertices  
        Zykov* contractChild = new Zykov(*_graph, nonNeighboringVertices, this);
        Graph<edge, vertex>& contractGraph = *contractChild->_graph;
        auto contractClique = contractGraph.findCliqueRandom({ contractGraph.getRoot(nonNeighboringVertices.first) });
        auto contractValidation = Validator::clique<edge, vertex>(contractClique, *contractChild->_graph, _checkerPtr->roomData);

        if (contractValidation == Validator::VALID) {
          searchQueue.push(contractChild);
        }
        else {
          _checkerPtr->processedNodeStatistics[contractValidation]++;
          delete contractChild;
        }
      }
      else _checkerPtr->processedNodeStatistics[Validator::INVALID_JOINING]++;

      return Validator::VALID;
    }

  private:

    void computeBounds(bool searchForNonJoinableVertices = false) {
      std::vector<std::pair<long, long>> nonNeighbors;
      Graph<edge, vertex>& graph = *_graph;
      long maxDegree = 0;
      long vertexCount = graph.getVertexCount();
      long bestVertex = 0;
      for (long i = 0; i < vertexCount; i++) {
        long rootI = graph.getRoot(i);
        if (rootI < i) 
          continue;
        long degree = 0;
        for (long j = 0; j < vertexCount; j++) {
          long rootJ = graph.getRoot(j);
          if (rootJ < j) continue;
          if (i != j && graph[rootI][rootJ]) degree++;
          else if (i != j && !graph[rootI][rootJ] && searchForNonJoinableVertices)
            nonNeighbors.push_back({ i,j });
        }
        if (degree > maxDegree) {
          maxDegree = degree;
          bestVertex = i;
        }
      }
      _upperBound = maxDegree + 1; // Degree-based upper bound

      // SET GRAPHDENSITY
      double V = (double)_currentVertexCount;
      double E = (double)(*_graph).getEdgeCount();

      _graphDensity = E / (V * (V - 1));


      if (searchForNonJoinableVertices) {
        long optimizationCounter = 0;
        for (auto& vertexPair : nonNeighbors) {
          auto v1 = vertexPair.first;
          auto v2 = vertexPair.second;

          // if can't be joined, place an edge.
          if (!graph.areJoinable(v1, v2)) {
            graph[v1][v2] = { 1 };
            graph[v2][v1] = { 1 };

            optimizationCounter++;
          }
        }
        if(optimizationCounter != 0)
        std::cout << "[OPTIMIZATION][NON_JOINABLE] OPTIMIZED " << optimizationCounter << " EDGES!\n";
      }
    }

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
      _lowerBound = parent->_lowerBound;
      _depth = parent->_depth + 1;
      computeBounds();
     
      _lastOperation = "added edge between " + std::to_string(vertices.first) + " and " + std::to_string(vertices.second);
    }

    // contract vertices constructor
    Zykov(Graph<edge, vertex>& g, std::pair<long, long> vertices, Zykov* parent) {
      _graph = g.cloneHeap();

      Graph<edge, vertex>& graph = *_graph;
      auto& roomData = parent->_checkerPtr->roomData;

      if (graph[vertices.first][vertices.second] || graph[vertices.second][vertices.first]) {
        std::cout << "BAD CONTRACTION!!!!! VERTICES ARE NEIGHBORS\n";
      }

      if (graph.getRoot(vertices.first) != graph.getRoot(vertices.first) || graph.getRoot(vertices.second) != graph.getRoot(vertices.second)) {
        std::cout << "BAD CONTRACTION!!!!! VERTICES ARE NOT ROOTS\n";
      }
      
      graph.joinVertices(vertices.first, vertices.second);
      _checkerPtr = parent->_checkerPtr;
      _currentVertexCount = parent->_currentVertexCount - 1;
      _bestAnswer = parent->_bestAnswer;
      _rightmostVertexIndex = parent->_rightmostVertexIndex;
      _lowerBound = parent->_lowerBound;
      _depth = parent->_depth + 1;
      computeBounds();

      _lastOperation = "joined " + std::to_string(vertices.first) + " and " + std::to_string(vertices.second);
    }

  }; // end class Coloring::Checker::Zykov
}// end namespace Coloring

