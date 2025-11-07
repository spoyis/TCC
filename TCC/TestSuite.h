#pragma once
#include <iostream>
#include "Graph.h"
#include <string>
#include "Coloring.h"
#include "InputParser.h"

namespace TestSuite { // begin namespace TestSuite

  template<int> auto constexpr defaultVal() { return 0; }

  namespace color {
    const std::string red("\033[0;31m");
    const std::string green("\033[1;32m");
    const std::string yellow("\033[1;33m");
    const std::string blue("\033[1;34m");
    const std::string cyan("\033[0;36m");
    const std::string magenta("\033[0;35m");
    const std::string reset("\033[0m");
  }
  
  static int totalFailss = 0;
  static int totalPasses = 0;

  static const std::string PASS = "- " + color::green + "[PASSED]: " + color::reset;
  static const std::string FAIL = "- " + color::red + "[FAILED]: " + color::reset;

   template<typename T>
  void run(const char* type)
  {
    std::cout << color::blue << "STARTING TEST SUITE FOR: ["
      << color::reset + color::magenta << type
      << color::reset + color::blue << ']'
      << color::reset << std::endl;

    testCompleteGraphValidation<T>(5);
    testCompleteGraphValidation<T>(10);
    testGraphCloning<T>();
    testEdgeJoining<T>();
    testJoinEdgeConsistency<T>();
    testVertexDegreeConsistency<T>();
  }


  

  template<typename T>
  void expect(const char* outputText, const T& expectedValue, const T& receivedValue) {
    if (receivedValue == expectedValue) {
      std::cout << PASS << outputText << std::endl;
      totalPasses++;
    }
    else {
      std::cout << FAIL << outputText
        << " (Expected: " << expectedValue << ", Got: " << receivedValue << ")"
        << std::endl;
      totalFailss++;
    }
  }


  void runChecker() {
    std::cout << color::blue << "STARTING CHECKER TESTS\n" << color::reset;

    {
      std::cout << color::cyan << "  -> TRIANGLE GRAPH 1:\n" << color::reset;
      auto parser = input::Parser("test_1.txt");
      parser.parse();

      auto checker = parser.getChecker();
      checker.setStrategy(Coloring::Heuristic::STRATEGY_LOWEST_DEGREE);
      //int result = checker.run();
      // rewrite this later.
      //expect("Test 1 should return a minimum coloring of 3", 3, result);
    }

    {
      std::cout << color::cyan << "  -> TEST 2:\n" << color::reset;
      auto parser = input::Parser("test_2.txt");
      parser.parse();

      auto checker = parser.getChecker();
      checker.setStrategy(Coloring::Heuristic::STRATEGY_LOWEST_DEGREE);
      //int result = checker.run();
      // rewrite this later.
      //expect("Test 2 should return a minimum coloring of 2", 2, result);
    }

    {
      std::cout << color::cyan << "  -> TEST 3:\n" << color::reset;
      auto parser = input::Parser("test_3.txt");
      parser.parse();

      auto checker = parser.getChecker();
      checker.setStrategy(Coloring::Heuristic::STRATEGY_LOWEST_DEGREE);
      //int result = checker.run();
      // rewrite this later.
      //expect("Test 3 should return a minimum coloring of 2", 2, result);
    }

    std::cout << color::blue << "CHECKER TESTS COMPLETE\n" << color::reset;
  }


  template<typename T>
  Graph<T, T> generateCompleteGraph(int size) {
    Graph<T, T> output(size);

    for (long i = 0; i < size; i++)
      for (long j = 0; j < size; j++)
        output[i][j] = truthyValue<T>();

    return output;
  }

  template<typename T>
  void testCompleteGraphValidation(int size) {
    std::cout << color::yellow << "--> [COMPLETE GRAPH SIZE " << size << "]\n" << color::reset;
    Graph<T, T> G = generateCompleteGraph<T>(size);

    expect("It should correctly validate complete graphs as true ", true, G.isDirectedComplete());
  }

  template<typename T>
  void testGraphCloning() {
    std::cout << color::yellow << "--> [ GRAPH DEEP COPY " << "]\n" << color::reset;
    
    Graph<T, T> G1 = generateCompleteGraph<T>(5);
    auto G2 = G1.clone();

    expect("It should correctly perform deep copies when called ", true, G1 == G2);
  }

  template<typename T>
  void testEdgeJoining() {
    std::cout << color::yellow << "--> [EDGE JOINING VALIDATION]\n" << color::reset;

    // Create 4-vertex undirected simple graph
    Graph<T, T> G(4);

    // Add symmetric edges: 0–1, 1–2, 2–3
    G[0][1] = truthyValue<T>();
    G[1][0] = truthyValue<T>();
    G[1][2] = truthyValue<T>();
    G[2][1] = truthyValue<T>();
    G[2][3] = truthyValue<T>();
    G[3][2] = truthyValue<T>();

    long beforeEdges = G.getEdgeCount();
    G.joinVertices(1, 2);

    // Determine surviving root
    int root1 = G.getRoot(1);
    int root2 = G.getRoot(2);
    expect("After join, roots should be identical", root1, root2);

    // Helpers
    auto isRoot = [&](int v) { return G.getRoot(v) == v; };

    // 1) Symmetry among roots only
    bool symmetricRoots = true;
    for (int i = 0; i < 4; ++i) {
      if (!isRoot(i)) continue;
      for (int j = 0; j < 4; ++j) {
        if (!isRoot(j)) continue;
        if (G[i][j] != G[j][i]) {
          symmetricRoots = false;
        }
      }
    }
    expect("Graph should remain symmetric among roots", true, symmetricRoots);

    // 2) No self-loops on roots
    bool noSelfLoops = true;
    for (int i = 0; i < 4; ++i) {
      if (!isRoot(i)) continue;
      if (G[i][i]) noSelfLoops = false;
    }
    expect("No self-loop should remain on any root", true, noSelfLoops);

    // 3) Non-root rows/cols should be cleared
    bool nonRootsCleared = true;
    for (int i = 0; i < 4; ++i) {
      if (isRoot(i)) continue;
      for (int j = 0; j < 4; ++j) {
        if (G[i][j] != 0 || G[j][i] != 0) {
          nonRootsCleared = false;
        }
      }
    }
    expect("Non-root rows/cols should be zeroed", true, nonRootsCleared);
  }

  template<typename T>
  void testJoinEdgeConsistency() {
    std::cout << color::yellow << "--> [CONSECUTIVE EDGE JOIN VALIDATION]\n" << color::reset;

    const int N = 6;
    Graph<T, T> G(N);

    // complete graph with no self loops
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        if (i == j) continue;
        G[i][j] = 1;
        G[j][i] = 1;
      }
    }

    auto checkInvariants = [&](const char* stage) {
      bool symmetricRoots = true;
      bool noSelfLoops = true;
      bool nonRootsCleared = true;
      bool rootsDontTouchNonRoots = true;
      bool edgeCountCorrect = true;

      // Count roots
      std::vector<int> roots;
      for (int i = 0; i < N; ++i)
        if (G.getRoot(i) == i)
          roots.push_back(i);

      int rootCount = static_cast<int>(roots.size());

      // Compute actual edge count among roots
      long actualEdges = 0;
      for (int i = 0; i < N; ++i) {
        bool isRootI = (G.getRoot(i) == i);

        // (1) symmetry among roots and count edges
        if (isRootI) {
          for (int j = i + 1; j < N; ++j) {
            bool isRootJ = (G.getRoot(j) == j);
            if (!isRootJ) continue;
            if (G[i][j] != G[j][i]) symmetricRoots = false;
            if (G[i][j]) actualEdges++;
          }
        }

        // (2) no selfloops on roots
        if (isRootI && G[i][i]) noSelfLoops = false;

        // (3) nonroot rows/cols zeroed
        if (!isRootI) {
          for (int j = 0; j < N; ++j) {
            if (G[i][j] != 0 || G[j][i] != 0) nonRootsCleared = false;
          }
        }

        // (4) roots don't connect to nonroots
        if (isRootI) {
          for (int j = 0; j < N; ++j) {
            if (G.getRoot(j) != j && (G[i][j] != 0 || G[j][i] != 0))
              rootsDontTouchNonRoots = false;
          }
        }
      }

      // (5) theoretical edge count
      long expectedEdges = (rootCount * (rootCount - 1)) / 2;
      if (actualEdges != expectedEdges) edgeCountCorrect = false;

      std::cout << color::cyan << "  -> " << stage << color::reset << "\n";
      expect("Symmetry among roots", true, symmetricRoots);
      expect("No self-loops on roots", true, noSelfLoops);
      expect("Non-root rows/cols cleared", true, nonRootsCleared);
      expect("Roots don't connect to non-roots", true, rootsDontTouchNonRoots);
      expect("Edge count matches K_n formula", true, edgeCountCorrect);

      if (!edgeCountCorrect) {
        std::cout << color::yellow
          << "    [DEBUG] Expected " << expectedEdges
          << " edges for K_" << rootCount
          << ", but found " << actualEdges << color::reset << "\n";
      }
      };


    checkInvariants("Initial (baseline)");

    // Perform consecutive joins that chain
    G.joinVertices(0, 1);
    checkInvariants("After join(0,1)");

    G.joinVertices(0, 2);
    checkInvariants("After join(0,2)");

    G.joinVertices(3, 4);
    checkInvariants("After join(3,4)");

    G.joinVertices(4, 5);
    checkInvariants("After join(4,5)");

    G.joinVertices(0, 3);
    checkInvariants("After join(0,3) - merge of two contracted components");
  }


  template<typename T>
  void testVertexDegreeConsistency() {
    std::cout << color::yellow << "--> [VERTEX DEGREE CONSISTENCY]\n" << color::reset;

    const int N = 6;
    Graph<T, T> G(N);

    // Build complete graph with no self-loops
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        if (i != j) G[i][j] = 1;
      }
    }

    auto checkDegrees = [&](const char* stage) {
      std::vector<int> roots;
      for (int i = 0; i < N; ++i)
        if (G.getRoot(i) == i)
          roots.push_back(i);

      int rootCount = static_cast<int>(roots.size());
      bool allMatch = true;

      for (auto r : roots) {
        long deg = G.getVertexDegree(r);
        long expected = rootCount - 1;
        if (deg != expected) {
          allMatch = false;
          std::cout << color::yellow << "    [DEBUG] root " << r
            << " degree " << deg << " (expected " << expected << ")\n"
            << color::reset;
        }
      }

      expect(stage, true, allMatch);
      };

    checkDegrees("Initial (K6)");
    G.joinVertices(0, 1);
    checkDegrees("After join(0,1)");
    G.joinVertices(0, 2);
    checkDegrees("After join(0,2)");
    G.joinVertices(3, 4);
    checkDegrees("After join(3,4)");
    G.joinVertices(0, 3);
    checkDegrees("After join(0,3)");
  }
}// end namespace TestSuite
