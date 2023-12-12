#include "TestSuite.h"
#include "UnionFind.h"
#include "Coloring.h"
#include "Graph.h"
#include "DIMACS.h"

int main() {
  TestSuite::run<int>("int");
  TestSuite::run<bool>("boolean");
  TestSuite::run<float>("float");

  Graph<int, int> g(2);
  g[0][1] = 0;
  g[1][1] = 0;
  g[1][1] = 0;
  g[1][0] = 0;

  Coloring::Checker<Graph<int,int>, int> checker(&g, 2);
  std::cout << checker.run() << '\n';

  try {
    DIMACS::Parser<Graph<int, int>> parser;
    Graph<int, int> g2 = parser.parse();
  }
  catch (const std::exception& e) {
    std::cerr << "CAUGHT: " << e.what() << std::endl;
  }
}