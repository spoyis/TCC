#include "TestSuite.h"
#include "UnionFind.h"
#include "Coloring.h"
#include "Graph.h"
#include "DIMACS.h"
#include "BipartiteMatching.h"
#include "InputParser.h"


int main() {
  TestSuite::run<int>("int");
  TestSuite::run<bool>("boolean");
  TestSuite::run<float>("float");
  //TestSuite::runChecker();
  
  
 // try {
 //   DIMACS::Parser<Graph<int, int>> parser;
  //  Graph<int, int> g2 = parser.parse();
 // }
 // catch (const std::exception& e) {
  //  std::cerr << "CAUGHT: " << e.what() << std::endl;
 // }
  /**/

  try {
   

    std::cout << "BEGIN PARSING INPUT\n";
    auto parser = input::Parser(std::string("base.txt"));
    parser.parse();
    std::cout << "INPUT SUCCESSFULLY PARSED\n";

    std::cout << "PROBLEM CHECKER BEGIN PROCESSING\n";
    auto checker = parser.getChecker();
    //checker.setStrategy(Coloring::Heuristic::STRATEGY_LOWEST_DEGREE);
    checker.setStrategy(Coloring::Heuristic::STRATEGY_RANDOM_FROM_CLIQUE);
    checker.setOptimizationStrategy(Coloring::ROOT_NODE_SAME_TIME_DIFFERENT_ROOMS, true);
    checker.setOptimizationStrategy(Coloring::VERTICES_WITH_NO_COLOR_INTERSECTION, true);
    auto result = checker.run();
    std::cout << "PROBLEM CHECKER FINISHED PROCESSING --- MIN COLORING VALUE IS " << result << "\n";
  }
  catch (const std::exception& e) {
    std::cerr << "\nCAUGHT: " << e.what() << std::endl;
  }
}

