#include "TestSuite.h"
#include "UnionFind.h"
#include "Coloring.h"
#include "Graph.h"
#include "DIMACS.h"
#include "BipartiteMatching.h"
#include "InputParser.h"
#include "OutputWriter.h"
#include <chrono>
#include <iostream>

std::mt19937 gen;
OutputWriter out;
unsigned long long GLOBAL_SEED = 0;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <inputfile> <outputfile>\n";
    return 1;
  }

  std::random_device rd;
  GLOBAL_SEED = rd();
  gen.seed(GLOBAL_SEED);
  out.seed = GLOBAL_SEED;

  std::cout << "Random seed: " << GLOBAL_SEED << "\n";

  std::string inputFile = argv[1];
  std::string outputFile = argv[2];

  try {
    std::cout << "BEGIN PARSING INPUT: " << inputFile << "\n";
    auto parser = input::Parser(inputFile);
    parser.parse();
    std::cout << "INPUT SUCCESSFULLY PARSED\n";

    std::cout << "PROBLEM CHECKER BEGIN PROCESSING\n";
    auto checker = parser.getChecker();

    

    checker.setStrategy(Coloring::Heuristic::STRATEGY_SHARED_NEIGHTBORS);
    checker.setOptimizationStrategy(Coloring::ROOT_NODE_SAME_TIME_DIFFERENT_ROOMS, true);
    checker.setOptimizationStrategy(Coloring::VERTICES_WITH_NO_COLOR_INTERSECTION, true);
    checker.setOptimizationStrategy(Coloring::ROOT_NODE_SINGLE_CONSTRAINT_PROPAGATION, true);

    auto start = std::chrono::high_resolution_clock::now();
    auto result = checker.run();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    out.executionTime = elapsed;

    std::cout << "PROBLEM CHECKER FINISHED PROCESSING\n";
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
   
    out.writeToFile(outputFile);
    std::cout << "Output written to: " << outputFile << "\n";
  }
  catch (const std::exception& e) {
    std::cerr << "\nCAUGHT: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
