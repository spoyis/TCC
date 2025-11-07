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

#ifdef _WIN32
#include <windows.h>
#endif

std::mt19937 gen;
OutputWriter out;
unsigned long long GLOBAL_SEED = 0;

int main(int argc, char* argv[]) {
  TestSuite::run<int>("int");
  std::string inputFile;
  std::string outputFile;

  if (argc >= 3) {
    inputFile = argv[1];
    outputFile = argv[2];
  }
  else {

    bool underDebugger = false;
    #ifdef _WIN32
        underDebugger = (IsDebuggerPresent() != 0);
    #endif

    if (underDebugger) {
      // Defaults to use while debugging in Visual Studio
      inputFile = "2025-1.txt";
      outputFile = "2025-1-debug-output.txt";
      std::cout << "[DEBUG MODE] No args supplied - using defaults:\n"
        << "  input:  " << inputFile << "\n"
        << "  output: " << outputFile << "\n";
    }
    else {
      std::cerr << "Usage: " << argv[0] << " <inputfile> <outputfile>\n";
      return 1;
    }
  }


  std::random_device rd;
  GLOBAL_SEED = 1287245168;
  gen.seed(GLOBAL_SEED);
  out.seed = GLOBAL_SEED;

  std::cout << "Random seed: " << GLOBAL_SEED << "\n";


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
