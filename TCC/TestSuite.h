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
      int result = checker.run();

      expect("Test 1 should return a minimum coloring of 3", 3, result);
    }

    {
      std::cout << color::cyan << "  -> TEST 2:\n" << color::reset;
      auto parser = input::Parser("test_2.txt");
      parser.parse();

      auto checker = parser.getChecker();
      checker.setStrategy(Coloring::Heuristic::STRATEGY_LOWEST_DEGREE);
      int result = checker.run();

      expect("Test 2 should return a minimum coloring of 2", 2, result);
    }

    {
      std::cout << color::cyan << "  -> TEST 3:\n" << color::reset;
      auto parser = input::Parser("test_3.txt");
      parser.parse();

      auto checker = parser.getChecker();
      checker.setStrategy(Coloring::Heuristic::STRATEGY_LOWEST_DEGREE);
      int result = checker.run();

      expect("Test 3 should return a minimum coloring of 2", 2, result);
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
}// end namespace TestSuite
