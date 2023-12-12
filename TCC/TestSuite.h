#pragma once
#include <iostream>
#include "Graph.h"
#include <string>

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
  

  void expect(const char* outputText, bool expectedValue, bool receivedValue) {
    if (receivedValue == expectedValue) {
      std::cout << PASS << outputText << std::endl;
    }
    else {
      std::cout << FAIL << outputText << std::endl;
    }
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
