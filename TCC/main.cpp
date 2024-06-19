#include "TestSuite.h"
#include "UnionFind.h"
#include "Coloring.h"
#include "Graph.h"
#include "DIMACS.h"
#include "BipartiteMatching.h"

int main() {
  TestSuite::run<int>("int");
  TestSuite::run<bool>("boolean");
  TestSuite::run<float>("float");

  
  Graph<int, int> g(2);
  g[0][1] = 0;
  g[1][1] = 0;
  g[1][1] = 1;
  g[1][0] = 1;
  std::vector<std::vector<long>> lol(2);
  Coloring::Checker<Graph<int,int>, int> checker(&g, lol);
  std::cout << checker.run() << '\n';


  
 // try {
 //   DIMACS::Parser<Graph<int, int>> parser;
  //  Graph<int, int> g2 = parser.parse();
 // }
 // catch (const std::exception& e) {
  //  std::cerr << "CAUGHT: " << e.what() << std::endl;
 // }
  /**/

  long n,m, t;
  std::cin >> t;
  std::cin >> n;
  
  std::vector<long> conjuntoN;
  std::vector<long> conjuntoM;

  for (long i = 0; i < n; i++)
  {
    long buffer;
    std::cin >> buffer;
    conjuntoN.push_back(buffer);
  }

  std::cin >> m;
  BipartiteMatching teste(n, m);
  for (long i = 0; i < m; i++) {
    long buffer;
    std::cin >> buffer;
    conjuntoM.push_back(buffer);
  }

  // cruzar os dois 

  for (long i = 0; i < n; i++) {
    teste.defineInitialCapacity(i, 1);
    for (long j = 0; j < m; j++) {
      if (conjuntoM[j] % conjuntoN[i] == 0) {
        teste.addEdge(i, j);
        std::cout << conjuntoN[i] << "  " << conjuntoM[j] << '\n';
      }
    }
  }

  auto resultado = teste.solve();
  std::cout << resultado << std::endl;
}

