#pragma once
#include "Coloring.h"
#include "Graph.h"
#include "BipartiteMatching.h"

namespace Coloring::Validator { // begin namespace Coloring::Validator


  //TODO: ENUM FOR DIFFERENT INVALID SITUATIONS
  using clique_t = std::vector<long>;

  // validates a clique
  template<typename edge, typename vertex>
  bool clique(clique_t clique, Graph<edge,vertex>& g) {
    vertex colors;
    
    for (auto& v : clique) {
      auto& colorArray = g.getVertexData(v);

      for (auto& color : colorArray) {
        bool flag = true;
        for (auto& processedColor : colors) {
          if (color == processedColor) {
            flag = false;
            break;
          }
        }
        if (flag) colors.push_back(color);
      }
    }

    // pigeonhole principle
    if (colors.size() < clique.size()) {


      /*std::cout << "====PIGEONHOLED====\n";

      for (auto v : clique) {
        std::cout << "PIGEON V: " << v << " LABEL === " << g.getVertexLabel(v) << '\n';
      }*/
      return false;

    }
 
    auto timeslotValidation = BipartiteMatching(clique.size(), colors.size());

    for (long v = 0; v < clique.size(); v++) {
      auto& colorArray = g.getVertexData(clique[v]);
      timeslotValidation.defineInitialCapacity(v, 1);
      for (long c = 0; c < colorArray.size(); c++) {
        auto color = colorArray[c];
        for (auto& processedColor : colors) {
          if (color == processedColor) {
            timeslotValidation.addEdge(v, c);
            break;
          }
        }

      }
    }


    long timeslotValidationOutput = timeslotValidation.solve();
    if (timeslotValidationOutput > clique.size()) {
      return false;
    }

    

    for (long v = 0; v < clique.size(); v++) {
    
    
    }

    // since this is a clique, every vertex will have a different room.
    // every vertex is actually a collection of different classrooms.
    // inside each vertex, run bipartite matching to make sure each of them can fit within a timeslot.
    return true;
  }

}// end namespace Coloring::Validator