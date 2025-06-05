#pragma once
#include "Coloring.h"
#include "Graph.h"
#include "BipartiteMatching.h"
#include <map>

namespace Coloring::Validator { // begin namespace Coloring::Validator

  enum ValidationState {
    VALID,
    INVALID_TIMESLOT_PIGEONHOLE,
    INVALID_TIMESLOT_COLORING,
    INVALID_ROOM_PIGEONHOLE,
    INVALID_ROOM_COLORING
  };
  
  using clique_t = std::vector<long>;

  // validates a clique
  template<typename edge, typename vertex>
  ValidationState clique(clique_t clique, Graph<edge,vertex>& g, std::vector<std::vector<long>> roomData) {
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
      return INVALID_TIMESLOT_PIGEONHOLE;

    }
 
    BipartiteMatching timeslotValidation(clique.size(), colors.size());
    //std::cout << "[DEBUG] VALIDATING TIMESLOTS\n";
    for (long v = 0; v < clique.size(); v++) {
      auto& colorArray = g.getVertexData(clique[v]);
      timeslotValidation.defineInitialCapacity(v, 1);
      for (auto color : colorArray) {
        for (long c = 0; c < colors.size(); c++) {
          if (color == colors[c]) {
            //std::cout << "[DEBUG] ADDING EDGE FROM " << v << " TO " << c << '\n';
            timeslotValidation.addEdge(v, c);
            break;
          }
        }

      }
    }


    long timeslotValidationOutput = timeslotValidation.solve();
    if (timeslotValidationOutput < clique.size()) {
      /*std::cout << "[DEBUG][TIME COLOR][PRUNING]: " << timeslotValidationOutput << " < " << clique.size() << '\n';

      for (auto& processedColor : colors) {
        std::cout << "[DEBUG][TIME COLOR] PROCESSED TIME COLOR: " << processedColor << '\n';
      }


      for (auto v : clique) {
        std::cout << "[DEBUG][TIME COLOR] V: " << v << " LABEL === " << g.getVertexLabel(v) << '\n';
        auto& colorArray = g.getVertexData(v);
        for (auto& c : colorArray) {
          std::cout << "[DEBUG][TIME COLOR] ORIGINAL TIME COLOR: " << c << '\n';
        }
      }
      */
      return INVALID_TIMESLOT_COLORING;
    }


    // BEGIN ROOM VALIDATION
    
    std::map<long, std::vector<long>> cliqueVertexMap;

    for (long v = 0; v < clique.size(); v++) {
      auto vertex = clique[v];
      cliqueVertexMap[vertex] = std::vector<long>();
    }

    for (long i = 0; i < g.getVertexCount(); i++) {
      auto root = g.getRoot(i);

      auto mappedVertex = cliqueVertexMap.find(root);
      if (mappedVertex != cliqueVertexMap.end()) {
        (*mappedVertex).second.push_back(i);
      }
    }

    // for each vertex in the clique
    for (auto& cliqueGroup : cliqueVertexMap) {
      auto& cliqueGroupArray = cliqueGroup.second; 
      auto vertexQuantity = cliqueGroupArray.size();
      std::vector<long> joinedRoomArray;

      // find all original vertexes that were merged into one.
      for (auto v : cliqueGroupArray) {
        auto& roomArray = roomData[v];

        // map out all possible rooms
        for (auto room : roomArray) {
          bool flag = true;
          for (auto filteredRoom : joinedRoomArray) {
            if (filteredRoom == room) {
              flag = false;
              break;
            }
          }

          if (flag) joinedRoomArray.push_back(room);
        }

      
      }

      // pigeonhole principle
      if (joinedRoomArray.size() < cliqueGroupArray.size()) {
          return INVALID_ROOM_PIGEONHOLE;
      }
      
      BipartiteMatching roomValidation(cliqueGroupArray.size(), joinedRoomArray.size());

      for (long v = 0; v < cliqueGroupArray.size(); v++) {
        // we know what colors(rooms) 'v' accepts, but we need to map its indexes to be the same as joinedRoomArray indexed them.
        auto& colorArray = roomData[cliqueGroupArray[v]];
        roomValidation.defineInitialCapacity(v, 1);

        for (auto color : colorArray) {
          for (long c = 0; c < joinedRoomArray.size(); c++){
            if (color == joinedRoomArray[c]) {
              
              roomValidation.addEdge(v, c);
              break;
            }
          }

        }
      }

      long roomValidationOutput = roomValidation.solve();
      if (roomValidationOutput < cliqueGroupArray.size()) {
        /*
        std::cout << "ROOM COLOR PRUNING: " << roomValidationOutput << " < " << cliqueGroupArray.size() << '\n';

        for (auto& processedColor : joinedRoomArray) {
          std::cout << "PROCESSED ROOM COLOR: " << processedColor << '\n';
        }

        for (auto v : cliqueGroupArray) {
          std::cout << "ROOM COLOR V: " << v << " LABEL === " << g.getVertexLabel(v) << '\n';
          auto& colorArray = roomData[v];
          for (auto& c : colorArray) {
            std::cout << "ORIGINAL ROOM COLOR: " << c << '\n';
          }
        }
        */

        return INVALID_ROOM_COLORING;
      }
    }

    // since this is a clique, every vertex will have a different room.
    // every vertex is actually a collection of different classrooms.
    // inside each vertex, run bipartite matching to make sure each of them can fit within a timeslot.
    return VALID;
  }

}// end namespace Coloring::Validator