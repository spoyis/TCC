#pragma once
#include "Coloring.h"
#include "Graph.h"
#include "BipartiteMatching.h"
#include <map>
#include <set>
namespace Coloring::Validator { // begin namespace Coloring::Validator

  namespace color {
    const std::string red("\033[0;31m");
    const std::string green("\033[1;32m");
    const std::string yellow("\033[1;33m");
    const std::string blue("\033[1;34m");
    const std::string cyan("\033[0;36m");
    const std::string magenta("\033[0;35m");
    const std::string reset("\033[0m");
  }

  static const std::string PASS = "- " + color::green + "[PASSED]: " + color::reset;
  static const std::string FAIL = "- " + color::red + "[FAILED]: " + color::reset;

  enum ValidationState {
    VALID,
    INVALID_TIMESLOT_PIGEONHOLE,
    INVALID_TIMESLOT_COLORING,
    INVALID_ROOM_PIGEONHOLE,
    INVALID_ROOM_COLORING,
    INVALID_JOINING,
    MISSING_COLORS,
  };
  
  using clique_t = std::vector<long>;

  // validates a clique
  template<typename edge, typename vertex>
  ValidationState clique(
    clique_t clique,
    Graph<edge, vertex>& g,
    std::vector<std::vector<long>> roomData,
    bool shouldSaveColoringData = false,                    // optional input
    std::vector<std::pair<long, long>>* timeslotMatches = nullptr,  // optional output
    std::map<long, std::vector<std::pair<long, long>>>* roomMatches = nullptr // optional output
  ) {
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
        //if (timeslotMatches)
         // std::cout << "CLIQUE ROOT " << clique[v] << " ACCEPTS " << color << '\n';
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

    // optionally store timeslot matching
    if (shouldSaveColoringData && timeslotMatches) {
      *timeslotMatches = timeslotValidation.getMatching(clique);
      for (long i = 0; i < timeslotMatches->size(); i++) {
        (*timeslotMatches)[i].second = colors[(*timeslotMatches)[i].second];
      }
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

      // optionally store room matching
      if (shouldSaveColoringData && roomMatches) {
        (*roomMatches)[cliqueGroup.first] = roomValidation.getMatching(cliqueGroupArray);
        // loop over the vector of pairs and remap the second element to the correct value
        for (long i = 0; i < (*roomMatches)[cliqueGroup.first].size(); i++) {
          long localIndex = (*roomMatches)[cliqueGroup.first][i].second;
          (*roomMatches)[cliqueGroup.first][i].second = joinedRoomArray[localIndex];
        }
      }
    }

    return VALID;
  }

  template<typename edge, typename vertex>
  bool coloringOutput(
    const std::vector<Coloring::ColoringOutput::SolutionData>& output,
    Graph<edge, vertex>& g,
    const std::vector<std::vector<long>>& roomData
  ) {
    std::cout << "[VALIDATION] BEGIN VALIDATION ON CURRENT SOLUTION\n";
    long vertexCount = g.getVertexCount();

    // --- [1] Every vertex should have a valid timeslot color ---
    for (long v = 0; v < vertexCount; ++v) {
      long root = g.getRoot(v);
      long timeslotColor = output[v].timeslotColor;

      const auto& allowedTimeslots = g.getVertexData(root);

      bool found = false;
      for (auto t : allowedTimeslots) {
        //std::cout << "ALLOWED TIMESLOTS FOR ROOT " << g.getRoot(v) << ": " << t << '\n';
        if (t == timeslotColor) {
          found = true;
          break;
        }
      }

      if (!found) {
        std::cerr << FAIL <<"[ERROR] Vertex " << v << " has invalid timeslot color " << timeslotColor
          << " not in root " << root << " allowed set.\n";
        return false;
      }
    }

    std::cout << PASS << "TIMESLOT VALIDATION " << std::endl;

    // --- [2] Every vertex should have a valid room color ---
    for (long v = 0; v < vertexCount; ++v) {
      long roomColor = output[v].roomColor;

      const auto& allowedRooms = roomData[v];
      bool found = false;
      for (auto r : allowedRooms) {
        if (r == roomColor) {
          found = true;
          break;
        }
      }

      if (!found) {
        std::cerr << FAIL << "[ERROR] Vertex " << v << " has invalid room color " << roomColor << ".\n";

        for (long i = 0; i < roomData[v].size(); i++) {
          std::cout << "ALLOWED COLORS ARE: " << roomData[v][i] << '\n';
        }
        return false;
      }
    }

    std::cout << PASS << "ROOM VALIDATION " << std::endl;

    // --- [3] No two (timeslot, room) pairs should repeat ---
    std::set<std::pair<long, long>> usedPairs;

    for (long v = 0; v < vertexCount; ++v) {
      long t = output[v].timeslotColor;
      long r = output[v].roomColor;
      auto pair = std::make_pair(t, r);

      if (usedPairs.find(pair) != usedPairs.end()) {
        std::cerr << FAIL << "[ERROR] Duplicate (timeslot, room) pair found: ("
          << t << ", " << r << ") for vertex " << v << ".\n";
        return false;
      }
      else {
        usedPairs.insert(pair);
      }
    }
    std::cout << PASS << "UNIQUE COLOR PAIR VALIDATION " << std::endl;
    // If we reach here, all checks passed
    return true;
  }

  template <typename edge, typename vertex>
  bool canJoinBasedOnRoomConstraint(
    long v1, long v2, 
    Graph<edge, vertex>& g, 
    const std::vector<std::vector<long>>& roomData
  ) {
    long r1 = g.getRoot(v1);
    long r2 = g.getRoot(v2);

    if (r1 == r2) return true; // already joined...
    if (g[r1][r2]) return false; // already adjacent...

    std::vector<long> vertexList;
    long vertexCount = g.getVertexCount();

    for (long i = 0; i < vertexCount; i++) {
      auto vertexRoot = g.getRoot(i);
      if (vertexRoot == r1 || vertexRoot == r2)
        vertexList.push_back(i);
    }

    std::set<long> joinedRoomsSet;
    for (auto v : vertexList) {
      for (auto room : roomData[v]) {
        joinedRoomsSet.insert(room);
      }
    }

    std::vector<long> joinedRooms(joinedRoomsSet.begin(), joinedRoomsSet.end());
    // pigeonhole violation — not enough unique rooms for distinct classes
    if (joinedRooms.size() < vertexList.size()) {
      return false;
    }

    // bipartite matching similar to joined clique vertex checking above
    BipartiteMatching roomValidation(vertexList.size(), joinedRooms.size());

    for (long v = 0; v < vertexList.size(); v++) {
      // we know what colors(rooms) 'v' accepts, but we need to map its indexes to be the same as roomData indexed them.
      auto& colorArray = roomData[vertexList[v]];
      roomValidation.defineInitialCapacity(v, 1);

      for (auto color : colorArray) {
        for (long c = 0; c < joinedRooms.size(); c++) {
          if (color == joinedRooms[c]) {
            roomValidation.addEdge(v, c);
            break;
          }
        }
      }
    }

    long roomValidationOutput = roomValidation.solve();
    return roomValidationOutput == vertexList.size();
  }

  // quickly checks if any vertices are missing colors, immediately invalidating the whole instance
  // this is used only on the root vertex of a Checker instance.
  template<typename edge, typename vertex>
  ValidationState naiveCheck(Graph<edge, vertex>& g, std::vector<std::vector<long>> roomData) {
    ValidationState output = VALID;
    auto vertexCount = g.getVertexCount();

    for (long i = 0; i < g.getVertexCount(); i++) {
      if (g.getVertexData(i).size() == 0 || roomData[i].size() == 0) {
        std::cout << "[VALIDATION][BAD GRAPH] VERTEX " << g.getVertexLabel(i) << " IS MISSING COLORS!\n";
        output = MISSING_COLORS;
      }
    }

    return output;
  }

}// end namespace Coloring::Validator