#pragma once
#include <chrono>
#include <fstream>
#include <string>

struct OutputWriter {

  struct GraphData {
    long vertex_count;
    long edge_count;
  };

  GraphData originalGraph{ 0,0 };
  GraphData preProcessedGraph{ 0,0 };
  long colors;
  std::chrono::duration<double> executionTime;
  long exploredVertices;
  long prunedTimeslot;
  long prunedRooms;

  void writeToFile(const std::string& outputName) const {
    std::ofstream out(outputName);
    if (!out.is_open()) {
      throw std::runtime_error("Failed to open output file: " + outputName);
    }

    out << originalGraph.vertex_count << "\n"
      << originalGraph.edge_count << "\n"
      << preProcessedGraph.vertex_count << "\n"
      << preProcessedGraph.edge_count << "\n"
      << colors << "\n"
      << executionTime.count() << "\n"
      << exploredVertices << "\n"
      << prunedTimeslot << "\n"
      << prunedRooms << "\n";

    out.close();
  }
};