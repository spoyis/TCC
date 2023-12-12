#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Graph.h"
// http://archive.dimacs.rutgers.edu/pub/challenge/graph/doc/ccformat.tex

namespace DIMACS {
  template<class graph_t>
  class Parser {
  private:
    enum operations {
      problem_op = 'p', node_op = 'n', edge_op = 'e', comment_op = 'c', solution_op = 's', bound_op = 'b'
    };

  public:
    using edge = typename graph_t::edge_t;
    using vertex = typename graph_t::vertex_t;

    Parser() :fileStream{std::cin}, inputFile("") {}

    Parser(const std::string& input) : fileStream(std::ifstream(input)), inputFile(input){
      if (!fileStream.is_open()) {
        throw std::ios_base::failure("Failed to open the file: " + input);
      }
    }

    #define PARSING_ERROR std::runtime_error("Parsing error: Unable to parse data.")
    Graph<edge,vertex> parse() noexcept(false) {
      if (!parseProblemHeader()) throw PARSING_ERROR;

      Graph<edge, vertex> parsedGraph(_numVertices);

      int lineIndex = 2;
      int remainingNodes = _numVertices;
      int remainingEdges = _numEdges;
      int lastCode;

      while (remainingNodes) {
        lastCode = parseNode(lineIndex++, parsedGraph);
        if (lastCode == 0) throw PARSING_ERROR;
        if (lastCode != comment_op) remainingNodes--;
      }

      while (remainingEdges) {
        lastCode = parseEdge(lineIndex++, parsedGraph);
        if (lastCode == 0) throw PARSING_ERROR;
        if (lastCode != comment_op) remainingEdges--;
      }

      return parsedGraph;
    }

  private:
    int _numEdges;
    int _numVertices;
    std::istream& fileStream;
    std::string inputFile;

    int parseProblemHeader() {
      std::string line;
      std::getline(fileStream, line);

      std::istringstream stream(line);

      std::string prefixOp, prefixWord;
      std::string expectedPrefixOp(1, static_cast<char>(problem_op));
      std::string expectedPrefixWord("edge");
      
      int edgesValue, nodesValue;

      if (!(stream >> prefixOp >> prefixWord >> nodesValue >> edgesValue) || 
            prefixOp != expectedPrefixOp ||
            prefixWord != expectedPrefixWord
        ) {
        std::cerr << "Error at line 1, bad formatting\n use: p edge NODES EDGES" << std::endl;
        return 0;
      }

      if (edgesValue < 0 || nodesValue < 0) {
        std::cerr << "Error at line 1, can't have negative values" << std::endl;
        return 0;
      }

      _numEdges = edgesValue;
      _numVertices = nodesValue;

      return 1;
    }

    int parseNode(int lineIndex, Graph<edge,vertex>& g) {
      std::string line;
      std::getline(fileStream, line);

      std::istringstream stream(line);

      std::string prefixOp;
      std::string expectedPrefixOp(1, static_cast<char>(node_op));
      std::string commentOp(1, static_cast<char>(comment_op));

      int indexValue;
      vertex nodesValue;

      stream >> prefixOp;
      if (prefixOp == commentOp) return comment_op;

      if (!(stream >> indexValue >> nodesValue) ||
          prefixOp != expectedPrefixOp && prefixOp != commentOp 
        ) {
        std::cerr << "Error at line " << lineIndex << ", bad formatting\n use : n ID VALUE" << std::endl;
        return 0;
      }

      if (indexValue < 0 || indexValue >= _numVertices) {
        std::cerr << "Error at line " << lineIndex << ", the given index is out of bounds" << std::endl;
        return 0;
      }

      g.setVertexData(indexValue, nodesValue);      

      return 1;
    }

    int parseEdge(int lineIndex, Graph<edge, vertex>& g) {
      std::string line;
      std::getline(fileStream, line);

      std::istringstream stream(line);

      std::string prefixOp;
      std::string expectedPrefixOp(1, static_cast<char>(edge_op));
      std::string commentOp(1, static_cast<char>(comment_op));

      int indexValue1;
      int indexValue2;
      edge edgeValue;

      stream >> prefixOp;
      if (prefixOp == commentOp) return comment_op;

      if (!(stream >> indexValue1 >> indexValue2 >> edgeValue) ||
        prefixOp != expectedPrefixOp && prefixOp != commentOp
        ) {
        std::cerr << "Error at line " << lineIndex << ", bad formatting\n use : e INDEX1 INDEX2 VALUE" << std::endl;
        return 0;
      }

      if (indexValue1 < 0 || indexValue1 >= _numVertices || indexValue2 < 0 || indexValue2 >= _numVertices) {
        std::cerr << "Error at line " << lineIndex << ", one or both indexes out of bounds" << std::endl;
        return 0;
      }

      g[indexValue1][indexValue2] = edgeValue;
      g[indexValue2][indexValue1] = edgeValue;

      return 1;
    }
  };
}