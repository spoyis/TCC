#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include "Graph.h"
#include "Coloring.h"
#include <memory>

namespace input {

  class Parser {
  public:
    Parser(const std::string& input) : fileStream(input.c_str()), inputFile(input), keywordString(11) {
      if (!fileStream.is_open()) {
        throw std::ios_base::failure("Failed to open the file: " + input);
      }

      initializeKeywordStrings();
    }

    #define PARSING_ERROR std::runtime_error("Parsing error: Unable to parse data.")
    #define MISSING_DATA std::runtime_error("Parsing error: There is missing data.")
    #define UNEXPECTED_DATA std::runtime_error("Parsing error: Theres more data than expected given the header.")
    void parse() noexcept(false) {
      if (!parseHeader()) throw PARSING_ERROR;
      
      //Graph<edge, vertex> outputGraph(_numVertices);
      std::string line;
      while (std::getline(fileStream, line)) {

        std::istringstream stream(line);
        std::string operation;
        if (line.size() == 0) {
          lineCount++;
          continue;
        }

        if (!(stream >> operation)) {
          std::cerr << "PARSING ERROR AT LINE " << lineCount << '\n';
          throw PARSING_ERROR;
        }

        if (operation == keywordString[label_kw]) {
          if (!parseLabel(stream))
            throw PARSING_ERROR;
        }
        else if (operation == keywordString[edge_kw]) {
          if (!parseEdge(stream))
            throw PARSING_ERROR;
        }
        else if (operation == keywordString[vertex_kw]) {
          if (!parseVertex(stream))
            throw PARSING_ERROR;
        }
        else if (operation.length() >= 2 && operation[0] == '/' && operation[1] == '/') { lineCount++; continue; }
        else {
          std::cerr << "UNEXPECTED TOKEN AT LINE  " << lineCount << " \"" << line << "\"\n";
          throw PARSING_ERROR;
        }

        lineCount++;
      }

      if (_parsedVertices != _numVertices) {
        std::cerr << "EXPECTED [" << _numVertices << "] VERTICES(CLASSES), BUT WAS GIVEN [" << _parsedVertices << "]\n";
        throw _parsedVertices < _numVertices ? MISSING_DATA : UNEXPECTED_DATA;
      }

      if (_parsedEdges != _numEdges) {
        std::cerr << "EXPECTED [" << _numEdges << "] EDGES(RESTRICTIONS), BUT WAS GIVEN [" << _parsedEdges << "]\n";
        throw _parsedEdges < _numEdges ? MISSING_DATA : UNEXPECTED_DATA;
      }

      if (_parsedClassrooms != _numClassrooms) {
        std::cerr << "EXPECTED [" << _numClassrooms << "] CLASSROOMS, BUT WAS GIVEN [" << _parsedClassrooms << "]\n";
        throw  _parsedClassrooms < _numClassrooms ? MISSING_DATA : UNEXPECTED_DATA;
      }

      if (_parsedTimeslots != _numTimeslots) {
        std::cerr << "EXPECTED [" << _numTimeslots << "] TIMESLOTS, BUT WAS GIVEN [" << _parsedTimeslots << "]\n";
        throw  _parsedTimeslots < _numTimeslots ? MISSING_DATA : UNEXPECTED_DATA;
      }

      insertColoringData();
    }

    Coloring::Checker<Graph<long, std::vector<long>>> getChecker() {
      return Coloring::Checker<Graph<long, std::vector<long>>>(_graph, *_roomData);
    }

  private:
    bool parseHeader() {
      std::string line;
      std::getline(fileStream, line);

      std::istringstream stream(line);
      std::string expectedWord;
      int parsedInteger;


      std::vector<std::string> expectedKeyword = {
        keywordString[header_kw],
        keywordString[room_count_kw],
        keywordString[timeslot_count_kw],
        keywordString[vertex_count_kw],
        keywordString[edge_count_kw]
      };

      std::string usage = "LINE 1 should be formatted as: HEADER classrooms <integer> timeslots <integer> classes <integer> restrictions <integer>";
      if(!(stream >> expectedWord) || (expectedWord != keywordString[header_kw])) {
        std::cerr << "Error at line 1, expected " << keywordString[header_kw] << " token\n" << usage;
        return false;
      }

      for (long i = 1; i <= 4; i++) {
        if (!(stream >> expectedWord >> parsedInteger) ||
            (expectedWord != expectedKeyword[i])
        ) {
          std::cerr << "Error at line 1, expected " + expectedKeyword[i] << " <integer>\n" << usage;
          return false;
        }
        
        switch (i) {
        case 1:
          _numClassrooms = parsedInteger;
          break;
        case 2:
          _numTimeslots = parsedInteger;
          break;
        case 3:
          _numVertices = parsedInteger;
          break;
        case 4:
          _numEdges = parsedInteger;
          break;
        }
      }

      _roomData = new std::vector<std::vector<long>>(_numVertices);
      _timeslotData = new std::vector<std::vector<long>>(_numVertices);
      _graph = new Graph<long, std::vector<long>>(_numVertices);
      Graph<long, std::vector<long>>& graphReference = *_graph;

      for (long i = 0; i < _numVertices; i++)
        for (long j = 0; j < _numVertices; j++)
          graphReference[i][j] = 0;

      return true;
    }

    void insertColoringData() {
      Graph<long, std::vector<long>>& graphReference = *_graph;

      for (long i = 0; i < _numVertices; i++) {
        graphReference.setVertexData(i, _timeslotData->operator[](i));
      }
    }

    bool parseLabel(std::istringstream& stream) {
      std::string keyword;
      std::string label;
      std::string expected = keywordString[label_kw] + " labelTarget label\n";
      Graph<long, std::vector<long>>& graphReference = *_graph;

      if (!(stream >> keyword >> label)) {
        std::cerr << "Error while parsing label at line " << lineCount << ", expected: " << expected;
        return false;
      }

      if (keyword == keywordString[room_kw]) {
        classroomLabelMap[label] = classroomLabelMap.size();
        _parsedClassrooms++;
      }
      else if (keyword == keywordString[vertex_kw]) {
        vertexLabelMap[label] = vertexLabelMap.size();
        _parsedVertices++;
        graphReference.setVertexLabel(label, vertexLabelMap.size() - 1);
      }
      else if (keyword == keywordString[timeslot_kw]) {
        timeslotLabelMap[label] = timeslotLabelMap.size();
        _parsedTimeslots++;
      }
      else {
        std::string validLabelTargets("[" + keywordString[room_kw] + " OR " + keywordString[vertex_kw] + " OR " +  keywordString[timeslot_kw] + "]");
        std::cerr << "Error while parsing label at line " << lineCount << ", invalid labelTarget.\n valid label targets are " << validLabelTargets << std::endl;
        return false;
      }
      return true;
    }

    bool parseEdge(std::istringstream& stream) {
      std::string add = "add";
      std::string inputLabel;
      std::string inputAction;
      std::string inputVertex1;
      std::string inputVertex2;
      std::string genericError = "Error while parsing edge data at line " + std::to_string(lineCount);

      if (!(stream >> inputLabel >> inputAction >> inputVertex1 >> inputVertex2)) {
        std::string expected = keywordString[edge_kw] + " [restrictionLabel] add [classLabel1] [classLabel2]\n";
        std::cerr << genericError << ", expected: " << expected;
        return false;
      }

      if (edgeLabelMap.find(inputLabel) != edgeLabelMap.end()) {
        std::cerr << genericError << ", given input label <" << inputLabel << "> was already previously defined.\n";
        return false;
      }

      edgeLabelMap[inputLabel] = edgeLabelMap.size();
      _parsedEdges++;

      if (inputAction != "add") {
        std::cerr << genericError << ", unknown action, expected \"add\"\n";
        return false;
      }

      auto vertex1Id = vertexLabelMap.find(inputVertex1);
      auto vertex2Id = vertexLabelMap.find(inputVertex2);

      if (vertex1Id == vertexLabelMap.end()) {
        std::cerr << genericError << ", given vertex label <" << inputVertex1 << "> not found on vertex label list\n";
        return false;
      }

      if (vertex2Id == vertexLabelMap.end()) {
        std::cerr << genericError << ", given vertex label <" << inputVertex2 << "> not found on vertex label list\n";
        return false;
      }

      Graph<long, std::vector<long>>& graph = *_graph;
      graph[vertex1Id->second][vertex2Id->second] = 1;
      graph[vertex2Id->second][vertex1Id->second] = 1;

      return true;
    }

    bool parseVertex(std::istringstream& stream) {
      std::string add = "add";
      std::string inputLabel;
      std::string inputAction;
      std::string inputTarget;
      std::string inputValue;
      std::string genericError = "Error while parsing vertex data at line " + std::to_string(lineCount);

      if (!(stream >> inputLabel >> inputAction >> inputTarget >> inputValue)) {
        std::string expected = keywordString[vertex_kw] + " [classLabel] add [" + keywordString[room_kw] + " OR " + keywordString[timeslot_kw] + " OR " + keywordString[clique_kw] + "] [targetLabelValue]\n";
        std::cerr << genericError << ", expected: " << expected;
        return false;
      }

      auto vertexId = vertexLabelMap.find(inputLabel);
      if (vertexId == vertexLabelMap.end()) {
        std::cerr << genericError << ", given input label <" << inputLabel << "> not found on vertex label list.\n";
        return false;
      }

      if (inputAction != "add") {
        std::cerr << genericError << ", unknown action, expected \"add\"\n";
        return false;
      }

      if (inputTarget == keywordString[room_kw]) {
        auto labelId = classroomLabelMap.find(inputValue);

        if (labelId == classroomLabelMap.end()) {
          std::cerr << genericError << ", given classroom label <" << inputValue << "> not found on classroom label list\n";
          return false;
        }

        _roomData->operator[](vertexId->second).push_back(labelId->second);
      }
      else if (inputTarget == keywordString[timeslot_kw]) {
        auto labelId = timeslotLabelMap.find(inputValue);

        if (labelId == timeslotLabelMap.end()) {
          std::cerr << genericError << ", given timeslot label <" << inputValue << "> not found on timeslot label list\n";
          return false;
        }

        _timeslotData->operator[](vertexId->second).push_back(labelId->second);
      }
      else if (inputTarget == keywordString[clique_kw]) {
        std::istringstream cliqueStream(inputValue);
        std::string cliqueVertex;
        std::vector<long> cliqueVertexIds = { vertexId->second };
        while (std::getline(cliqueStream, cliqueVertex, '+')) {
          auto cliqueVertexId = vertexLabelMap.find(cliqueVertex);

          if (cliqueVertexId == vertexLabelMap.end()) {
            std::cerr << genericError << ", given clique vertex label <" << cliqueVertex << "> not found on vertex label list.\n";
            return false;
          }

          cliqueVertexIds.push_back(cliqueVertexId->second);
        }

        Graph<long, std::vector<long>>& graph = *_graph;
        
        for (long i = 0; i < cliqueVertexIds.size(); i++) {
          for (long j = i + 1; j < cliqueVertexIds.size(); j++) {
            auto v1 = cliqueVertexIds[i];
            auto v2 = cliqueVertexIds[j];

            graph[v1][v2] = 1;
            graph[v2][v1] = 1;
          }
        }
      }
      else {
        std::cerr << genericError << ", unknown target value. accepted values are [" << keywordString[room_kw] << " OR " << keywordString[timeslot_kw] << " OR "  << keywordString[clique_kw] << "]\n";
        return false;
      }

      return true;
    }

  private:
    int _numEdges;
    int _numVertices;
    int _numClassrooms;
    int _numTimeslots;
    long lineCount{2};
    std::ifstream fileStream;
    std::string inputFile;
    std::map<std::string, int> edgeLabelMap;
    std::map<std::string, int> vertexLabelMap;
    std::map<std::string, int> classroomLabelMap;
    std::map<std::string, int> timeslotLabelMap;
    std::vector<std::vector<long>>* _roomData;
    std::vector<std::vector<long>>* _timeslotData;
    Graph<long, std::vector<long>>* _graph;
    long _parsedEdges{ 0 };
    long _parsedVertices{ 0 };
    long _parsedClassrooms{ 0 };
    long _parsedTimeslots{ 0 };

    enum keywords {
      header_kw,
      room_count_kw,
      timeslot_count_kw,
      vertex_count_kw,
      edge_count_kw,
      edge_kw,
      vertex_kw,
      room_kw,
      timeslot_kw,
      label_kw,
      clique_kw
    };
    std::vector<std::string> keywordString;

    void initializeKeywordStrings() {
      keywordString[header_kw] = "HEADER";
      keywordString[room_count_kw] = "classrooms";
      keywordString[timeslot_count_kw] = "timeslots";
      keywordString[edge_count_kw] = "restrictions";
      keywordString[vertex_count_kw] = "classes";

      keywordString[edge_kw] = "restriction";
      keywordString[vertex_kw] = "class";
      keywordString[room_kw] = "room";
      keywordString[timeslot_kw] = "timeslot";
      keywordString[label_kw] = "label";
      keywordString[clique_kw] = "clique";
    }
  };
}