
template <typename edge, typename vertex>
class Graph {
public:
  class AdjacencyMatrix {
  public:
    edge* data;       // flat array of size vertexCount * vertexCount
    int* degrees;     // pointer to Graph's degree array
    int _vertexCount;

    // ----- Constructor -----
    AdjacencyMatrix(int vertexCount, int* degArray)
      : _vertexCount(vertexCount), degrees(degArray)
    {
      data = new edge[vertexCount * vertexCount](); 
    }

    // ----- Copy / clone -----
    AdjacencyMatrix* clone() {
      AdjacencyMatrix* copy = new AdjacencyMatrix(_vertexCount, degrees);
      std::copy(data, data + (_vertexCount * _vertexCount), copy->data);
      return copy;
    }

    // ----- Row proxy -----
    struct RowProxy {
      AdjacencyMatrix* mat;
      int row;

      struct EdgeProxy {
        RowProxy* rowProxy;
        int col;

        EdgeProxy& operator=(edge value) {
          edge& old = *(rowProxy->mat->data + rowProxy->row * rowProxy->mat->_vertexCount + col);
          if (old == 0 && value != 0) {        // edge added
            rowProxy->mat->degrees[rowProxy->row]++;
            rowProxy->mat->degrees[col]++;
          }
          else if (old != 0 && value == 0) { // edge removed
            rowProxy->mat->degrees[rowProxy->row]--;
            rowProxy->mat->degrees[col]--;
          }
          old = value;
          return *this;
        }

        operator edge() const {
          return *(rowProxy->mat->data + rowProxy->row * rowProxy->mat->_vertexCount + col);
        }
      };

      EdgeProxy operator[](int col) {
        return EdgeProxy{ this, col };
      }
    };

    RowProxy operator[](int row) {
      return RowProxy{ this, row };
    }

    // ----- Destructor -----
    ~AdjacencyMatrix() {
      delete[] data;
    }
  };
};