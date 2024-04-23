#pragma once
#include <iostream>

template<typename T>
class UnionFind
{ // begin class UnionFind
public:
  UnionFind(long n) : _size(n){
		parent = new long[_size];
    setSize = new long[_size];
		data = new T[_size];

    for (long i = 0; i < _size; i++) { parent[i] = i; setSize[i] = 1; }
	}

	// Move constructor
	UnionFind(UnionFind&& other) noexcept {
		this->parent = other.parent;
		this->data = other.data;
		this->_size = other._size;
    this->setSize = other.setSize;

    other.parent = nullptr;
    other.data = nullptr;
    other.setSize = nullptr;
	}

  ~UnionFind() {
    delete[] data;
    delete[] parent;
    delete[] setSize;
  }

  // Overloaded for lvalue
  T& operator[](std::size_t index) {
    return findOp(index);
  }

  // Overloaded for rvalue
  T operator[](std::size_t index) const {
    return findOp(index);
  }

  // Union operation
  void unionOp(std::size_t index1, std::size_t index2) {
    long root1 = findRoot(index1);
    long root2 = findRoot(index2);

    if (root1 != root2) {
      // Prioritize the smallest index as the root.
      if (root1 < root2) {
        parent[root2] = root1;
        setSize[root2] += setSize[root1];
      }
      else {
        parent[root1] = root2;
        setSize[root1] += setSize[root2];
      }
    }
  }

  long findRoot(long index) {
    while (parent[index] != index) {
      // Path compression
      parent[index] = parent[parent[index]];
      index = parent[index];
    }
    return index;
  }

  long getSize(long index) {
    return setSize[findRoot(index)];
  }

  UnionFind* clone() {
    UnionFind* copy = new UnionFind(_size);

    for (long i = 0; i < _size; i++) {
      copy->data[i] = this->data[i];
      copy->parent[i] = this->parent[i];
      copy->setSize[i] = this->setSize[i];
    }

    return copy;
  }

private:
  long _size;
  long* parent;
  long* setSize;
  T* data;

  T findOp(std::size_t index) const {
    return data[findRoot(index)];
  }

  T& findOp(std::size_t index) {
    return data[findRoot(index)];
  }
 
}; // end class UnionFind