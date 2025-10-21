#pragma once
#include <vector>

namespace Coloring{
  struct ColoringOutput {
    struct SolutionData {
      long roomColor;
      long timeslotColor;
    };
    enum status {
      NO_SOLUTION, // zero
      HAS_SOLUTION // one
    };
    std::vector<SolutionData> solutionData;
    status solutionStatus;
  };

  constexpr ColoringOutput DEFAULT_NO_SOLUTION = {
      {},
      ColoringOutput::NO_SOLUTION
};
}