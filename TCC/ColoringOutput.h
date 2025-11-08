#pragma once
#include <array>

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

  inline const ColoringOutput DEFAULT_NO_SOLUTION = {
      {},
      ColoringOutput::NO_SOLUTION
  };
}