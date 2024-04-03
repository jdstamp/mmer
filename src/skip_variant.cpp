//
// Created by Julian Stamp on 4/3/24.
//

#include "skip_variant.h"
bool skip_variant(const std::vector<int> &ind, int i) {
  if (ind.empty()) {
    // If the vector is empty, there are no elements to skip
    return false;
  }

  // Look for i+1 because R uses 1-based indexing
  // If there is no match, std::find returns ind.end()
  return std::find(ind.begin(), ind.end(), i + 1) == ind.end();
}