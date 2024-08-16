#include <iostream>
#include <iomanip>
#include <algorithm>
#include "PrintUtils.h"

display_table_t & display_table_t::begin_row() {
  _table.push_back(std::vector<std::string>());
  _current_row++;
  return *this;
}

void display_table_t::display() {
  // Dimensions of tables
  const size_t num_rows = _table.size();
  const size_t num_cols = ([&] () {
    size_t largest = 0;
    for (auto & v : _table) {
      largest = std::max(largest, v.size());
    }
    return largest;
  })();

  // Pad
  for (auto & row : _table) {
    while (row.size() < num_cols) {
      row.emplace_back(" ");
    }
  }

  // Get the largest string for each col
  std::vector<size_t> col_sizes(num_cols);
  for (size_t j = 0; j < num_cols; ++j) {
    size_t largest = 0;
    for (size_t i = 0; i < num_rows; ++i) {
      largest = std::max(largest, _table[i][j].size());
    }
    col_sizes[j] = largest;
  }


  // construct the output
  std::stringstream ss;
  ss << std::left;
  auto insert_hr = [&]() {
    for (size_t j = 0; j < num_cols; ++j) {
      ss << std::setw(col_sizes[j] + 3)
      << std::setfill('-') << "+";
    }
    ss << "+" << std::setfill(' ') << std::endl;
  };
  {
    int i = 0, j = 0;
    for (auto & row : _table) {
      insert_hr();
      ss << "| ";
      for (auto & col : row) {
        ss << std::setw(col_sizes[j])
           << col << " | ";
        j++;
      }
      ss << std::endl;
      i++; j = 0;
    }
    insert_hr();
  }
  std::cout << ss.str();
}
