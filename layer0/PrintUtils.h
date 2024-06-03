#pragma once
#include <type_traits>
#include <vector>
#include <string>
#include <sstream>

/***********************************************************************
 * Table pretty printing utility class
 * ---------------------------------------------------------------------
 * Will autofill to the row with the largest column
 ***********************************************************************/
class display_table_t {
public:
  // both begins and ends a row
  display_table_t & begin_row();

  // inserts a new cell with some printable type in it
  template <typename T>
  display_table_t & insert_cell(T s) {
    std::stringstream ss;
    ss << s;
    _table[_current_row].emplace_back(ss.str());
    return *this;
  }

  // finalizes the table and prints it to std out
  void display();
private:
  int _current_row { -1 };
  std::vector<std::vector<std::string>> _table;
};
