/*
 * Utility functions
 *
 * (c) 2015 Schrodinger, Inc.
 */

#include <string>
#include <vector>

#include <string.h>

std::vector<std::string> strsplit(const std::string &s, char delim=0);

bool cstrlessnat(const char * a, const char * b);
bool strlessnat(const std::string& a, const std::string& b);

/*
 * C string comparison class
 */
struct cstrless_t {
  bool operator()(const char * a, const char * b) const {
    return strcmp(a, b) < 0;
  }
};
