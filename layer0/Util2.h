/*
 * Utility functions
 *
 * (c) 2015 Schrodinger, Inc.
 */

#include <string>
#include <vector>

std::vector<std::string> strsplit(const std::string &s, char delim=0);

bool cstrlessnat(const char * a, const char * b);
bool strlessnat(const std::string& a, const std::string& b);
