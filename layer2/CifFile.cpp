/*
 * CIF tokenizer
 *
 * All keys are canonicalized to lowercase
 *
 * (c) 2014 Schrodinger, Inc.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include <variant>

#include "CifFile.h"
#include "File.h"
#include "MemoryDebug.h"
#include "strcasecmp.h"

#if !defined(_PYMOL_NO_MSGPACKC)
#include <msgpack.hpp>
#endif

namespace pymol {
namespace _cif_detail {

template <> const char* raw_to_typed(const char* s) { return s; }
template <> std::string raw_to_typed(const char* s) { return s; }
template <> char        raw_to_typed(const char* s) { return s[0]; }
template <> int         raw_to_typed(const char* s) { return atoi(s); }

/**
 * Convert to floating point number, ignores uncertainty notation
 * 1.23(45)e2 -> 1.23e2
 */
template <> double raw_to_typed(const char* s)
{
  const char *close, *open = strchr(s, '(');
  if (open && (close = strchr(open, ')'))) {
    return atof(std::string(s, open - s).append(close + 1).c_str());
  }
  return atof(s);
}

template <> float raw_to_typed(const char* s)
{
  return static_cast<float>(raw_to_typed<double>(s));
}

} // namespace _cif_detail

// basic IO and string handling

// Return true if "c" is whitespace or null
static bool iswhitespace0(char c) {
  return strchr(" \t\r\n", c) ? true : false;
}

// Return true if "c" is whitespace
static bool iswhitespace(char c) {
  return (c && iswhitespace0(c));
}

// Return true if "c" is line feed or carriage return
static bool islinefeed(char c) {
  return (c == '\r' || c == '\n');
}

// Return true if "c" is line feed or carriage return or null
static bool islinefeed0(char c) {
  return (!c || islinefeed(c));
}

// Return true if "c" is double or single quote
static bool isquote(char c) {
  return (c == '"' || c == '\'');
}

// FreeBSD name conflict
#ifdef isspecial
#undef isspecial
#endif

// Return true if token is a STAR keyword
static bool isspecial(const char *token) {
  return (token[0] == '_'
      || strncasecmp("data_", token, 5) == 0
      || strncasecmp("save_", token, 5) == 0
      || strcasecmp("loop_", token) == 0
      || strcasecmp("stop_", token) == 0
      || strcasecmp("global_", token) == 0);
}

// convert all chars to lowercase
static void tolowerinplace(char *p) {
  for (; *p; p++) {
    if (*p <= 'Z' && *p >= 'A')
      *p -= 'Z' - 'z';
  }
}

// CIF stuff

static const cif_array EMPTY_ARRAY(nullptr);

/*
 * Class to store CIF loops. Only for parsing, do not use in any higher level
 * reading functions.
 */
class cif_loop {
public:
  int ncols;
  int nrows;
  const char **values;

  // methods
  const char * get_value_raw(int row, int col) const;
};

// get table value, return nullptr if indices out of bounds
const char * cif_loop::get_value_raw(int row, int col) const {
  if (row >= nrows)
    return nullptr;
  return values[row * ncols + col];
}

// get the number of elements in this array
unsigned cif_array::size() const {
  if (auto arr = std::get_if<cif_detail::cif_str_array>(&m_array)) {
    return (arr->col == cif_detail::cif_str_array::NOT_IN_LOOP)
               ? 1
               : arr->pointer.loop->nrows;
  } else if (auto arr = std::get_if<cif_detail::bcif_array>(&m_array)) {
    return arr->m_arr.size();
  }
  return 0;
}

/// Get array value, return nullptr if `pos >= size()` or value in ['.', '?']
const char* cif_detail::cif_str_array::get_value_raw(unsigned pos) const
{
  if (col == NOT_IN_LOOP)
    return (pos > 0) ? nullptr : pointer.value;
  return pointer.loop->get_value_raw(pos, col);
}

// true if all values in ['.', '?']
bool cif_array::is_missing_all() const {
  for (unsigned i = 0, n = size(); i != n; ++i) {
    if (!is_missing(i))
      return false;
  }

  return true;
}

/**
 * Get a pointer to array or nullptr if not found
 *
 * Can lookup different aliases, the first one found is returned.
 * Also supports an alias shortcut for the trivial case where mmCIF uses
 * a colon and CIF uses an underscore: `get_arr("_foo?bar")` is identical to
 * `get_arr("_foo.bar", "_foo_bar")`
 *
 * @param key data name, must be lower case
 */
const cif_array * cif_data::get_arr(const char * key) const {
  if (auto data = std::get_if<pymol::cif_detail::cif_str_data>(&m_data)) {
    const auto& dict = data->m_dict;
    const char* p = strchr(key, '?');
    std::remove_reference_t<decltype(dict)>::const_iterator it;

#ifndef NDEBUG
    for (const char* q = key; *q; ++q) {
      assert("key must be lower case" && !('Z' >= *q && *q >= 'A'));
    }
#endif

    // support alias shortcut: '?' matches '.' and '_'
    if (p != nullptr) {
      std::string tmp(key);
      // replace '?' by '.' or '_'
      tmp[p - key] = '.';
      if ((it = dict.find(tmp.c_str())) != dict.end())
        return &it->second;
      tmp[p - key] = '_';
      if ((it = dict.find(tmp.c_str())) != dict.end())
        return &it->second;
    } else {
      if ((it = dict.find(key)) != dict.end())
        return &it->second;
    }
  } else if (auto data = std::get_if<pymol::cif_detail::bcif_data>(&m_data)) {

    const auto& dict = data->m_dict;

    std::string_view keyView(key);
    auto split_key = [](const char c) {
      return c == '.' /*|| c == '_'*/ || c == '?';
    };
    auto splitTokenIt = std::find_if(keyView.begin(), keyView.end(), split_key);
    if (splitTokenIt == keyView.end()) {
      return nullptr;
    }
    auto dist = std::distance(keyView.begin(), splitTokenIt);
    auto categoryView = keyView.substr(0, dist);
    auto categoryStr = std::string(categoryView);
    auto categoryIt = dict.find(categoryStr.c_str());
    if (categoryIt == dict.end()) {
      return nullptr;
    }
    auto& category = categoryIt->second;
    auto columnView = keyView.substr(dist + 1);
    auto columnStr = std::string(columnView);
    auto columnIt = category.find(columnStr.c_str());
    if (columnIt == category.end()) {
      return nullptr;
    }
    return &columnIt->second;
  }

  return nullptr;
}

const char* cif_data::code() const
{
  if (auto data = std::get_if<pymol::cif_detail::cif_str_data>(&m_data)) {
    return data->m_code ? data->m_code : "";
  }
  return "";
}

const cif_array* cif_data::empty_array() {
  return &EMPTY_ARRAY;
}

const cif_detail::cif_str_data* cif_data::get_saveframe(const char* code) const {
  if (auto data = std::get_if<pymol::cif_detail::cif_str_data>(&m_data)) {
    const auto& saveframes = data->m_saveframes;
    auto it = saveframes.find(code);
    if (it != saveframes.end())
      return &it->second;
  }
  return nullptr;
}

bool cif_file::parse_file(const char* filename) {
  char* contents = FileGetContents(filename, nullptr);

  if (!contents) {
    error(std::string("failed to read file ").append(filename).c_str());
    return false;
  }

  return parse(std::move(contents));
}

bool cif_file::parse_string(const char* contents) {
  return parse(std::move(mstrdup(contents)));
}

void cif_file::error(const char* msg) {
  std::cout << "ERROR " << msg << std::endl;
}

// constructor
cif_file::cif_file(const char* filename, const char* contents_) {
  if (contents_) {
    parse_string(contents_);
  } else if (filename) {
    parse_file(filename);
  }
}

// constructor
cif_file::cif_file() = default;
cif_file::cif_file(cif_file&&) = default;

// move assignment
cif_file& cif_file::operator=(cif_file&&) = default;

// destructor
cif_file::~cif_file() = default;

bool cif_file::parse(char*&& p) {
  m_datablocks.clear();
  m_tokens.clear();
  m_contents.reset(p);

  if (!p) {
    error("parse(nullptr)");
    return false;
  }

  auto& tokens = m_tokens;
  char quote;
  char prev = '\0';

  std::vector<bool> keypossible;

  // tokenize
  while (true) {
    while (iswhitespace(*p))
      prev = *(p++);

    if (!*p)
      break;

    if (*p == '#') {
      while (!(islinefeed0(*++p)));
      prev = *p;
    } else if (isquote(*p)) { // will nullptr the closing quote
      quote = *p;
      keypossible.push_back(false);
      tokens.push_back(p + 1);
      while (*++p && !(*p == quote && iswhitespace0(p[1])));
      if (*p)
        *(p++) = 0;
      prev = *p;
    } else if (*p == ';' && islinefeed(prev)) {
      // multi-line tokens start with ";" and end with "\n;"
      // multi-line tokens cannot be keys, only values.
      keypossible.push_back(false);
      tokens.push_back(p + 1);
      // advance until `\n;`
      while (*++p && !(islinefeed(*p) && p[1] == ';'));
      // step to next line and null the line feed
      if (*p) {
        *p = 0;
        // \r\n on Windows)
        if (p - 1 > tokens.back() && *(p - 1) == '\r') {
          *(p - 1) = 0;
        }
        p += 2;
      }
      prev = ';';
    } else { // will null the whitespace
      char * q = p++;
      while (!iswhitespace0(*p)) ++p;
      prev = *p;
      if (p - q == 1 && (*q == '?' || *q == '.')) {
        // store values '.' (inapplicable) and '?' (unknown) as null-pointers
        q = nullptr;
        keypossible.push_back(false);
      } else {
        if (*p)
          *(p++) = 0;
        keypossible.push_back(true);
      }
      tokens.push_back(q);
    }
  }

  cif_detail::cif_str_data* current_frame = nullptr;
  std::vector<cif_detail::cif_str_data*> frame_stack;
  std::unique_ptr<cif_data> global_block;
  decltype(m_datablocks) datablocksnew;

  // parse into dictionary
  for (unsigned int i = 0, n = tokens.size(); i < n; i++) {
    if (!keypossible[i]) {
      error("expected key (1)");
      return false;
    } else if (tokens[i][0] == '_') {
      if (!current_frame) {
        error("missing data_ (unexpected data name)");
        return false;
      }

      if (i + 1 == n) {
        error("truncated");
        return false;
      }

      tolowerinplace(tokens[i]);
      current_frame->m_dict[tokens[i]].m_array = cif_detail::cif_str_array{};
      auto& cif_arr = std::get<cif_detail::cif_str_array>(
          current_frame->m_dict[tokens[i]].m_array);
      cif_arr.set_value(tokens[i + 1]);

      i++;
    } else if (strcasecmp("loop_", tokens[i]) == 0) {
      if (!current_frame) {
        error("missing data_ (unexpected loop)");
        return false;
      }

      int ncols = 0;
      int nrows = 0;
      cif_loop *loop = nullptr;

      // loop data
      loop = new cif_loop;
      current_frame->m_loops.emplace_back(loop);

      // columns
      while (++i < n && keypossible[i] && tokens[i][0] == '_') {
        tolowerinplace(tokens[i]);
        current_frame->m_dict[tokens[i]].m_array = cif_detail::cif_str_array{};
        auto& cif_arr = std::get<cif_detail::cif_str_array>(
            current_frame->m_dict[tokens[i]].m_array);
        cif_arr.set_loop(loop, ncols);

        ncols++;
      }

      if (loop) {
        // loop data
        loop->values = (const char **) &tokens[i];
        loop->ncols = ncols;
      }

      // rows
      while (i < n && !(keypossible[i] && isspecial(tokens[i]))) {
        i += ncols;

        if (i > n) {
          error("truncated loop");
          return false;
        }

        nrows++;
      }

      // loop data
      if (loop) {
        loop->nrows = nrows;
      }

      i--;

    } else if (strncasecmp("data_", tokens[i], 5) == 0) {
      auto& new_data = datablocksnew[tokens[i] + 5];
      new_data.m_data = cif_detail::cif_str_data();
      current_frame = &std::get<cif_detail::cif_str_data>(new_data.m_data);
      current_frame->m_code = tokens[i] + 5;
      frame_stack = {current_frame};

    } else if (strncasecmp("global_", tokens[i], 5) == 0) {
      // STAR feature, not supported in CIF
      auto new_data = new cif_data;
      new_data->m_data = cif_detail::cif_str_data{};
      current_frame = &std::get<cif_detail::cif_str_data>(new_data->m_data);
      global_block.reset(new_data);
      frame_stack = {current_frame};

    } else if (strncasecmp("save_", tokens[i], 5) == 0) {
      if (tokens[i][5]) {
        // begin
        if (!current_frame) {
          error("top-level save_");
          return false;
        }

        const char * key(tokens[i] + 5);
        current_frame = &current_frame->m_saveframes[key];
        frame_stack.push_back(current_frame);
      } else {
        // end
        if (frame_stack.size() < 2) {
          error("unexpected save_");
          return false;
        }

        frame_stack.pop_back();
        current_frame = frame_stack.back();
      }
    } else {
      error("expected key (2)");
      return false;
    }
  }

  m_datablocks = std::move(datablocksnew);

  return true;
}


#if !defined(_PYMOL_NO_MSGPACKC)
enum class DataTypes
{
  Int8 = 1,
  Int16 = 2,
  Int32 = 3,
  UInt8 = 4,
  UInt16 = 5,
  UInt32 = 6,
  Float32 = 32,
  Float64 = 33,
};

template <typename T>
void decodeAndPushBack(const std::vector<unsigned char>& bytes, std::size_t& i,
    std::size_t size, std::vector<CifArrayElement>& result)
{
  T value;
  std::memcpy(&value, &bytes[i], size);
  result.push_back(value);
}

static std::vector<CifArrayElement> byte_array_decode(const std::vector<unsigned char>& bytes, DataTypes dataType)
{
  std::vector<CifArrayElement> result;
  std::unordered_map<DataTypes, std::size_t> dataTypeSize = {
    {DataTypes::Int8, sizeof(std::int8_t)},
    {DataTypes::Int16, sizeof(std::int16_t)},
    {DataTypes::Int32, sizeof(std::int32_t)},
    {DataTypes::UInt8, sizeof(std::uint8_t)},
    {DataTypes::UInt16, sizeof(std::uint16_t)},
    {DataTypes::UInt32, sizeof(std::uint32_t)},
    {DataTypes::Float32, sizeof(float)},
    {DataTypes::Float64, sizeof(double)},
  };

  auto size = dataTypeSize[dataType];
  for (std::size_t i = 0; i < bytes.size(); i += size) {
    CifArrayElement valueVar;
    switch (dataType) {
    case DataTypes::Int8:
      decodeAndPushBack<std::int8_t>(bytes, i, size, result);
      break;
    case DataTypes::Int16:
      decodeAndPushBack<std::int16_t>(bytes, i, size, result);
      break;
    case DataTypes::Int32:
      decodeAndPushBack<std::int32_t>(bytes, i, size, result);
      break;
    case DataTypes::UInt8:
      decodeAndPushBack<std::uint8_t>(bytes, i, size, result);
      break;
    case DataTypes::UInt16:
      decodeAndPushBack<std::uint16_t>(bytes, i, size, result);
      break;
    case DataTypes::UInt32:
      decodeAndPushBack<std::uint32_t>(bytes, i, size, result);
      break;
    case DataTypes::Float32:
      decodeAndPushBack<float>(bytes, i, size, result);
      break;
    case DataTypes::Float64:
      decodeAndPushBack<double>(bytes, i, size, result);
      break;
    }
  }
  return result;
}

static std::vector<CifArrayElement> integer_packing_decode(
    const std::vector<CifArrayElement>& packedInts, int byteCount, int srcSize,
    bool isUnsigned)
{
  std::vector<CifArrayElement> result(srcSize);
  std::int32_t upperLimit;
  if (isUnsigned) {
    upperLimit = byteCount == 1 ? std::numeric_limits<std::uint8_t>::max()
                                : std::numeric_limits<std::uint16_t>::max();
  } else {
    upperLimit = byteCount == 1 ? std::numeric_limits<std::int8_t>::max()
                                : std::numeric_limits<std::int16_t>::max();
  }
  std::int32_t lowerLimit = -upperLimit - 1;

  auto as_int = [isUnsigned, byteCount](auto&& elem) -> std::int32_t {
    if (isUnsigned) {
      return byteCount == 1 ? static_cast<std::int32_t>(std::get<std::uint8_t>(elem))
                            : static_cast<std::int32_t>(std::get<std::uint16_t>(elem));
    } else {
      return byteCount == 1 ? static_cast<std::int32_t>(std::get<std::int8_t>(elem))
                            : static_cast<std::int32_t>(std::get<std::int16_t>(elem));
    }
  };

  auto at_limit = [isUnsigned, upperLimit, lowerLimit](std::int32_t t) -> bool {
    return isUnsigned ? (t == upperLimit)
                      : (t == upperLimit || t == lowerLimit);
  };

  for (int i = 0, j = 0; i < packedInts.size(); ++i, ++j) {
    std::int32_t value = 0;
    std::int32_t t = as_int(packedInts[i]);
    while (at_limit(t)) {
      value += t;
      t = as_int(packedInts[++i]);
    }
    value += t;
    result[j] = value;
  }
  return result;
}

static std::vector<CifArrayElement> delta_decode(
    std::vector<CifArrayElement>& data, std::int32_t origin, DataTypes srcType)
{
  std::vector<CifArrayElement> result = data;
  result[0] = origin;
  auto add_int32_t = [](auto&& a, auto&& b) -> std::int32_t {
    return std::get<std::int32_t>(a) + std::get<std::int32_t>(b);
  };
  std::inclusive_scan(result.begin(), result.end(), result.begin(), add_int32_t);
  return result;
}

static std::vector<CifArrayElement> run_length_decode(
    std::vector<CifArrayElement>& data, DataTypes srcType, int srcSize)
{
  std::vector<CifArrayElement> result;
  for (std::size_t i = 0; i < data.size(); i += 2) {
    auto item = std::get<std::int32_t>(data[i]);
    auto count = std::get<std::int32_t>(data[i + 1]);
    for (std::int32_t j = 0; j < count; j++) {
      result.push_back(item);
    }
  }
  return result;
}

static std::vector<CifArrayElement> fixed_array_decode(
    std::vector<CifArrayElement>& data, int factor, DataTypes srcType)
{
  std::vector<CifArrayElement> result = data;
  auto div_int32_t = [factor, srcType](auto&& a) -> auto {
    return srcType == DataTypes::Float32
               ? std::get<std::int32_t>(a) / static_cast<float>(factor)
               : std::get<std::int32_t>(a) / static_cast<double>(factor);
  };
  std::transform(data.begin(), data.end(), result.begin(), div_int32_t);
  return result;
}

static std::vector<CifArrayElement> interval_quant_decode(
    std::vector<CifArrayElement>& data, double min, double max, int numSteps,
    DataTypes srcType)
{
  std::vector<CifArrayElement> result = data;
  auto delta = (max - min) / (numSteps - 1);
  std::transform(data.begin(), data.end(), result.begin(),
      [min, delta](auto&& a) -> double {
        return min + std::get<std::int32_t>(a) * delta;
      });
  return result;
}

static std::vector<CifArrayElement> parse_bcif_decode(
    const std::vector<unsigned char>& rawData,
    std::vector<std::map<std::string, msgpack::object>>& dataEncoding);

static std::vector<CifArrayElement> string_array_decode(
    const std::vector<unsigned char>& data,
    std::vector<std::map<std::string, msgpack::object>>& indicesEncoding,
    const std::string& stringData, const std::vector<unsigned char>& offsets,
    std::vector<std::map<std::string, msgpack::object>>& offsetEncoding)
{
  auto decodedOffsets = parse_bcif_decode(offsets, offsetEncoding);
  auto indices = parse_bcif_decode(data, indicesEncoding);

  std::vector<CifArrayElement> result;
  result.reserve(indices.size());

  std::vector<std::string> strings = {""};
  strings.reserve(decodedOffsets.size());
  for (int i = 1; i < decodedOffsets.size(); i++) {
    auto start = std::get<std::int32_t>(decodedOffsets[i - 1]);
    auto end = std::get<std::int32_t>(decodedOffsets[i]);
    auto str = stringData.substr(start, end - start);
    strings.push_back(str);
  }

  for (int i = 0; i < indices.size(); i++) {
    auto index = std::get<std::int32_t>(indices[i]);
    result.push_back(strings[index + 1]);
  }
  return result;
}

static void parse_bcif_decode_kind(const std::string& kind,
    const std::vector<unsigned char>& rawData,
    std::vector<CifArrayElement>& result,
    std::map<std::string, msgpack::object>& dataEncoding)
{
  if (kind == "ByteArray") {
    auto type = dataEncoding["type"].as<int>();
    result = byte_array_decode(rawData, static_cast<DataTypes>(type));
  } else if (kind == "FixedPoint") {
    auto factor = dataEncoding["factor"].as<int>();
    auto srcType = dataEncoding["srcType"].as<int>();
    result = fixed_array_decode(result, factor, static_cast<DataTypes>(srcType));
  } else if (kind == "IntervalQuantization") {
    auto min = dataEncoding["min"].as<float>();
    auto max = dataEncoding["max"].as<float>();
    auto numSteps = dataEncoding["numSteps"].as<float>();
    auto srcType = dataEncoding["srcType"].as<int>();
    result = interval_quant_decode(result, min, max, numSteps, static_cast<DataTypes>(srcType));
  } else if (kind == "RunLength") {
    auto srcType = dataEncoding["srcType"].as<int>();
    auto srcSize = dataEncoding["srcSize"].as<int>();
    result = run_length_decode(result, static_cast<DataTypes>(srcType), srcSize);
  } else if (kind == "Delta") {
    auto origin = dataEncoding["origin"].as<int>();
    auto srcType = dataEncoding["srcType"].as<int>();
    result = delta_decode(result, origin, static_cast<DataTypes>(srcType));
  } else if (kind == "IntegerPacking") {
    auto byteCount = dataEncoding["byteCount"].as<int>();
    auto srcSize = dataEncoding["srcSize"].as<int>();
    auto isUnsigned = dataEncoding["isUnsigned"].as<bool>();
    result = integer_packing_decode(result, byteCount, srcSize, isUnsigned);
  } else if (kind == "StringArray") {
    auto indicesEncoding = dataEncoding["dataEncoding"].as<std::vector<std::map<std::string, msgpack::object>>>();
    auto stringData = dataEncoding["stringData"].as<std::string>();
    auto offsets = dataEncoding["offsets"].as<std::vector<unsigned char>>();
    auto offsetEncoding = dataEncoding["offsetEncoding"].as<std::vector<std::map<std::string, msgpack::object>>>();
    result = string_array_decode(rawData, indicesEncoding, stringData, offsets, offsetEncoding);
  }
}

static std::vector<CifArrayElement> parse_bcif_decode(const std::vector<unsigned char>& rawData,
    std::vector<std::map<std::string, msgpack::object>>& dataEncoding)
{
  std::vector<CifArrayElement> result;
  for (auto it = std::rbegin(dataEncoding); it != std::rend(dataEncoding); ++it) {
    auto& dataEncode = *it;
    parse_bcif_decode_kind(
        dataEncode["kind"].as<std::string>(), rawData, result, dataEncode);
  }
  return result;
}


bool cif_file::parse_bcif(const char* bytes, std::size_t size)
{
  m_datablocks.clear();
  m_tokens.clear();

  auto oh = msgpack::unpack(bytes, size);
  auto msgobj = oh.get();
  auto dict = msgobj.as<std::map<std::string, msgpack::object>>();

  auto dataBlocksRaw = dict["dataBlocks"].as<std::vector<msgpack::object>>();
  for (const auto& block : dataBlocksRaw) {
    auto blockMap = block.as<std::map<std::string, msgpack::object>>();
    auto header = blockMap["header"].as<std::string>();
    auto categoriesRaw = blockMap["categories"].as<std::vector<msgpack::object>>();
    auto& categoriesData = m_datablocks[header].m_data.emplace<pymol::cif_detail::bcif_data>();
    for (const auto& category : categoriesRaw) {
      auto categoryMap = category.as<std::map<std::string, msgpack::object>>();
      auto categoryName = categoryMap["name"].as<std::string>();
      std::transform(categoryName.begin(), categoryName.end(),
          categoryName.begin(), ::tolower);
      auto columnsRaw = categoryMap["columns"].as<std::vector<msgpack::object>>();
      auto& columns = categoriesData.m_dict[categoryName];
      for (const auto& column : columnsRaw) {
        auto columnMap = column.as<std::map<std::string, msgpack::object>>();
        auto columnName = columnMap["name"].as<std::string>();
        std::transform(columnName.begin(), columnName.end(),
          columnName.begin(), ::tolower);
        auto dataRaw = columnMap["data"].as<std::map<std::string, msgpack::object>>();
        auto dataData = dataRaw["data"].as<std::vector<unsigned char>>();
        auto dataEncoding = dataRaw["encoding"].as<std::vector<std::map<std::string, msgpack::object>>>();
        columns[columnName] = parse_bcif_decode(dataData, dataEncoding);
      }
    }
  }
  return true;
}
#else
bool cif_file::parse_bcif(const char* bytes, std::size_t size)
{
  return false;
}
#endif // !defined(_PYMOL_NO_MSGPACKC)

} // namespace pymol

// vi:sw=2:ts=2
