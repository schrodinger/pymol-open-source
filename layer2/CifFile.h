/*
 * CIF tokenizer
 *
 * (c) 2014 Schrodinger, Inc.
 */

#ifndef _H_CIFFILE
#define _H_CIFFILE

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <map>
#include <memory>
#include <vector>
#include <string>
#include <variant>

// for pymol::default_free
#include "MemoryDebug.h"

template<class... Ts>
struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

namespace pymol {
namespace _cif_detail {

/**
 * Null-terminated string view.
 */
class zstring_view {
  const char* m_data;

public:
  zstring_view(const char* s) : m_data(s) {}

  bool operator<(zstring_view rhs) const {
    return std::strcmp(m_data, rhs.m_data) < 0;
  }
};

/**
 * Convert a raw cif data value to a typed value
 */
template <typename T> T raw_to_typed(const char*);

} // namespace _cif_detail

// cif data types
class cif_data;
class cif_loop;
class cif_array;
namespace cif_detail {
  struct cif_str_data;
  struct bcif_data;
};
using CIFData = std::variant<cif_detail::cif_str_data, cif_detail::bcif_data>;

/**
 * Class for reading CIF files.
 * Parses the entire file and exposes its data blocks.
 *
 * Read CIF file:
 * @verbatim auto cf = cif_file("file.cif"); @endverbatim
 *
 * Read CIF string:
 * @verbatim auto cf = cif_file(nullptr, cifstring); @endverbatim
 *
 * Iterate over data blocks:
 * @verbatim
   for (auto& [code, block] : cf.datablocks()) {
     // data_<code>
     const char* code = block->code();

     // get data item pointer, or nullptr if name not found
     auto* dataitem1 = block->get_arr("_some.name1");
     auto* dataitem2 = block->get_arr("_some.name2", "_alternate.name");

     // get data item pointer, or default item if name not found
     auto* dataitem3 = block->get_opt("_some.name3");

     // value of data item
     const char* stringvalue = dataitem3->as_s();
     int intvalue = dataitem3->as_i();

     // values of looped data item
     for (unsigned i = 0, i_end = dataitem3->size(); i != i_end; ++i) {
       const char* stringvalue = dataitem3->as_s(i);
     }
   }
   @endverbatim
 */
class cif_file {
  std::vector<char*> m_tokens;
  std::map<std::string, cif_data> m_datablocks;
  std::unique_ptr<char, pymol::default_free> m_contents;

  /**
   * Parse CIF string
   * @param p CIF string (takes ownership)
   * @post datablocks() is valid
   */
  bool parse(char*&&);

public:
  /// Parse CIF file
  bool parse_file(const char*);

  /// Parse CIF string
  bool parse_string(const char*);

  /**
   * Parse BinaryCIF blob
   * @param bytes BinaryCIF blob
   * @param size Blob size
   * @post datablocks() is valid
  */
  bool parse_bcif(const char* bytes, std::size_t size);

protected:
  /// Report a parsing error
  virtual void error(const char*);

public:
  cif_file();
  cif_file(cif_file&&);
  cif_file(const cif_file&) = delete;
  cif_file& operator=(cif_file&&);
  cif_file& operator=(const cif_file&) = delete;
  virtual ~cif_file();

  /// Construct from file name or buffer
  cif_file(const char* filename, const char* contents = nullptr);

  /// Data blocks
  const std::map<std::string, cif_data>& datablocks() const { return m_datablocks; }
};


using CifArrayElement = std::variant<std::int8_t, std::int16_t, std::int32_t,
    std::uint8_t, std::uint16_t, std::uint32_t, float, double, std::string>;

namespace cif_detail {
  struct cif_str_array {
    enum { NOT_IN_LOOP = -1 };

    // column index, -1 if not in loop
    short col;

    // pointer to either loop or single value
    union {
      const cif_loop * loop;
      const char * value;
    } pointer;

    // Raw data value or NULL for unknown/inapplicable and `pos >= size()`
    const char* get_value_raw(unsigned pos = 0) const;

    // point this array to a loop (only for parsing)
    void set_loop(const cif_loop * loop, short col_) {
      col = col_;
      pointer.loop = loop;
    };

    // point this array to a single value (only for parsing)
    void set_value(const char * value) {
      col = NOT_IN_LOOP;
      pointer.value = value;
    };
  };
  struct bcif_array {
    std::vector<CifArrayElement> m_arr{};
  };

  /**
   * Returns a typed value from a CIF data element.
   * If the element is missing or inapplicable, return `d`.
   * @param var CIF data element
   * @param d default value
   * @return typed value
   */
  template <typename T> T var_to_typed(const CifArrayElement& var, const T& d)
  {
    if constexpr (std::is_same_v<T, const char*>) {
      auto& str = std::get<std::string>(var);
      return !str.empty() ? str.c_str() : d;
    } else {
      if (auto ptr = std::get_if<std::string>(&var); ptr && ptr->empty()) {
        return d;
      }
      if constexpr (!std::is_same_v<T, std::string>) {
        return std::visit(overloaded{[](const std::string& s) -> T {
                                       return _cif_detail::raw_to_typed<T>(
                                           s.c_str());
                                     },
                              [](const auto& v) -> T { return v; }},
            var);
      }
    }
    return d;
  }
}

/**
 * View on a CIF data array. The viewed data is owned by the cif_file
 */
class cif_array {
  friend class cif_file;

private:
  mutable std::string m_internal_str_cache;
  std::variant<cif_detail::cif_str_array, cif_detail::bcif_array> m_array;

public:
  // constructor
  cif_array() = default;

  // constructor (only needed for EMPTY_ARRAY)
  cif_array(std::nullptr_t) { 
    if (auto arr = std::get_if<cif_detail::cif_str_array>(&m_array)) {
      arr->set_value(nullptr);
    } else if (auto arr = std::get_if<cif_detail::bcif_array>(&m_array)) {
      arr->m_arr.clear();
    }
  }

  cif_array(std::vector<CifArrayElement>&& arr) {
    m_array = cif_detail::bcif_array{std::move(arr)};
  }

  /// Number of elements in this array (= number of rows in loop)
  unsigned size() const;

  /// True if value in ['.', '?']
  bool is_missing(unsigned pos = 0) const {
    if (auto arr = std::get_if<cif_detail::cif_str_array>(&m_array)) {
      return !arr->get_value_raw(pos);
    } else {
      return false;
    }
  }

  /// True if all values in ['.', '?']
  bool is_missing_all() const;

  /**
   * Get element as type T. If `pos >= size()` then return `d`.
   * @param pos element index (= row index in loop)
   * @param d default value for unknown/inapplicable elements
   */
  template <typename T> T as(unsigned pos = 0, T d = T()) const {
    if (auto arr = std::get_if<cif_detail::cif_str_array>(&m_array)) {
      const char* s = arr->get_value_raw(pos);
      return s ? _cif_detail::raw_to_typed<T>(s) : d;
    } else if (auto arr = std::get_if<cif_detail::bcif_array>(&m_array)) {
      if (pos >= arr->m_arr.size())
        return d;
      auto& var = arr->m_arr[pos];
      return cif_detail::var_to_typed<T>(var, d);
    }
    return d;
  }

  /**
   * Get element as null-terminated string. The default value is the empty
   * string, unlike as<const char*>() which returns nullptr as the default value.
   * If `pos >= size()` then return `d`.
   * @param pos element index (= row index in loop)
   * @param d default value for unknown/inapplicable elements
   */
  const char* as_s(unsigned pos = 0, const char* d = "") const {
    if (std::get_if<cif_detail::cif_str_array>(&m_array)) {
      return as(pos, d);
    } else if (auto arr = std::get_if<cif_detail::bcif_array>(&m_array)) {
      if (pos >= arr->m_arr.size())
        return d;
      if (auto str_ptr = std::get_if<std::string>(&arr->m_arr[pos])) {
        return str_ptr->c_str();
      }
      m_internal_str_cache = std::visit([](auto&& arg) -> std::string {
        if constexpr (std::is_same_v<std::decay_t<decltype(arg)>,
                          std::string>) {
          return arg;
        } else {
          return std::to_string(arg);
        }
      }, arr->m_arr[pos]);
      return m_internal_str_cache.c_str();
    }
    return d;
  }

  /// Alias for as<int>()
  int as_i(unsigned pos = 0, int d = 0) const { return as(pos, d); }

  /// Alias for as<double>()
  double as_d(unsigned pos = 0, double d = 0.) const { return as(pos, d); }

  /**
   * Get a copy of the array.
   * @param d default value for unknown/inapplicable elements
   */
  template <typename T> std::vector<T> to_vector(T d = T()) const {
    auto n = size();
    std::vector<T> v;
    v.reserve(n);
    for (unsigned i = 0; i < n; ++i)
      v.push_back(as<T>(i, d));
    return v;
  }
};

/**
 * CIF data block. The viewed data is owned by the cif_file.
 */

namespace cif_detail {
  struct cif_str_data {
    // data_<code>
    const char* m_code = nullptr;

    std::map<_cif_detail::zstring_view, cif_array> m_dict;
    std::map<std::string, cif_array> m_dict_str;
    std::map<_cif_detail::zstring_view, cif_detail::cif_str_data> m_saveframes;

    // only needed for freeing
    std::vector<std::unique_ptr<cif_loop>> m_loops;
  };

  using ColumnMap = std::map<std::string, std::vector<CifArrayElement>>;
  using CategoryMap = std::map<std::string, ColumnMap>;
  using DataBlockMap = std::map<std::string, CategoryMap>;
  struct bcif_data {
    std::string m_code;
    std::map<std::string, std::map<std::string, cif_array>> m_dict;
  };
}

class cif_data {
  friend class cif_file;

  CIFData m_data;

  // generic default value
  static const cif_array* empty_array();

public:

  cif_data() = default;
  cif_data(const cif_data&) = delete;
  cif_data(cif_data&&) = default;
  cif_data& operator=(const cif_data&) = delete;
  cif_data& operator=(cif_data&&) = default;

  /// Block code (never nullptr)
  const char* code() const;

  // Get a pointer to array or nullptr if not found
  const cif_array* get_arr(const char* key) const;
  template <typename... Args>
  const cif_array* get_arr(const char* key, Args... aliases) const
  {
    auto arr = get_arr(key);
    return arr ? arr : get_arr(aliases...);
  }

  /// Like get_arr() but return a default value instead of nullptr if not found
  template <typename... Args> const cif_array* get_opt(Args... keys) const
  {
    auto arr = get_arr(keys...);
    return arr ? arr : empty_array();
  }

  /// Get a pointer to a save frame or nullptr if not found
  const cif_detail::cif_str_data* get_saveframe(const char* code) const;
};

} // namespace pymol

#endif
// vi:sw=2:ts=2
