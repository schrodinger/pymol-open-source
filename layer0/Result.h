#pragma once

#include "Util2.h"
#include <string>
#include <type_traits>

namespace pymol
{
struct Void {
};

/**
 * Expressive error handling alternative
 */
class Error
{
public:
  Error() {}

  /**
   * Stores Error message
   * @param ts collection of printable values to be joined
   * into an error message
   */
  template <typename... PrintableTs>
  Error(PrintableTs&&... ts)
      : m_errMsg{pymol::join_to_string(std::forward<PrintableTs>(ts)...)}
  {
  }

  /**
   * Retrieves error message
   * @return error message
   */

  const std::string& what() const noexcept { return m_errMsg; }

private:
  std::string m_errMsg;
};

/**
 * Class that expresses the expected result of a function
 */

template <typename ResultT=Void> class Result
{

public:
  using type = ResultT;

  Result() = default;

  /**
   * Constructor alternative that allows for convertible types
   * param r result returned from function with a type convertible to ResultT
   */

  template <typename FwdResultT> Result(FwdResultT&& r)
  {
    m_result = std::forward<FwdResultT>(r);
  }

  /**
   * Constructor alternative that takes in pymol::Error. Value of expected type
   * should not be taken at this point.
   * @param e error object to express why value should not be used
   */

  Result(Error&& e) : m_error{std::move(e)}, m_valid{false} {}

  /**
   * Determines whether the value of the expected type can be used.
   */

  explicit operator bool() const noexcept { return m_valid; }

  /**
   * Retrieves the underlying error object
   */

  const Error& error() const noexcept { return m_error; }

  /**
   * Retrieves the value of the expected object
   */

  ResultT& result() { return m_result; }

private:
  ResultT m_result;
  Error m_error;
  bool m_valid{true};
};
} // namespace pymol

