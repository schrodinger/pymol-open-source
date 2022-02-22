#pragma once

#include <functional>
#include <vector>
#include "pymol/functional.h"

namespace pymol
{
template <typename... ParamTs> class Event
{
public:
  /**
   * Invokes all listener callbacks with parameters.
   * @param params parameters for callback
   */
  void invoke(ParamTs... params)
  {
    if (m_callbacks.empty()) {
      return;
    }
    for (const auto& callback : m_callbacks) {
      pymol::invoke(callback, params...);
    }
  }

  /**
   * @func function to be called when event occurs
   */
  void add_listener(std::function<void(ParamTs...)> func)
  {
    m_callbacks.push_back(std::move(func));
  }

  /**
   * Removes all listeners
   */
  void reset() { m_callbacks.clear(); }

  /**
   * @return the number of listeners to this event
   */
  std::size_t size() const { return m_callbacks.size(); }

private:
  std::vector<std::function<void(ParamTs...)>> m_callbacks;
};
} // namespace pymol
