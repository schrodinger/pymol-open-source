#pragma once

namespace pymol
{

/**
 * C++17's invoke
 * @param f Function to be invoked
 * @param args arguments to supply to function
 * Note: Does not return reference until C++14
 */

template <typename Func, typename... FuncArgs>
/*decltype(auto)*/ auto invoke(Func&& f, FuncArgs&&... args)
    -> decltype(std::forward<Func>(f)(std::forward<FuncArgs>(args)...))
{
  return std::forward<Func>(f)(std::forward<FuncArgs>(args)...);
}

} // namespace pymol
