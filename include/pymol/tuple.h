#pragma once

#include "pymol/functional.h"
#include "pymol/type_traits.h"
#include "pymol/utility.h"

namespace pymol
{

/*
 * C++17's std::apply
 * requires pymol::(or std::)index_sequence
 * Takes a function and a tuple of lvalue/rvalue references to arguments
 * and applys them to a function
 */

namespace apply_detail
{
template <typename Func, typename TupleArgs, std::size_t... I>
/*decltype(auto)*/ auto apply_impl(Func&& f, TupleArgs&& t,
    index_sequence<I...>) -> decltype(pymol::invoke(std::forward<Func>(f),
    std::get<I>(std::forward<TupleArgs>(t))...))
{
  return pymol::invoke(
      std::forward<Func>(f), std::get<I>(std::forward<TupleArgs>(t))...);
}
} // namespace apply_detail

template <class F, class TupleArgs>
/*decltype(auto)*/ auto apply(F&& f, TupleArgs&& t) -> decltype(
    apply_detail::apply_impl(std::forward<F>(f), std::forward<TupleArgs>(t),
        make_index_sequence<
            std::tuple_size<remove_reference_t<TupleArgs>>::value>{}))
{
  return apply_detail::apply_impl(std::forward<F>(f),
      std::forward<TupleArgs>(t),
      make_index_sequence<
          std::tuple_size<remove_reference_t<TupleArgs>>::value>{});
}

} // namespace pymol
