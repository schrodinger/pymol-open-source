#pragma once

namespace pymol
{

/**
 * C++14's std::integer_sequnce & helpers
 */
template <typename IntegralT, IntegralT... Ints> struct integer_sequence {
  constexpr std::size_t size() noexcept { return sizeof...(Ints); }
};
namespace integer_sequence_detail
{
template <typename IntegerSequenceT, std::size_t NewInt>
struct integer_sequence_pushback;

template <typename IntegralT, std::size_t NewInt, IntegralT... Ints>
struct integer_sequence_pushback<integer_sequence<IntegralT, Ints...>, NewInt> {
  using type = integer_sequence<IntegralT, Ints..., sizeof...(Ints)>;
};

template <typename IntegralT, IntegralT I, IntegralT N>
struct make_integer_sequence_impl;

template <typename IntegralT, IntegralT I, IntegralT N>
struct make_integer_sequence_impl {
  using type = typename integer_sequence_pushback<
      typename make_integer_sequence_impl<IntegralT, I + 1, N>::type,
      I + 1>::type;
};

template <typename IntegralT, IntegralT N>
struct make_integer_sequence_impl<IntegralT, N, N> {
  using type = integer_sequence<IntegralT>;
};

} // namespace integer_sequence_detail

template <typename IntegralT, IntegralT N>
using make_integer_sequence =
    typename integer_sequence_detail::make_integer_sequence_impl<IntegralT,
        IntegralT{0}, N>::type;

template <std::size_t... Ints>
using index_sequence = integer_sequence<std::size_t, Ints...>;

template <std::size_t N>
using make_index_sequence = make_integer_sequence<std::size_t, N>;

template <class... T>
using index_sequence_for = make_index_sequence<sizeof...(T)>;

} // namespace pymol

