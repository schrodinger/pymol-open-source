#include "Result.h"
#include "Test.h"
using namespace pymol::test;

static
pymol::Result<int> sum(int i, int j)
{
  return i + j;
}

TEST_CASE("Result simple", "[Result]")
{
  auto ans = sum(3, 5);
  static_assert(
      std::is_same<int, decltype(ans)::type>::value, "Not correct result type");
  static_assert(std::is_same<decltype(ans), pymol::Result<int>>::value,
      "Not correct result type");
  REQUIRE(ans);
  REQUIRE(ans.result() == 8);
}

template <typename T, typename U>
static
pymol::Result<pymol::common_type_t<T, U>> sum(T i, U j)
{
  return i + j;
}

TEST_CASE("Result template", "[Result]")
{
  auto ans = sum(3, 5.);
  static_assert(std::is_same<double, decltype(ans)::type>::value,
      "Not correct result type");
  static_assert(std::is_same<decltype(ans), pymol::Result<double>>::value,
      "Not correct result type");
  REQUIRE(ans);
  REQUIRE(pymol::almost_equal(ans.result(), 8.));
}

static
pymol::Result<int> sumError(int i, int j)
{
  return pymol::Error{"Values cannot be summed."};
}

TEST_CASE("Err Message", "[Result]")
{
  auto ans = sumError(4, 5);
  static_assert(
      std::is_same<int, decltype(ans)::type>::value, "Not correct result type");
  static_assert(std::is_same<decltype(ans), pymol::Result<int>>::value,
      "Not correct result type");
  REQUIRE(!ans);
  REQUIRE(ans.error().what() == "Values cannot be summed.");
}

static pymol::Result<> returnEmptyBrackets()
{
  return {};
}

TEST_CASE("Default constructor", "[Result]")
{
  auto ans = returnEmptyBrackets();
  REQUIRE(ans);
}

TEST_CASE("Error codes", "[Result]")
{
  auto error1 = pymol::Error();
  REQUIRE(error1.code() == pymol::Error::DEFAULT);
  auto error2 = pymol::Error::make<pymol::Error::QUIET>();
  REQUIRE(error2.code() == pymol::Error::QUIET);
}

TEST_CASE("Result from error code", "[Result]")
{
  // error
  pymol::Result<int> result1 = pymol::Error::QUIET;
  REQUIRE(!bool(result1));
  REQUIRE(result1.error().code() == pymol::Error::QUIET);

  // no error
  pymol::Result<int> result2 = static_cast<int>(pymol::Error::QUIET);
  REQUIRE(bool(result2));
  REQUIRE(result2.result() == static_cast<int>(pymol::Error::QUIET));
}
