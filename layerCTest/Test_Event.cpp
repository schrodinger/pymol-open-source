#include "Test.h"

#include "Event.h"

static int x = 2;

static void some_func(int inc)
{
  x += inc;
  x += 100;
}

TEST_CASE("Event Test", "[Event]")
{
  auto callback = [&](int inc) {
    x += inc;
    x += 10;
  };
  pymol::Event<int> event_publisher{};
  REQUIRE(event_publisher.size() == 0);
  event_publisher.add_listener(callback);
  REQUIRE(event_publisher.size() == 1);
  event_publisher.invoke(5); // 2 + 5 + 10
  REQUIRE(x == 17);
  event_publisher.add_listener(some_func);
  event_publisher.invoke(-5);
  REQUIRE(x == 117); // 17 + (- 5 + 10) + (- 5 + 100)
}
