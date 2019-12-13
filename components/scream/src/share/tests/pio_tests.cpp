#include <catch2/catch.hpp>

#include "share/scream_config.hpp"
#include "share/scream_pack.hpp"
#include "share/util/scream_utils.hpp"
#include "share/pio/scorpio/scream_scorpio_utils.hpp"


namespace {

TEST_CASE("pio_test", "test_createfile") {

  scream::pio::hello_world();

  scream::pio::scorpio_createfile();


} // TEST_CASE
} // empty namespace
