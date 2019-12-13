#include <iostream>
#include "share/scream_assert.hpp"
#include "share/pio/scorpio/scream_scorpio_utils.hpp"

namespace scream {
namespace pio {

void hello_world () {

  std::cout << "Hello, hello, hello, is there anybody out there?\n" << std::flush;

}

void scorpio_createfile() {
  scream_scorpio_createfile();
}




} // namespace pio
} // namespace scream
