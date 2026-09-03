#include <exception>
#include <iostream>
#include <string_view>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/version.hpp"

namespace {

void print_usage(std::ostream& output) {
  output << "Usage: my_planet_sim [--help] [--config PATH]\n";
}

}  // namespace

int main(const int argc, const char* const argv[]) {
  static_cast<void>(mps::version());
  try {
    if (argc == 1) {
      return 0;
    }
    if (argc == 2 && std::string_view(argv[1]) == "--help") {
      print_usage(std::cout);
      return 0;
    }
    if (argc == 3 && std::string_view(argv[1]) == "--config") {
      static_cast<void>(mps::load_experiment_config(argv[2]));
      return 0;
    }
    print_usage(std::cerr);
    return 2;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
