#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "myplanetsim/core/field_views.hpp"

namespace {

using Clock = std::chrono::steady_clock;
volatile double benchmark_sink = 0.0;

struct Options {
  std::size_t cells = 6 * 48 * 48;
  std::size_t levels = 40;
  std::size_t components = 4;
  std::size_t repeats = 5;
};

[[nodiscard]] Options parse_options(const int argc, char** argv) {
  Options result{};
  for (int i = 1; i < argc; i += 2) {
    if (i + 1 >= argc)
      throw std::invalid_argument(
          "usage: benchmark_field_layout [--cells C] [--levels K] "
          "[--components Q] [--repeats R]");
    const auto key = std::string_view(argv[i]);
    const auto value = static_cast<std::size_t>(std::stoull(argv[i + 1]));
    if (key == "--cells")
      result.cells = value;
    else if (key == "--levels")
      result.levels = value;
    else if (key == "--components")
      result.components = value;
    else if (key == "--repeats")
      result.repeats = value;
    else
      throw std::invalid_argument("unknown benchmark_field_layout option");
  }
  if (result.cells == 0 || result.levels == 0 || result.components == 0 ||
      result.repeats == 0)
    throw std::invalid_argument("field layout benchmark sizes must be positive");
  return result;
}

template <typename Function>
[[nodiscard]] double time_seconds(Function&& function) {
  const auto start = Clock::now();
  benchmark_sink = function();
  return std::chrono::duration<double>(Clock::now() - start).count();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto options = parse_options(argc, argv);
    const auto logical_size = options.components * options.cells * options.levels;
    std::vector<double> canonical(logical_size);
    for (std::size_t i = 0; i < canonical.size(); ++i)
      canonical[i] = std::sin(static_cast<double>(i) * 1.0e-4);
    const auto current = mps::make_cell_column_field_view<const double>(
        std::span<const double>(canonical), options.components, options.cells,
        options.levels);

    std::cout << "cells=" << options.cells << '\n'
              << "levels=" << options.levels << '\n'
              << "components=" << options.components << '\n'
              << "repeats=" << options.repeats << '\n';
    const auto current_horizontal_s = time_seconds([&] {
      double sum = 0.0;
      for (std::size_t repeat = 0; repeat < options.repeats; ++repeat)
        for (std::size_t component = 0; component < options.components; ++component)
          for (std::size_t level = 0; level < options.levels; ++level)
            for (std::size_t cell = 0; cell < options.cells; ++cell)
              sum += current(component, cell, level);
      return sum;
    });
    const auto current_column_s = time_seconds([&] {
      double sum = 0.0;
      for (std::size_t repeat = 0; repeat < options.repeats; ++repeat)
        for (std::size_t component = 0; component < options.components; ++component)
          for (std::size_t cell = 0; cell < options.cells; ++cell) {
            const auto column = current.column(component, cell);
            for (std::size_t level = 0; level < options.levels; ++level)
              sum += column[level];
          }
      return sum;
    });
    std::cout << "current_horizontal_s=" << current_horizontal_s << '\n'
              << "current_column_s=" << current_column_s << '\n';

    for (const std::size_t lanes : {4U, 8U, 16U}) {
      const auto blocks = (options.cells + lanes - 1) / lanes;
      std::vector<double> storage(options.components * blocks * options.levels * lanes);
      mps::BlockedField3DView<double> blocked(storage, options.components,
                                              options.cells, options.levels, lanes);
      const auto pack_s = time_seconds([&] {
        for (std::size_t component = 0; component < options.components; ++component)
          for (std::size_t cell = 0; cell < options.cells; ++cell)
            for (std::size_t level = 0; level < options.levels; ++level)
              blocked(component, cell, level) = current(component, cell, level);
        return storage.front();
      });
      const auto horizontal_s = time_seconds([&] {
        double sum = 0.0;
        for (std::size_t repeat = 0; repeat < options.repeats; ++repeat)
          for (std::size_t component = 0; component < options.components; ++component)
            for (std::size_t block = 0; block < blocks; ++block)
              for (std::size_t level = 0; level < options.levels; ++level)
                for (std::size_t lane = 0; lane < lanes; ++lane) {
                  const auto cell = block * lanes + lane;
                  if (cell < options.cells) sum += blocked(component, cell, level);
                }
        return sum;
      });
      const auto column_s = time_seconds([&] {
        double sum = 0.0;
        for (std::size_t repeat = 0; repeat < options.repeats; ++repeat)
          for (std::size_t component = 0; component < options.components; ++component)
            for (std::size_t cell = 0; cell < options.cells; ++cell)
              for (std::size_t level = 0; level < options.levels; ++level)
                sum += blocked(component, cell, level);
        return sum;
      });
      std::cout << "blocked_lanes=" << lanes << " pack_s=" << pack_s
                << " horizontal_s=" << horizontal_s << " column_s=" << column_s
                << " padded_cells=" << blocked.padded_cells() << '\n';
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
