#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "myplanetsim/core/types.hpp"

namespace mps {

template <typename T>
class Field2D {
 public:
  Field2D(const Index nx, const Index ny, const Index halo = 0)
      : nx_(nx), ny_(ny), halo_(halo) {
    if (nx <= 0 || ny <= 0) {
      throw std::invalid_argument("Field2D interior extents must be positive");
    }
    if (halo < 0) {
      throw std::invalid_argument("Field2D halo width must be non-negative");
    }

    const auto index_max = std::numeric_limits<Index>::max();
    if (halo > (index_max - nx) / 2 || halo > (index_max - ny) / 2) {
      throw std::length_error("Field2D extent overflows Index");
    }
    storage_nx_ = nx + 2 * halo;
    storage_ny_ = ny + 2 * halo;

    const auto width = static_cast<std::size_t>(storage_nx_);
    const auto height = static_cast<std::size_t>(storage_ny_);
    if (height > std::numeric_limits<std::size_t>::max() / width) {
      throw std::length_error("Field2D storage size overflows size_t");
    }
    data_.resize(width * height);
  }

  [[nodiscard]] Index nx() const noexcept { return nx_; }
  [[nodiscard]] Index ny() const noexcept { return ny_; }
  [[nodiscard]] Index halo() const noexcept { return halo_; }
  [[nodiscard]] Index storage_nx() const noexcept { return storage_nx_; }
  [[nodiscard]] Index storage_ny() const noexcept { return storage_ny_; }
  [[nodiscard]] std::size_t size() const noexcept { return data_.size(); }

  [[nodiscard]] T& at(const Index i, const Index j) {
    check_index(i, j);
    return data_[linear_index_unchecked(i, j)];
  }

  [[nodiscard]] const T& at(const Index i, const Index j) const {
    check_index(i, j);
    return data_[linear_index_unchecked(i, j)];
  }

  [[nodiscard]] T& operator()(const Index i, const Index j) noexcept {
    assert(contains(i, j));
    return data_[linear_index_unchecked(i, j)];
  }

  [[nodiscard]] const T& operator()(const Index i, const Index j) const noexcept {
    assert(contains(i, j));
    return data_[linear_index_unchecked(i, j)];
  }

  [[nodiscard]] std::span<T> storage() noexcept { return data_; }
  [[nodiscard]] std::span<const T> storage() const noexcept { return data_; }

  [[nodiscard]] std::span<T> interior_row(const Index j) {
    if (j < 0 || j >= ny_) {
      throw std::out_of_range("Field2D interior row is out of range");
    }
    return std::span<T>(&data_[linear_index_unchecked(0, j)],
                        static_cast<std::size_t>(nx_));
  }

  [[nodiscard]] std::span<const T> interior_row(const Index j) const {
    if (j < 0 || j >= ny_) {
      throw std::out_of_range("Field2D interior row is out of range");
    }
    return std::span<const T>(&data_[linear_index_unchecked(0, j)],
                              static_cast<std::size_t>(nx_));
  }

  void fill(const T& value) { std::fill(data_.begin(), data_.end(), value); }

 private:
  [[nodiscard]] bool contains(const Index i, const Index j) const noexcept {
    return i >= -halo_ && i < nx_ + halo_ && j >= -halo_ && j < ny_ + halo_;
  }

  void check_index(const Index i, const Index j) const {
    if (!contains(i, j)) {
      throw std::out_of_range("Field2D index (" + std::to_string(i) + ", " +
                              std::to_string(j) + ") is outside storage");
    }
  }

  [[nodiscard]] std::size_t linear_index_unchecked(const Index i,
                                                   const Index j) const noexcept {
    const auto storage_i = static_cast<std::size_t>(i + halo_);
    const auto storage_j = static_cast<std::size_t>(j + halo_);
    return storage_j * static_cast<std::size_t>(storage_nx_) + storage_i;
  }

  Index nx_;
  Index ny_;
  Index halo_;
  Index storage_nx_;
  Index storage_ny_;
  std::vector<T> data_;
};

}  // namespace mps
