#pragma once

#include <cstddef>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace mps {

template <typename T>
class ColumnView {
 public:
  using value_type = std::remove_cv_t<T>;

  constexpr ColumnView() = default;
  constexpr ColumnView(T* data, const std::size_t levels,
                       const std::ptrdiff_t level_stride) noexcept
      : data_(data), levels_(levels), level_stride_(level_stride) {}

  [[nodiscard]] constexpr std::size_t size() const noexcept { return levels_; }
  [[nodiscard]] constexpr std::ptrdiff_t stride() const noexcept {
    return level_stride_;
  }
  [[nodiscard]] constexpr T& operator[](const std::size_t level) const noexcept {
    return data_[static_cast<std::ptrdiff_t>(level) * level_stride_];
  }
  [[nodiscard]] T& at(const std::size_t level) const {
    if (level >= levels_) throw std::out_of_range("column level is out of range");
    return (*this)[level];
  }

 private:
  T* data_ = nullptr;
  std::size_t levels_ = 0;
  std::ptrdiff_t level_stride_ = 1;
};

template <typename T>
class Field3DView {
 public:
  using value_type = std::remove_cv_t<T>;

  constexpr Field3DView() = default;
  Field3DView(std::span<T> storage, const std::size_t components,
              const std::size_t cells, const std::size_t levels,
              const std::ptrdiff_t component_stride, const std::ptrdiff_t cell_stride,
              const std::ptrdiff_t level_stride)
      : storage_(storage),
        components_(components),
        cells_(cells),
        levels_(levels),
        component_stride_(component_stride),
        cell_stride_(cell_stride),
        level_stride_(level_stride) {
    if (component_stride < 0 || cell_stride < 0 || level_stride < 0)
      throw std::invalid_argument("field strides must be non-negative");
    if (components_ != 0 && cells_ != 0 && levels_ != 0 &&
        offset(components_ - 1, cells_ - 1, levels_ - 1) >= storage_.size())
      throw std::invalid_argument("field storage is smaller than its view");
  }

  [[nodiscard]] constexpr std::size_t components() const noexcept {
    return components_;
  }
  [[nodiscard]] constexpr std::size_t cells() const noexcept { return cells_; }
  [[nodiscard]] constexpr std::size_t levels() const noexcept { return levels_; }
  [[nodiscard]] constexpr std::ptrdiff_t component_stride() const noexcept {
    return component_stride_;
  }
  [[nodiscard]] constexpr std::ptrdiff_t cell_stride() const noexcept {
    return cell_stride_;
  }
  [[nodiscard]] constexpr std::ptrdiff_t level_stride() const noexcept {
    return level_stride_;
  }

  [[nodiscard]] T& operator()(const std::size_t component, const std::size_t cell,
                              const std::size_t level) const noexcept {
    return storage_[offset(component, cell, level)];
  }
  [[nodiscard]] T& at(const std::size_t component, const std::size_t cell,
                      const std::size_t level) const {
    if (component >= components_ || cell >= cells_ || level >= levels_)
      throw std::out_of_range("field index is out of range");
    return (*this)(component, cell, level);
  }
  [[nodiscard]] ColumnView<T> column(const std::size_t component,
                                     const std::size_t cell) const {
    if (component >= components_ || cell >= cells_)
      throw std::out_of_range("field column is out of range");
    return {storage_.data() + offset(component, cell, 0), levels_, level_stride_};
  }

 private:
  [[nodiscard]] constexpr std::size_t offset(const std::size_t component,
                                             const std::size_t cell,
                                             const std::size_t level) const noexcept {
    return static_cast<std::size_t>(static_cast<std::ptrdiff_t>(component) *
                                        component_stride_ +
                                    static_cast<std::ptrdiff_t>(cell) * cell_stride_ +
                                    static_cast<std::ptrdiff_t>(level) * level_stride_);
  }

  std::span<T> storage_{};
  std::size_t components_ = 0;
  std::size_t cells_ = 0;
  std::size_t levels_ = 0;
  std::ptrdiff_t component_stride_ = 0;
  std::ptrdiff_t cell_stride_ = 0;
  std::ptrdiff_t level_stride_ = 0;
};

template <typename T>
[[nodiscard]] Field3DView<T> make_cell_column_field_view(const std::span<T> storage,
                                                         const std::size_t components,
                                                         const std::size_t cells,
                                                         const std::size_t levels) {
  return {storage,
          components,
          cells,
          levels,
          static_cast<std::ptrdiff_t>(cells * levels),
          static_cast<std::ptrdiff_t>(levels),
          1};
}

template <typename T>
class BlockedField3DView {
 public:
  BlockedField3DView(std::span<T> storage, const std::size_t components,
                     const std::size_t cells, const std::size_t levels,
                     const std::size_t lanes)
      : storage_(storage),
        components_(components),
        cells_(cells),
        levels_(levels),
        lanes_(lanes),
        blocks_(lanes == 0 ? 0 : (cells + lanes - 1) / lanes) {
    if (lanes_ == 0)
      throw std::invalid_argument("blocked field lanes must be positive");
    if (storage_.size() < required_size())
      throw std::invalid_argument("blocked field storage is smaller than its view");
  }

  [[nodiscard]] std::size_t required_size() const noexcept {
    return components_ * blocks_ * levels_ * lanes_;
  }
  [[nodiscard]] std::size_t padded_cells() const noexcept { return blocks_ * lanes_; }
  [[nodiscard]] T& operator()(const std::size_t component, const std::size_t cell,
                              const std::size_t level) const noexcept {
    const auto block = cell / lanes_;
    const auto lane = cell % lanes_;
    return storage_[((component * blocks_ + block) * levels_ + level) * lanes_ + lane];
  }

 private:
  std::span<T> storage_{};
  std::size_t components_ = 0;
  std::size_t cells_ = 0;
  std::size_t levels_ = 0;
  std::size_t lanes_ = 0;
  std::size_t blocks_ = 0;
};

}  // namespace mps
