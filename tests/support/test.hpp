#pragma once

#include <cmath>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mps::test {

using TestFunction = void (*)();

struct TestCase {
  std::string_view name;
  TestFunction function;
};

inline std::vector<TestCase>& registry() {
  static std::vector<TestCase> tests;
  return tests;
}

class Registrar {
 public:
  Registrar(const std::string_view name, const TestFunction function) {
    registry().push_back(TestCase{name, function});
  }
};

inline int& failure_count() {
  static int count = 0;
  return count;
}

inline void record_failure(const char* file, const int line,
                           const std::string_view expression,
                           const std::string_view details = {}) {
  ++failure_count();
  std::cerr << file << ':' << line << ": check failed: " << expression;
  if (!details.empty()) {
    std::cerr << " (" << details << ')';
  }
  std::cerr << '\n';
}

template <typename Actual, typename Expected>
void check_equal(const Actual& actual, const Expected& expected, const char* file,
                 const int line, const std::string_view expression) {
  if (!(actual == expected)) {
    std::ostringstream details;
    details << "actual=" << actual << ", expected=" << expected;
    record_failure(file, line, expression, details.str());
  }
}

inline void check_near(const double actual, const double expected,
                       const double tolerance, const char* file, const int line,
                       const std::string_view expression) {
  if (!std::isfinite(actual) || !std::isfinite(expected) || tolerance < 0.0 ||
      std::abs(actual - expected) > tolerance) {
    std::ostringstream details;
    details << "actual=" << actual << ", expected=" << expected
            << ", tolerance=" << tolerance;
    record_failure(file, line, expression, details.str());
  }
}

inline int run_all() {
  const int initial_failures = failure_count();
  for (const auto& test : registry()) {
    try {
      test.function();
    } catch (const std::exception& error) {
      record_failure(__FILE__, __LINE__, test.name, error.what());
    } catch (...) {
      record_failure(__FILE__, __LINE__, test.name, "unknown exception");
    }
  }

  const int failures = failure_count() - initial_failures;
  if (failures == 0) {
    std::cout << registry().size() << " test case(s) passed\n";
    return 0;
  }
  std::cerr << failures << " check(s) failed\n";
  return 1;
}

}  // namespace mps::test

#define MPS_TEST_CONCAT_INNER(lhs, rhs) lhs##rhs
#define MPS_TEST_CONCAT(lhs, rhs) MPS_TEST_CONCAT_INNER(lhs, rhs)

#define MPS_TEST_CASE(name) MPS_TEST_CASE_IMPL(name, __LINE__)
#define MPS_TEST_CASE_IMPL(name, line)                                            \
  static void MPS_TEST_CONCAT(mps_test_function_, line)();                        \
  static const ::mps::test::Registrar MPS_TEST_CONCAT(mps_test_registrar_, line)( \
      name, &MPS_TEST_CONCAT(mps_test_function_, line));                          \
  static void MPS_TEST_CONCAT(mps_test_function_, line)()

#define MPS_CHECK(expression)                                       \
  do {                                                              \
    if (!(expression)) {                                            \
      ::mps::test::record_failure(__FILE__, __LINE__, #expression); \
    }                                                               \
  } while (false)

#define MPS_CHECK_EQ(actual, expected)                                 \
  do {                                                                 \
    ::mps::test::check_equal((actual), (expected), __FILE__, __LINE__, \
                             #actual " == " #expected);                \
  } while (false)

#define MPS_CHECK_NEAR(actual, expected, tolerance)                                \
  do {                                                                             \
    ::mps::test::check_near((actual), (expected), (tolerance), __FILE__, __LINE__, \
                            #actual " ~= " #expected);                             \
  } while (false)

#define MPS_CHECK_THROWS_AS(expression, exception_type)                    \
  do {                                                                     \
    bool mps_caught_expected_exception = false;                            \
    try {                                                                  \
      static_cast<void>(expression);                                       \
    } catch (const exception_type&) {                                      \
      mps_caught_expected_exception = true;                                \
    } catch (...) {                                                        \
    }                                                                      \
    if (!mps_caught_expected_exception) {                                  \
      ::mps::test::record_failure(__FILE__, __LINE__,                      \
                                  #expression " throws " #exception_type); \
    }                                                                      \
  } while (false)
