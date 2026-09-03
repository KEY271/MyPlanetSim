#include "support/test.hpp"

MPS_TEST_CASE("intentional failure") { MPS_CHECK(false); }

int main() { return mps::test::run_all(); }
