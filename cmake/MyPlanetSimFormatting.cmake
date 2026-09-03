find_program(MPS_CLANG_FORMAT_EXECUTABLE NAMES clang-format)

file(
  GLOB_RECURSE MPS_FORMAT_SOURCES
  CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/apps/*.cpp"
  "${PROJECT_SOURCE_DIR}/include/*.hpp"
  "${PROJECT_SOURCE_DIR}/src/*.cpp"
  "${PROJECT_SOURCE_DIR}/tests/*.cpp"
  "${PROJECT_SOURCE_DIR}/tests/*.hpp"
)

if(MPS_CLANG_FORMAT_EXECUTABLE)
  add_custom_target(
    format-check
    COMMAND
      "${MPS_CLANG_FORMAT_EXECUTABLE}" --dry-run --Werror ${MPS_FORMAT_SOURCES}
    COMMENT "Checking C++ formatting with ${MPS_CLANG_FORMAT_EXECUTABLE}"
    VERBATIM
  )
else()
  add_custom_target(
    format-check
    COMMAND
      "${CMAKE_COMMAND}" -E echo
      "clang-format was not found; install it to run format-check"
    COMMAND "${CMAKE_COMMAND}" -E false
    VERBATIM
  )
endif()
