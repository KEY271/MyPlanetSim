cmake_minimum_required(VERSION 3.25)

project(
  MyPlanetSim
  VERSION 0.1.0
  DESCRIPTION "A staged global atmospheric simulator for fictional planets"
  LANGUAGES CXX
)

include(cmake/MyPlanetSimWarnings.cmake)
include(cmake/MyPlanetSimSanitizers.cmake)
include(cmake/MyPlanetSimFormatting.cmake)

find_package(Git QUIET)
set(MPS_GIT_REVISION "unknown")
set(MPS_GIT_DIRTY_VALUE 0)
if(Git_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse --short=12 HEAD
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    RESULT_VARIABLE MPS_GIT_REVISION_RESULT
    OUTPUT_VARIABLE MPS_GIT_REVISION_OUTPUT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(MPS_GIT_REVISION_RESULT EQUAL 0)
    set(MPS_GIT_REVISION "${MPS_GIT_REVISION_OUTPUT}")
  endif()
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" status --porcelain
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE MPS_GIT_STATUS
    ERROR_QUIET
  )
  if(NOT MPS_GIT_STATUS STREQUAL "")
    set(MPS_GIT_DIRTY_VALUE 1)
  endif()
endif()

configure_file(
  cmake/build_info.hpp.in
  "${CMAKE_CURRENT_BINARY_DIR}/generated/myplanetsim/build_info.hpp"
  @ONLY
)

add_library(
  myplanetsim_core
  src/config/experiment_config.cpp
  src/core/planet_parameters.cpp
  src/core/validation.cpp
  src/core/version.cpp
  src/diagnostics/reductions.cpp
  src/geometry/cubed_sphere.cpp
  src/grid/cubed_sphere_field.cpp
  src/grid/cubed_sphere_grid.cpp
  src/grid/panel_topology.cpp
  src/io/checkpoint.cpp
  src/io/run_metadata.cpp
  src/phase0/ode_experiment.cpp
  src/numerics/spherical_operators.cpp
  src/transport/spherical_transport.cpp
)
target_compile_features(myplanetsim_core PUBLIC cxx_std_20)
target_include_directories(
  myplanetsim_core
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
)
target_link_libraries(
  myplanetsim_core
  PRIVATE myplanetsim_warnings myplanetsim_sanitizers
)

add_executable(my_planet_sim apps/my_planet_sim.cpp)
target_link_libraries(
  my_planet_sim
  PRIVATE myplanetsim_core myplanetsim_warnings myplanetsim_sanitizers
)

include(CTest)
if(BUILD_TESTING)
  add_test(NAME smoke.help COMMAND my_planet_sim --help)
  add_test(
    NAME smoke.phase0_ode
    COMMAND my_planet_sim --config ${CMAKE_CURRENT_SOURCE_DIR}/configs/phase0_ode.cfg
  )
  add_test(
    NAME smoke.phase1_transport
    COMMAND my_planet_sim --config ${CMAKE_CURRENT_SOURCE_DIR}/configs/phase1_solid_body.cfg
  )
  add_test(
    NAME smoke.phase1_deformational
    COMMAND my_planet_sim --config ${CMAKE_CURRENT_SOURCE_DIR}/configs/phase1_deformational.cfg
  )
  set_tests_properties(
    smoke.phase1_deformational
    PROPERTIES PASS_REGULAR_EXPRESSION "result.status = complete"
  )
  set_tests_properties(
    smoke.phase1_transport
    PROPERTIES PASS_REGULAR_EXPRESSION "result.status = complete"
  )
  set_tests_properties(
    smoke.phase0_ode
    PROPERTIES PASS_REGULAR_EXPRESSION "result.status = complete"
  )
  add_test(NAME smoke.invalid_cli COMMAND my_planet_sim --unknown value)
  set_tests_properties(smoke.invalid_cli PROPERTIES WILL_FAIL TRUE)
  add_subdirectory(tests)
endif()
