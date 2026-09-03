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

add_library(
  myplanetsim_core
  src/config/experiment_config.cpp
  src/core/planet_parameters.cpp
  src/core/validation.cpp
  src/core/version.cpp
  src/diagnostics/reductions.cpp
)
target_compile_features(myplanetsim_core PUBLIC cxx_std_20)
target_include_directories(
  myplanetsim_core
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
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
  add_test(NAME smoke.empty_program COMMAND my_planet_sim)
  add_test(NAME smoke.help COMMAND my_planet_sim --help)
  add_test(
    NAME smoke.config
    COMMAND my_planet_sim --config ${CMAKE_CURRENT_SOURCE_DIR}/configs/phase0_ode.cfg
  )
  add_subdirectory(tests)
endif()
