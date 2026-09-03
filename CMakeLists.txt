cmake_minimum_required(VERSION 3.25)

project(
  MyPlanetSim
  VERSION 0.1.0
  DESCRIPTION "A staged global atmospheric simulator for fictional planets"
  LANGUAGES CXX
)

add_library(myplanetsim_core src/core/version.cpp)
target_compile_features(myplanetsim_core PUBLIC cxx_std_20)
target_include_directories(
  myplanetsim_core
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

add_executable(my_planet_sim apps/my_planet_sim.cpp)
target_link_libraries(my_planet_sim PRIVATE myplanetsim_core)

if(MSVC)
  target_compile_options(myplanetsim_core PRIVATE /W4 /permissive-)
  target_compile_options(my_planet_sim PRIVATE /W4 /permissive-)
else()
  target_compile_options(myplanetsim_core PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(my_planet_sim PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(CTest)
if(BUILD_TESTING)
  add_test(NAME smoke.empty_program COMMAND my_planet_sim)
endif()
