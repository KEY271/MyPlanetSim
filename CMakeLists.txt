cmake_minimum_required(VERSION 3.25)

project(
  MyPlanetSim
  VERSION 0.1.0
  DESCRIPTION "A staged global atmospheric simulator for fictional planets"
  LANGUAGES CXX
)

add_executable(my_planet_sim src/main.cpp)
target_compile_features(my_planet_sim PRIVATE cxx_std_20)

if(MSVC)
  target_compile_options(my_planet_sim PRIVATE /W4 /permissive-)
else()
  target_compile_options(my_planet_sim PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(CTest)
if(BUILD_TESTING)
  add_test(NAME smoke.empty_program COMMAND my_planet_sim)
endif()
