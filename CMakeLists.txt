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
  src/control/control_request.cpp
  src/core/planet_parameters.cpp
  src/core/orbit_parameters.cpp
  src/core/validation.cpp
  src/core/version.cpp
  src/diagnostics/reductions.cpp
  src/diagnostics/climate_statistics.cpp
  src/diagnostics/shallow_water_diagnostics.cpp
  src/diagnostics/vertical_column_diagnostics.cpp
  src/diagnostics/dry_hydrostatic_diagnostics.cpp
  src/diagnostics/planetary_scales.cpp
  src/dynamics/shallow_water_state.cpp
  src/dynamics/dry_hydrostatic_state.cpp
  src/dynamics/dry_hydrostatic_flux.cpp
  src/dynamics/dry_hydrostatic_sources.cpp
  src/dynamics/dry_hydrostatic_coupling.cpp
  src/dynamics/dry_hydrostatic_reconstruction.cpp
  src/dynamics/dry_hydrostatic_diffusion.cpp
  src/dynamics/dry_hydrostatic_driver.cpp
  src/dynamics/dry_hydrostatic_benchmarks.cpp
  src/dynamics/surface_orography.cpp
  src/dynamics/surface_boundary.cpp
  src/dynamics/vertical_column_state.cpp
  src/dynamics/vertical_column_driver.cpp
  src/dynamics/shallow_water_compatible.cpp
  src/dynamics/shallow_water_diffusion.cpp
  src/dynamics/shallow_water_initial_edits.cpp
  src/dynamics/shallow_water_flux.cpp
  src/dynamics/shallow_water_driver.cpp
  src/dynamics/shallow_water_benchmarks.cpp
  src/dynamics/shallow_water_rhs.cpp
  src/dynamics/shallow_water_reconstruction.cpp
  src/geometry/cubed_sphere.cpp
  src/grid/cubed_sphere_grid.cpp
  src/grid/dual_topology.cpp
  src/io/checkpoint.cpp
  src/io/frame.cpp
  src/io/run_metadata.cpp
  src/phase0/ode_experiment.cpp
  src/physics/held_suarez.cpp
  src/physics/planetary_newtonian.cpp
  src/physics/surface_energy_balance.cpp
  src/numerics/compatible_operators.cpp
  src/numerics/spherical_operators.cpp
  src/transport/spherical_transport.cpp
  src/thermodynamics/dry_thermodynamics.cpp
  src/vertical/hybrid_pressure_coordinate.cpp
  src/vertical/hydrostatic_column.cpp
  src/vertical/vertical_transport.cpp
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

add_executable(benchmark_dry_core tools/benchmark_dry_core.cpp)
target_compile_definitions(
  benchmark_dry_core
  PRIVATE
    MPS_DRY_CORE_BENCHMARK_CONFIG="${CMAKE_CURRENT_SOURCE_DIR}/configs/phase7_held_suarez.cfg"
)
target_link_libraries(
  benchmark_dry_core
  PRIVATE myplanetsim_core myplanetsim_warnings myplanetsim_sanitizers
)
set_target_properties(
  benchmark_dry_core PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools"
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
    NAME smoke.phase2_shallow_water
    COMMAND my_planet_sim --config ${CMAKE_CURRENT_SOURCE_DIR}/configs/phase2_rest.cfg
  )
  add_test(
    NAME smoke.phase4_vertical_column
    COMMAND my_planet_sim --config ${CMAKE_CURRENT_SOURCE_DIR}/configs/phase4_isothermal.cfg
  )
  add_test(
    NAME smoke.phase4_dry_adiabatic
    COMMAND my_planet_sim --config
            ${CMAKE_CURRENT_SOURCE_DIR}/configs/phase4_dry_adiabatic.cfg
  )
  add_test(
    NAME smoke.phase4_moving_surface_pressure
    COMMAND my_planet_sim --config
            ${CMAKE_CURRENT_SOURCE_DIR}/configs/phase4_moving_surface_pressure.cfg
  )
  add_test(
    NAME smoke.phase4_manufactured_transport
    COMMAND my_planet_sim --config
            ${CMAKE_CURRENT_SOURCE_DIR}/configs/phase4_manufactured_transport.cfg
  )
  add_test(
    NAME smoke.phase1_deformational
    COMMAND my_planet_sim --config ${CMAKE_CURRENT_SOURCE_DIR}/configs/phase1_deformational.cfg
  )
  set_tests_properties(
    smoke.phase1_deformational
    PROPERTIES PASS_REGULAR_EXPRESSION "result.status = complete"
  )
  add_test(
    NAME phase4.csv_output
    COMMAND
      ${CMAKE_COMMAND} -DMPS_EXECUTABLE=$<TARGET_FILE:my_planet_sim>
      -DMPS_CONFIG=${CMAKE_CURRENT_SOURCE_DIR}/configs/phase4_manufactured_transport.cfg
      -DMPS_WORKING_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR} -P
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_phase4_csv.cmake
  )
  add_test(
    NAME phase4.checkpoint_restart
    COMMAND
      ${CMAKE_COMMAND} -DMPS_EXECUTABLE=$<TARGET_FILE:my_planet_sim>
      -DMPS_CONFIG=${CMAKE_CURRENT_SOURCE_DIR}/configs/phase4_manufactured_transport.cfg
      -DMPS_WORKING_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR} -P
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_phase4_checkpoint.cmake
  )
  set_tests_properties(
    phase4.csv_output
    phase4.checkpoint_restart
    PROPERTIES LABELS phase4_gate
  )
  set_tests_properties(
    smoke.phase1_transport
    PROPERTIES PASS_REGULAR_EXPRESSION "result.status = complete"
  )
  set_tests_properties(
    smoke.phase0_ode
    PROPERTIES PASS_REGULAR_EXPRESSION "result.status = complete"
  )
  set_tests_properties(
    smoke.phase2_shallow_water
    PROPERTIES PASS_REGULAR_EXPRESSION "result.status = complete"
  )
  set_tests_properties(
    smoke.phase4_vertical_column
    smoke.phase4_dry_adiabatic
    smoke.phase4_moving_surface_pressure
    smoke.phase4_manufactured_transport
    PROPERTIES PASS_REGULAR_EXPRESSION "result.status = complete"
               LABELS phase4_gate
  )
  # The unit gates drive the driver directly, so a shipped preset can regress without
  # any of them failing. This one carries a tracer through the full linear-reconstruction
  # path and caught exactly that.
  add_test(
    NAME smoke.phase5_solid_body_transport
    COMMAND my_planet_sim --config
            ${CMAKE_CURRENT_SOURCE_DIR}/configs/phase5_solid_body_transport.cfg
  )
  set_tests_properties(
    smoke.phase5_solid_body_transport
    PROPERTIES PASS_REGULAR_EXPRESSION "result.status = complete" LABELS phase5_gate
  )
  add_test(
    NAME phase7.csv_output
    COMMAND
      ${CMAKE_COMMAND} -DMPS_EXECUTABLE=$<TARGET_FILE:my_planet_sim>
      -DMPS_CONFIG=${CMAKE_CURRENT_SOURCE_DIR}/configs/phase7_held_suarez_short.cfg
      -DMPS_WORKING_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR} -P
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_phase7_csv.cmake
  )
  set_tests_properties(phase7.csv_output PROPERTIES LABELS phase7_gate)
  add_test(
    NAME phase8.csv_output
    COMMAND
      ${CMAKE_COMMAND} -DMPS_EXECUTABLE=$<TARGET_FILE:my_planet_sim>
      -DMPS_CONFIG=${CMAKE_CURRENT_SOURCE_DIR}/configs/phase8_earth_geography.cfg
      -DMPS_WORKING_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR} -P
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_phase8_csv.cmake
  )
  set_tests_properties(phase8.csv_output PROPERTIES LABELS phase8_gate)
  add_test(NAME smoke.invalid_cli COMMAND my_planet_sim --unknown value)
  set_tests_properties(smoke.invalid_cli PROPERTIES WILL_FAIL TRUE)
  add_test(
    NAME smoke.invalid_progress_interval
    COMMAND my_planet_sim --config ${CMAKE_CURRENT_SOURCE_DIR}/configs/phase7_held_suarez_short.cfg
            --progress-interval-s invalid
  )
  set_tests_properties(smoke.invalid_progress_interval PROPERTIES WILL_FAIL TRUE)
  add_subdirectory(tests)
endif()
