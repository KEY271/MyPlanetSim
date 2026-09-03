function(mps_add_test target source)
  add_executable(${target} ${source})
  target_compile_features(${target} PRIVATE cxx_std_20)
  target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(
    ${target}
    PRIVATE myplanetsim_core myplanetsim_warnings myplanetsim_sanitizers
  )
  add_test(NAME ${target} COMMAND ${target})
endfunction()

mps_add_test(unit.test_support test_support.cpp)
mps_add_test(unit.test_support_failure test_support_failure.cpp)
set_tests_properties(unit.test_support_failure PROPERTIES WILL_FAIL TRUE)
mps_add_test(unit.validation test_validation.cpp)
mps_add_test(unit.planet_parameters test_planet_parameters.cpp)
mps_add_test(unit.experiment_config test_experiment_config.cpp)
mps_add_test(unit.field2d test_field2d.cpp)
mps_add_test(unit.periodic_boundary test_periodic_boundary.cpp)
