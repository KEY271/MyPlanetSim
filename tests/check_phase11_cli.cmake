function(check_radiation_restart)
set(continuous_checkpoint "${MPS_WORKING_DIRECTORY}/phase11-continuous.chk")
set(partial_checkpoint "${MPS_WORKING_DIRECTORY}/phase11-partial.chk")
set(restarted_checkpoint "${MPS_WORKING_DIRECTORY}/phase11-restarted.chk")

execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}" --progress-interval-s 0
          --checkpoint "${continuous_checkpoint}"
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE stdout
  ERROR_VARIABLE stderr
)
if(NOT result EQUAL 0 OR NOT stdout MATCHES "result.status = complete")
  message(FATAL_ERROR "Phase 11 CLI run failed: ${stderr}\n${stdout}")
endif()

set(output_directory "${MPS_WORKING_DIRECTORY}/phase11-diagnostics-output")
file(READ "${output_directory}/radiation_diagnostics.csv" diagnostics)
file(STRINGS "${output_directory}/radiation_diagnostics.csv" continuous_rows)
file(READ "${output_directory}/radiation_column.csv" column)
file(COPY_FILE "${output_directory}/radiation_column.csv"
     "${MPS_WORKING_DIRECTORY}/phase11-continuous-column.csv")
file(READ "${output_directory}/surface_state.csv" surface)
if(NOT diagnostics MATCHES "toa_incoming_sw_w_m2,toa_reflected_sw_w_m2,olr_w_m2")
  message(FATAL_ERROR "radiation diagnostics rate columns are missing")
endif()
if(NOT diagnostics MATCHES "segment_start_time_s,segment_start_step")
  message(FATAL_ERROR "radiation segment-combination metadata is missing")
endif()
if(NOT diagnostics MATCHES
   "interval_surface_storage_change_j,interval_surface_time_integration_residual_j")
  message(FATAL_ERROR "radiation diagnostics energy columns are missing")
endif()
if(NOT diagnostics MATCHES
   "cfl_retry_count,invariant_retry_count,solver_retry_count")
  message(FATAL_ERROR "radiation retry counters are missing")
endif()
if(NOT column MATCHES "record_kind,index,pressure_pa,shortwave_optical_depth" OR
   NOT column MATCHES "interface" OR NOT column MATCHES "full")
  message(FATAL_ERROR "radiation column schema is incomplete")
endif()
if(NOT surface MATCHES "surface_temperature_k")
  message(FATAL_ERROR "gray radiation surface state is missing")
endif()
if(diagnostics MATCHES "nan|inf" OR column MATCHES "nan|inf" OR
   surface MATCHES "nan|inf")
  message(FATAL_ERROR "Phase 11 CSV contains a non-finite value")
endif()

execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}" --progress-interval-s 0
          --stop-after-step 1 --checkpoint "${partial_checkpoint}"
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE partial_result
  OUTPUT_QUIET
  ERROR_VARIABLE partial_stderr
)
if(NOT partial_result EQUAL 0)
  message(FATAL_ERROR "Phase 11 partial run failed: ${partial_stderr}")
endif()
file(STRINGS "${output_directory}/radiation_diagnostics.csv" partial_rows)
execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}" --progress-interval-s 0
          --restart "${partial_checkpoint}" --checkpoint "${restarted_checkpoint}"
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE restart_result
  OUTPUT_QUIET
  ERROR_VARIABLE restart_stderr
)
if(NOT restart_result EQUAL 0)
  message(FATAL_ERROR "Phase 11 restart failed: ${restart_stderr}")
endif()
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E compare_files "${continuous_checkpoint}"
          "${restarted_checkpoint}"
  RESULT_VARIABLE checkpoint_difference
)
if(NOT checkpoint_difference EQUAL 0)
  message(FATAL_ERROR "Phase 11 restarted checkpoint differs from continuous run")
endif()
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E compare_files
          "${MPS_WORKING_DIRECTORY}/phase11-continuous-column.csv"
          "${output_directory}/radiation_column.csv"
  RESULT_VARIABLE column_difference
)
if(NOT column_difference EQUAL 0)
  message(FATAL_ERROR "Phase 11 restarted orbit or flux differs from continuous run")
endif()

# Compare each accepted interval, including all energy fields and retry counts.
# Restart's initial row is instantaneous and has no integrated contribution.
file(STRINGS "${output_directory}/radiation_diagnostics.csv" restarted_rows)
list(REMOVE_AT partial_rows 0 1)
list(REMOVE_AT restarted_rows 0 1)
list(APPEND partial_rows ${restarted_rows})
list(REMOVE_AT continuous_rows 0 1)
list(LENGTH continuous_rows expected_count)
list(LENGTH partial_rows actual_count)
if(NOT expected_count EQUAL actual_count)
  message(FATAL_ERROR "restart changed the accepted interval count")
endif()
math(EXPR last "${expected_count} - 1")
foreach(index RANGE 0 ${last})
  list(GET continuous_rows ${index} expected)
  list(GET partial_rows ${index} actual)
  string(REPLACE "," ";" expected "${expected}")
  string(REPLACE "," ";" actual "${actual}")
  # Segment metadata, cached RHS call counts and wall times legitimately differ.
  list(SUBLIST expected 4 35 expected_values)
  list(SUBLIST actual 4 35 actual_values)
  if(NOT expected_values STREQUAL actual_values)
    message(FATAL_ERROR "restart changed flux, integrated energy or counters at ${index}")
  endif()
endforeach()

endfunction()

file(READ "${MPS_CONFIG}" base_config)
# Keep a nonzero orbit origin across checkpoint dispatch.
string(REPLACE "run.start_time_s = 0" "run.start_time_s = 12345" base_config "${base_config}")
string(REPLACE "run.end_time_s = 20" "run.end_time_s = 12385" base_config "${base_config}")
foreach(integrator IN ITEMS ssprk3 semi_implicit)
  set(MPS_CONFIG "${MPS_WORKING_DIRECTORY}/phase11-${integrator}.cfg")
  set(case_config "${base_config}")
  if(integrator STREQUAL "semi_implicit")
    string(APPEND case_config "
dry_hydrostatic.time_integrator = semi_implicit
dry_hydrostatic.advective_cfl = 0.45
semi_implicit.reference_surface_pressure_pa = 100000
semi_implicit.reference_temperature_k = 280
semi_implicit.reference_update = fixed
semi_implicit.implicit_weight = 0.5
semi_implicit.wave_cfl_threshold = 0.45
semi_implicit.maximum_implicit_modes = 2
semi_implicit.nonlinear_iterations = 5
semi_implicit.linear_relative_tolerance = 1e-8
semi_implicit.linear_absolute_tolerance = 1e-12
semi_implicit.linear_maximum_iterations = 80
semi_implicit.gmres_restart = 20
semi_implicit.minimum_time_step_s = 0.01
")
  endif()
  file(WRITE "${MPS_CONFIG}" "${case_config}")
  check_radiation_restart()
endforeach()
