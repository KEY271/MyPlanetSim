execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}" --progress-interval-s 0
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE stdout
  ERROR_VARIABLE stderr
)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Phase 10 diagnostics run failed: ${stderr}")
endif()
if(NOT stdout MATCHES "result.status = complete" OR
   NOT stdout MATCHES "semi_implicit.vertical_mode_count = 4")
  message(FATAL_ERROR "Phase 10 metadata is incomplete: ${stdout}")
endif()

set(output_directory "${MPS_WORKING_DIRECTORY}/phase10-diagnostics-output")
file(READ "${output_directory}/semi_implicit_diagnostics.csv" diagnostics)
if(NOT diagnostics MATCHES
   "requested_dt_s,accepted_dt_s,advective_cfl,implicit_wave_cfl,vertical_cfl")
  message(FATAL_ERROR "semi-implicit diagnostics header is incomplete")
endif()
if(NOT diagnostics MATCHES "linear_iterations_total,linear_iterations_maximum")
  message(FATAL_ERROR "semi-implicit solver diagnostics header is incomplete")
endif()
if(diagnostics MATCHES "nan|inf")
  message(FATAL_ERROR "semi-implicit diagnostics contain a non-finite value")
endif()
string(REGEX MATCHALL "\n" diagnostic_newlines "${diagnostics}")
list(LENGTH diagnostic_newlines diagnostic_line_count)
if(NOT diagnostic_line_count EQUAL 3)
  message(FATAL_ERROR "semi-implicit CSV must contain one header and two samples")
endif()
