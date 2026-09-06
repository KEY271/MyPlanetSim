execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}"
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE stdout
  ERROR_VARIABLE stderr
)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Phase 7 short run failed: ${stderr}")
endif()
if(NOT stdout MATCHES "result.status = complete")
  message(FATAL_ERROR "Phase 7 short run did not report completion")
endif()
if(NOT stderr MATCHES "progress status=running" OR
   NOT stderr MATCHES "progress status=complete")
  message(FATAL_ERROR "Phase 7 short run did not report progress: ${stderr}")
endif()

set(output_directory "${MPS_WORKING_DIRECTORY}/phase7-short-output")
file(READ "${output_directory}/physics_diagnostics.csv" physics)
file(READ "${output_directory}/climate_statistics.csv" climate)
if(NOT physics MATCHES "thermal_energy_rate_w,rayleigh_drag_work_w")
  message(FATAL_ERROR "physics diagnostics header is incomplete")
endif()
if(NOT climate MATCHES "latitude_lower_deg,latitude_upper_deg,level")
  message(FATAL_ERROR "climate statistics header is incomplete")
endif()
if(physics MATCHES "nan|inf" OR climate MATCHES "nan|inf")
  message(FATAL_ERROR "Phase 7 CSV contains a non-finite value")
endif()
string(REGEX MATCHALL "\n" climate_newlines "${climate}")
list(LENGTH climate_newlines climate_line_count)
if(NOT climate_line_count EQUAL 9)
  message(FATAL_ERROR "climate CSV must contain one header and 8 rows")
endif()

execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}" --stop-after-step 1
          --progress-interval-s 0
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE stopped_result
  OUTPUT_VARIABLE stopped_stdout
  ERROR_VARIABLE stopped_stderr
)
if(NOT stopped_result EQUAL 0 OR
   NOT stopped_stdout MATCHES "result.status = stopped" OR
   NOT stopped_stdout MATCHES "result.step = 1")
  message(FATAL_ERROR "Phase 7 controlled stop failed: ${stopped_stdout}${stopped_stderr}")
endif()
if(stopped_stderr MATCHES "progress status=")
  message(FATAL_ERROR "--progress-interval-s 0 did not disable progress")
endif()
