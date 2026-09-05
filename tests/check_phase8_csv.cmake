execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}"
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE stdout
  ERROR_VARIABLE stderr
)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Phase 8 surface run failed: ${stderr}")
endif()
if(NOT stdout MATCHES "result.status = complete")
  message(FATAL_ERROR "Phase 8 surface run did not report completion")
endif()

set(output_directory "${MPS_WORKING_DIRECTORY}/output/phase8_earth_geography")
file(READ "${output_directory}/surface_state.csv" surface_state)
file(READ "${output_directory}/surface_diagnostics.csv" surface_diagnostics)
if(NOT surface_state MATCHES "panel,i,j,longitude_deg,latitude_deg,height_m,land_fraction")
  message(FATAL_ERROR "surface state header is incomplete")
endif()
if(NOT surface_diagnostics MATCHES
   "surface_storage_rate_w,interval_absorbed_stellar_energy_j")
  message(FATAL_ERROR "surface diagnostics header is incomplete")
endif()
if(NOT surface_diagnostics MATCHES "interval_surface_budget_residual_j")
  message(FATAL_ERROR "surface interval budget is missing")
endif()
if(surface_state MATCHES "nan|inf" OR surface_diagnostics MATCHES "nan|inf")
  message(FATAL_ERROR "Phase 8 CSV contains a non-finite value")
endif()
string(REGEX MATCHALL "\n" state_newlines "${surface_state}")
list(LENGTH state_newlines state_line_count)
if(NOT state_line_count EQUAL 97)
  message(FATAL_ERROR "surface state CSV must contain one header and 96 cells")
endif()
