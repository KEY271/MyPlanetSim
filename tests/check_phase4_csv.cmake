if(NOT DEFINED MPS_EXECUTABLE OR NOT DEFINED MPS_CONFIG OR
   NOT DEFINED MPS_WORKING_DIRECTORY)
  message(FATAL_ERROR "phase4 CSV check requires executable, config, and working directory")
endif()

execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}"
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE run_result
  OUTPUT_VARIABLE run_output
  ERROR_VARIABLE run_error
)
if(NOT run_result EQUAL 0)
  message(FATAL_ERROR "phase4 CSV run failed: ${run_error}")
endif()
if(NOT run_output MATCHES "result.status = complete")
  message(FATAL_ERROR "phase4 CSV run did not complete")
endif()

set(output_directory
    "${MPS_WORKING_DIRECTORY}/output/phase4_manufactured_transport")
file(READ "${output_directory}/column_profile.csv" profile)
file(READ "${output_directory}/column_diagnostics.csv" diagnostics)
if(NOT profile MATCHES "location,index,pressure_pa,geopotential_m2_s2,height_m")
  message(FATAL_ERROR "column profile header is missing location or physical units")
endif()
if(NOT profile MATCHES "interface,0," OR NOT profile MATCHES "full,0,")
  message(FATAL_ERROR "column profile does not distinguish interface and full levels")
endif()
if(NOT diagnostics MATCHES "dry_mass_budget_residual_kg_m2" OR
   NOT diagnostics MATCHES "continuity_residual_pa_s" OR
   NOT diagnostics MATCHES "hydrostatic_l2_residual_m2_s2" OR
   NOT diagnostics MATCHES "non_finite_count")
  message(FATAL_ERROR "column diagnostics CSV is missing Phase 4 budget fields")
endif()
