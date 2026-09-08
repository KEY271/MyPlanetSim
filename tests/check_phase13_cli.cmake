file(READ "${MPS_CONFIG}" case_config)
string(REPLACE "output.directory = output/phase13_moist_ci"
               "output.directory = ${MPS_WORKING_DIRECTORY}/phase13-moist-output"
               case_config "${case_config}")
set(case_path "${MPS_WORKING_DIRECTORY}/phase13-moist-ci.cfg")
set(continuous_checkpoint
    "${MPS_WORKING_DIRECTORY}/phase13-moist-continuous.chk")
set(partial_checkpoint "${MPS_WORKING_DIRECTORY}/phase13-moist-partial.chk")
set(restarted_checkpoint "${MPS_WORKING_DIRECTORY}/phase13-moist-restarted.chk")
file(WRITE "${case_path}" "${case_config}")

execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${case_path}" --progress-interval-s 0
          --checkpoint "${continuous_checkpoint}"
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE result
  ERROR_VARIABLE stderr
)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Phase 13 moist CLI failed: ${stderr}")
endif()

set(output_directory "${MPS_WORKING_DIRECTORY}/phase13-moist-output")
foreach(file IN ITEMS tracer_diagnostics.csv moisture_diagnostics.csv
                      moist_convection_diagnostics.csv moist_column.csv
                      moist_surface_state.csv)
  if(NOT EXISTS "${output_directory}/${file}")
    message(FATAL_ERROR "Phase 13 output is missing ${file}")
  endif()
  file(READ "${output_directory}/${file}" contents)
  if(contents MATCHES "(^|,)(nan|inf)(,|$)")
    message(FATAL_ERROR "Phase 13 ${file} contains non-finite data")
  endif()
endforeach()
file(READ "${output_directory}/tracer_diagnostics.csv" tracer)
file(READ "${output_directory}/moisture_diagnostics.csv" moisture)
file(READ "${output_directory}/moist_convection_diagnostics.csv" convection)
file(READ "${output_directory}/moist_column.csv" column)
if(NOT tracer MATCHES "interval_transport_diffusion_residual_kg")
  message(FATAL_ERROR "Phase 13 tracer diagnostics schema is incomplete")
endif()
if(NOT moisture MATCHES
   "precipitable_water_kg_m2,area_mean_relative_humidity" OR
   NOT moisture MATCHES "dry_latent_surface_energy_j" OR
   NOT moisture MATCHES "water_budget_residual_kg")
  message(FATAL_ERROR "Phase 13 moisture diagnostics schema is incomplete")
endif()
if(NOT convection MATCHES "area_mean_cape_j_kg" OR
   NOT convection MATCHES "deep_area_fraction")
  message(FATAL_ERROR "Phase 13 convection diagnostics schema is incomplete")
endif()
if(NOT column MATCHES
   "vapor_mixing_ratio,saturation_mixing_ratio,relative_humidity")
  message(FATAL_ERROR "Phase 13 moist column schema is incomplete")
endif()

execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${case_path}" --progress-interval-s 0
          --stop-after-step 1 --checkpoint "${partial_checkpoint}"
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE partial_result
  ERROR_VARIABLE partial_stderr
)
if(NOT partial_result EQUAL 0)
  message(FATAL_ERROR "Phase 13 partial run failed: ${partial_stderr}")
endif()
execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${case_path}" --progress-interval-s 0
          --restart "${partial_checkpoint}" --checkpoint "${restarted_checkpoint}"
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE restart_result
  ERROR_VARIABLE restart_stderr
)
if(NOT restart_result EQUAL 0)
  message(FATAL_ERROR "Phase 13 restart failed: ${restart_stderr}")
endif()
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E compare_files "${continuous_checkpoint}"
          "${restarted_checkpoint}"
  RESULT_VARIABLE checkpoint_difference
)
if(NOT checkpoint_difference EQUAL 0)
  message(FATAL_ERROR "Phase 13 restart changed the final checkpoint")
endif()
