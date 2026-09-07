function(check_dry_mixing_cli integrator)
  file(READ "${MPS_CONFIG}" case_config)
  string(REPLACE "surface.air_exchange_coefficient_w_m2_k = 8"
                 "surface.air_exchange_coefficient_w_m2_k = 0"
                 case_config "${case_config}")
  string(REPLACE "output.directory = phase11-diagnostics-output"
                 "output.directory = phase12-diagnostics-output"
                 case_config "${case_config}")
  string(APPEND case_config "
convection.kind = dry_adjustment
convection.stability_tolerance_k = 1e-10
boundary_layer.kind = bulk_k_profile
boundary_layer.integrator = backward_euler
boundary_layer.critical_richardson = 1
boundary_layer.turbulent_prandtl = 1
boundary_layer.gustiness_m_s = 1
surface.land_roughness_momentum_m = 0.1
surface.land_roughness_heat_m = 0.01
surface.ocean_roughness_momentum_m = 0.001
surface.ocean_roughness_heat_m = 0.0001
")
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
  set(case_path "${MPS_WORKING_DIRECTORY}/phase12-${integrator}.cfg")
  set(continuous_checkpoint
      "${MPS_WORKING_DIRECTORY}/phase12-${integrator}-continuous.chk")
  set(partial_checkpoint
      "${MPS_WORKING_DIRECTORY}/phase12-${integrator}-partial.chk")
  set(restarted_checkpoint
      "${MPS_WORKING_DIRECTORY}/phase12-${integrator}-restarted.chk")
  file(WRITE "${case_path}" "${case_config}")
  execute_process(
    COMMAND "${MPS_EXECUTABLE}" --config "${case_path}" --progress-interval-s 0
            --checkpoint "${continuous_checkpoint}"
    WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
    RESULT_VARIABLE result
    ERROR_VARIABLE stderr
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Phase 12 ${integrator} CLI failed: ${stderr}")
  endif()
  set(output_directory "${MPS_WORKING_DIRECTORY}/phase12-diagnostics-output")
  file(READ "${output_directory}/convection_diagnostics.csv" convection)
  file(READ "${output_directory}/boundary_layer_diagnostics.csv" boundary_layer)
  file(READ "${output_directory}/mixing_column.csv" mixing_column)
  if(NOT convection MATCHES
     "minimum_theta_difference_before_k,minimum_theta_difference_after_k")
    message(FATAL_ERROR "Phase 12 convection diagnostics schema is incomplete")
  endif()
  if(NOT boundary_layer MATCHES
     "physical_shear_dissipation_j,physical_surface_drag_dissipation_j")
    message(FATAL_ERROR "Phase 12 boundary-layer diagnostics schema is incomplete")
  endif()
  if(NOT mixing_column MATCHES "record_kind,index,pressure_pa,height_m" OR
     NOT mixing_column MATCHES "final_diagnosed_surface" OR
     NOT mixing_column MATCHES "state_after_mixing")
    message(FATAL_ERROR "Phase 12 mixing column schema is incomplete")
  endif()
  if(convection MATCHES "nan|inf" OR boundary_layer MATCHES "nan|inf" OR
     mixing_column MATCHES "nan|inf")
    message(FATAL_ERROR "Phase 12 mixing CSV contains non-finite data")
  endif()

  execute_process(
    COMMAND "${MPS_EXECUTABLE}" --config "${case_path}" --progress-interval-s 0
            --stop-after-step 1 --checkpoint "${partial_checkpoint}"
    WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
    RESULT_VARIABLE partial_result
    ERROR_VARIABLE partial_stderr
  )
  if(NOT partial_result EQUAL 0)
    message(FATAL_ERROR "Phase 12 partial run failed: ${partial_stderr}")
  endif()
  execute_process(
    COMMAND "${MPS_EXECUTABLE}" --config "${case_path}" --progress-interval-s 0
            --restart "${partial_checkpoint}" --checkpoint "${restarted_checkpoint}"
    WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
    RESULT_VARIABLE restart_result
    ERROR_VARIABLE restart_stderr
  )
  if(NOT restart_result EQUAL 0)
    message(FATAL_ERROR "Phase 12 restart failed: ${restart_stderr}")
  endif()
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E compare_files "${continuous_checkpoint}"
            "${restarted_checkpoint}"
    RESULT_VARIABLE checkpoint_difference
  )
  if(NOT checkpoint_difference EQUAL 0)
    message(FATAL_ERROR "Phase 12 restart changed the final checkpoint")
  endif()
endfunction()

check_dry_mixing_cli(ssprk3)
check_dry_mixing_cli(semi_implicit)
