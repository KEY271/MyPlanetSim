if(NOT DEFINED MPS_EXECUTABLE OR NOT DEFINED MPS_CONFIG OR
   NOT DEFINED MPS_WORKING_DIRECTORY)
  message(FATAL_ERROR
          "phase4 checkpoint check requires executable, config, and working directory")
endif()

set(test_directory "${MPS_WORKING_DIRECTORY}/phase4_checkpoint_roundtrip")
set(partial_checkpoint "${test_directory}/partial.chk")
set(restarted_checkpoint "${test_directory}/restarted.chk")
set(direct_checkpoint "${test_directory}/direct.chk")
file(MAKE_DIRECTORY "${test_directory}")
file(REMOVE "${partial_checkpoint}" "${restarted_checkpoint}" "${direct_checkpoint}")

execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}" --stop-after-step 20
          --checkpoint "${partial_checkpoint}"
  WORKING_DIRECTORY "${test_directory}"
  RESULT_VARIABLE stopped_result
  OUTPUT_VARIABLE stopped_output
  ERROR_VARIABLE stopped_error
)
if(NOT stopped_result EQUAL 0)
  message(FATAL_ERROR "phase4 stopped run failed: ${stopped_error}")
endif()
if(NOT stopped_output MATCHES "result.status = stopped" OR
   NOT stopped_output MATCHES "result.step = 20")
  message(FATAL_ERROR "phase4 stopped run did not stop at step 20")
endif()
if(NOT EXISTS "${partial_checkpoint}")
  message(FATAL_ERROR "phase4 stopped run did not write a checkpoint file")
endif()
file(READ "${partial_checkpoint}" partial_contents)
if(NOT partial_contents MATCHES "version = 2" OR
   NOT partial_contents MATCHES "layout_id = vertical_column_hybrid_v1")
  message(FATAL_ERROR "phase4 checkpoint has the wrong version or state layout")
endif()

execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}" --restart
          "${partial_checkpoint}" --checkpoint "${restarted_checkpoint}"
  WORKING_DIRECTORY "${test_directory}"
  RESULT_VARIABLE restarted_result
  OUTPUT_VARIABLE restarted_output
  ERROR_VARIABLE restarted_error
)
if(NOT restarted_result EQUAL 0)
  message(FATAL_ERROR "phase4 restarted run failed: ${restarted_error}")
endif()
if(NOT restarted_output MATCHES "result.status = complete")
  message(FATAL_ERROR "phase4 restarted run did not complete")
endif()

execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}" --checkpoint
          "${direct_checkpoint}"
  WORKING_DIRECTORY "${test_directory}"
  RESULT_VARIABLE direct_result
  OUTPUT_VARIABLE direct_output
  ERROR_VARIABLE direct_error
)
if(NOT direct_result EQUAL 0)
  message(FATAL_ERROR "phase4 direct run failed: ${direct_error}")
endif()
if(NOT direct_output MATCHES "result.status = complete")
  message(FATAL_ERROR "phase4 direct run did not complete")
endif()

file(SHA256 "${restarted_checkpoint}" restarted_sha256)
file(SHA256 "${direct_checkpoint}" direct_sha256)
if(NOT restarted_sha256 STREQUAL direct_sha256)
  message(FATAL_ERROR
          "phase4 restarted checkpoint differs from the uninterrupted result")
endif()
