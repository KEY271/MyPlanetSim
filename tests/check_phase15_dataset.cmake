set(dataset_root "${MPS_WORKING_DIRECTORY}/output/phase15-viewer-rest-n4/viewer")
set(continuous_checkpoint "${MPS_WORKING_DIRECTORY}/phase15-continuous.chk")
set(split_checkpoint "${MPS_WORKING_DIRECTORY}/phase15-split.chk")
set(split_final_checkpoint "${MPS_WORKING_DIRECTORY}/phase15-split-final.chk")
file(REMOVE_RECURSE "${MPS_WORKING_DIRECTORY}/output/phase15-viewer-rest-n4")
file(REMOVE "${continuous_checkpoint}" "${continuous_checkpoint}.statistics"
  "${split_checkpoint}" "${split_checkpoint}.statistics"
  "${split_final_checkpoint}" "${split_final_checkpoint}.statistics")

execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}" --progress-interval-s 0
    --checkpoint "${continuous_checkpoint}"
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE error
)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Phase 15 fixture failed (${result}):\n${output}\n${error}")
endif()

set(manifest_path "${dataset_root}/manifest.json")
set(terrain_path "${dataset_root}/terrain.bin")
foreach(required IN ITEMS "${manifest_path}" "${terrain_path}")
  if(NOT EXISTS "${required}")
    message(FATAL_ERROR "Phase 15 dataset is missing ${required}")
  endif()
endforeach()

file(SHA256 "${continuous_checkpoint}" continuous_checkpoint_hash)
set(continuous_period_hashes)
foreach(index RANGE 0 1)
  file(SHA256 "${dataset_root}/means/period_00000${index}.bin" period_hash)
  list(APPEND continuous_period_hashes "${period_hash}")
endforeach()

file(REMOVE_RECURSE "${MPS_WORKING_DIRECTORY}/output/phase15-viewer-rest-n4")
execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}" --progress-interval-s 0
    --stop-after-step 6 --checkpoint "${split_checkpoint}"
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE split_result
  OUTPUT_VARIABLE split_output
  ERROR_VARIABLE split_error
)
if(NOT split_result EQUAL 0)
  message(FATAL_ERROR "Phase 15 split fixture failed (${split_result}):\n${split_output}\n${split_error}")
endif()
execute_process(
  COMMAND "${MPS_EXECUTABLE}" --config "${MPS_CONFIG}" --progress-interval-s 0
    --restart "${split_checkpoint}" --checkpoint "${split_final_checkpoint}"
  WORKING_DIRECTORY "${MPS_WORKING_DIRECTORY}"
  RESULT_VARIABLE restart_result
  OUTPUT_VARIABLE restart_output
  ERROR_VARIABLE restart_error
)
if(NOT restart_result EQUAL 0)
  message(FATAL_ERROR "Phase 15 restart fixture failed (${restart_result}):\n${restart_output}\n${restart_error}")
endif()

file(SHA256 "${split_final_checkpoint}" split_checkpoint_hash)
if(NOT split_checkpoint_hash STREQUAL continuous_checkpoint_hash)
  message(FATAL_ERROR "Phase 15 restart changed the final model checkpoint")
endif()
foreach(index RANGE 0 1)
  list(GET continuous_period_hashes ${index} continuous_period_hash)
  file(SHA256 "${dataset_root}/means/period_00000${index}.bin" split_period_hash)
  if(NOT split_period_hash STREQUAL continuous_period_hash)
    message(FATAL_ERROR "Phase 15 restart changed period ${index}")
  endif()
endforeach()

file(READ "${manifest_path}" manifest)
string(JSON schema_version GET "${manifest}" schemaVersion)
string(JSON period_count LENGTH "${manifest}" periods)
if(NOT schema_version EQUAL 1 OR NOT period_count EQUAL 2)
  message(FATAL_ERROR "unexpected Phase 15 manifest schema or period count")
endif()

foreach(index RANGE 0 1)
  string(JSON complete GET "${manifest}" periods ${index} complete)
  string(JSON relative_path GET "${manifest}" periods ${index} path)
  string(JSON declared_size GET "${manifest}" periods ${index} byteLength)
  set(period_path "${dataset_root}/${relative_path}")
  if(NOT complete OR NOT EXISTS "${period_path}")
    message(FATAL_ERROR "Phase 15 period ${index} is incomplete or missing")
  endif()
  file(SIZE "${period_path}" actual_size)
  if(NOT actual_size EQUAL declared_size)
    message(FATAL_ERROR "Phase 15 period ${index} byte length mismatch")
  endif()
endforeach()
