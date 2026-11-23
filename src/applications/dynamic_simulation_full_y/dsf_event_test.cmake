# -------------------------------------------------------------
# Run dsf.x and dsf2.x on the same input and compare the watch files.
# Arguments: DSF, DSF2, DIR (folder with the inputs), INPUT, RAW, DYR, WATCH
# -------------------------------------------------------------
foreach(tag dsf dsf2)
  set(work "${DIR}/event_test/${tag}")
  file(REMOVE_RECURSE "${work}")
  file(MAKE_DIRECTORY "${work}")
  foreach(f ${INPUT} ${RAW} ${DYR})
    file(COPY "${DIR}/${f}" DESTINATION "${work}")
  endforeach()
endforeach()

execute_process(COMMAND "${DSF}" ${INPUT}
  WORKING_DIRECTORY "${DIR}/event_test/dsf"
  RESULT_VARIABLE status OUTPUT_QUIET)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "dsf.x failed: ${status}")
endif()

execute_process(COMMAND "${DSF2}" ${INPUT}
  WORKING_DIRECTORY "${DIR}/event_test/dsf2"
  RESULT_VARIABLE status OUTPUT_QUIET)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "dsf2.x failed: ${status}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files
  "${DIR}/event_test/dsf/${WATCH}" "${DIR}/event_test/dsf2/${WATCH}"
  RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "dsf.x and dsf2.x watch files differ")
endif()
