# NOTE: if you have baseline-type tests, add the subdirectory OUTSIDE the following if statement
if (NOT ${SCREAM_BASELINES_ONLY})
  if ("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
    add_subdirectory(homme)
  endif()

  add_subdirectory(p3)
  add_subdirectory(shoc)
  add_subdirectory(rrtmgp)
  if (SCREAM_ENABLE_EXPERIMENTAL)
    add_subdirectory(zm)
  endif()
endif()
