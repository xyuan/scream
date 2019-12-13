# Allow a few variants for specifying the scorpio installation directory
if (SCORPIO_DIR)
  SET (SCORPIO_INSTALL_DIR ${SCORPIO_DIR})
elseif (scorpio_DIR)
  SET (SCORPIO_INSTALL_DIR ${scorpio_DIR})
else ()
  # Build scorpio submodule if user did not specify scorpio_DIR.
  SET(SCORPIO_SRC    ${CMAKE_SOURCE_DIR}/../../externals/scorpio)
  SET(SCORPIO_BINARY ${CMAKE_BINARY_DIR}/scorpio/build)
  SET(PIO_ENABLE_TESTS OFF)
  SET(PIO_TEST_CLOSE_OPEN_FOR_SYNC OFF)

  SET(GENF90_PATH ${CMAKE_SOURCE_DIR}/../../cime/src/externals/genf90)
  SET(PIO_ENABLE_FORTRAN ON)

  SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${SCORPIO_SRC}/cmake")
  add_subdirectory(${SCORPIO_SRC} ${SCORPIO_BINARY})

  SET (SCORPIO_F90_DIR ${SCORPIO_BINARY}/src/flib)

  # ADonahue 12-13-2019:
  # These would only be needed if we wanted to link directly to the C files in scorpio.
  # But there are a number of issues associated with this, so we will stick to using
  # the Fortran modules which accompany SCORPIO.
  SET (SCORPIO_INCLUDE_DIRS ${SCORPIO_BINARY}/src/clib ${SCORPIO_SRC}/src/clib)
  SET (SCORPIO_LIBRARY_DIR ${SCORPIO_SRC}/src/clib)
  SET (SCORPIO_LIB pioc)

endif ()
