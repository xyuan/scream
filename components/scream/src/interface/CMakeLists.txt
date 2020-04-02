set(ATM_SRC
  scream_cxx_f90_interface.cpp
  atm_comp_mct.F90
  dead_mct_mod.F90
  dead_data_mod.F90
  dead_mod.F90
  spmd_utils.F90
  mpishorthand.F
  units.F90
  examplePio.f90 
)

if (USE_ESMF_LIB)
  set(ESMFDIR "esmf")
else()
  set(ESMFDIR "noesmf")
endif()

set(CSM_SHARE "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/csm_share")
#set(PIO_LIBDIR "${INSTALL_SHAREDPATH}/lib")
#if (PIO_VERSION STREQUAL 2)
#  # This is a pio2 library
#  set(PIOLIB "${PIO_LIBDIR}/libpiof.a ${PIO_LIBDIR}/libpioc.a")
#  set(PIOLIBNAME "-lpiof -lpioc")
#  set(PIO_SRC_DIR "${CIMEROOT}/src/externals/scorpio")
#else()
#  # This is a pio1 library
#  set(PIOLIB "${PIO_LIBDIR}/libpio.a")
#  set(PIOLIBNAME "-lpio")
#  if (NOT EXISTS "${CIMEROOT}/src/externals/scorpio_classic/pio")
#    set(PIO_SRC_DIR "${CIMEROOT}/src/externals/scorpio_classic")
#  else()
#    set(PIO_SRC_DIR "${CIMEROOT}/src/externals/scorpio_classic/pio")
#  endif()
#endif()
message("ASD - PIO_SRC_DIR: ${PIO_SRC_DIR}")
message("ASD - PIOLIB: ${PIOLIB}")
message("ASD - PIO_MODDIR: ${PIO_MODDIR}")
link_directories(${SCREAM_TPL_LIBRARY_DIRS})

add_library(atm ${ATM_SRC})
target_include_directories(atm PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(atm SYSTEM PUBLIC ${SCREAM_TPL_INCLUDE_DIRS} ${CSM_SHARE} ${INSTALL_SHAREDPATH}/include ${PIO_MODDIR} ${INCLUDE_NETCDF_FORTRAN})
set_target_properties(atm PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/cmake/atm)
target_link_libraries(atm scream_share scream_control p3 shoc ${SCREAM_TPL_LIBRARIES} ${PIOLIB})
