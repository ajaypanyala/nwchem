
set(IPCCSD_SOURCE_FILES
tce_print_ipx1.F
tce_print_ipx2.F
tce_ipx1_offset.F
ipccsd_x1.F
tce_ipx2_offset.F
tce_jacobi_ipx1.F
tce_jacobi_ipx2.F
ipccsd_x2.F
tce_eom_ipxguess.F
)

add_library(tce_ipccsd STATIC ${IPCCSD_SOURCE_FILES})
#add_dependencies(tce_ipccsd ${NWC_DEPENDENCIES})
target_include_directories(tce_ipccsd PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_ipccsd PUBLIC DEBUG_PRINT)
if(IPCCSD) 
  target_compile_definitions(tce_ipccsd PUBLIC IPCCSD)
endif()

target_compile_options(tce_ipccsd PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_ipccsd PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
