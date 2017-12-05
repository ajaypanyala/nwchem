
set(EACCSD_SOURCE_FILES
tce_jacobi_eax2.F
eaccsd_x2.F
tce_print_eax2.F
tce_jacobi_eax1.F
tce_print_eax1.F
eaccsd_x1.F
tce_eom_eaxguess.F
tce_eax2_offset.F
tce_eax1_offset.F
)

add_library(tce_eaccsd STATIC ${EACCSD_SOURCE_FILES})
#add_dependencies(tce_eaccsd ${NWC_DEPENDENCIES})
target_include_directories(tce_eaccsd PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_eaccsd PUBLIC DEBUG_PRINT)
if(EACCSD) 
  target_compile_definitions(tce_eaccsd PUBLIC EACCSD)
endif()

target_compile_options(tce_eaccsd PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_eaccsd PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
