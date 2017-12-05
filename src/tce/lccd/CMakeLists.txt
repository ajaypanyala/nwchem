
set(LCCD_SOURCE_FILES
lccd_t2.F
lccd_e.F
)

add_library(tce_lccd STATIC ${LCCD_SOURCE_FILES})
#add_dependencies(tce_lccd ${NWC_DEPENDENCIES})
target_include_directories(tce_lccd PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_lccd PUBLIC DEBUG_PRINT)
target_compile_options(tce_lccd PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_lccd PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
