
set(CCD_SOURCE_FILES
ccd_e.F
ccd_t2.F
)

add_library(tce_ccd STATIC ${CCD_SOURCE_FILES})
#add_dependencies(tce_ccd ${NWC_DEPENDENCIES})
target_include_directories(tce_ccd PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_ccd PUBLIC DEBUG_PRINT)


target_compile_options(tce_ccd PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_ccd PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
