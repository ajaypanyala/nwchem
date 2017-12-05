
set(LCCSD_SOURCE_FILES
lccsd_t2.F
lccsd_e.F
lccsd_t1.F
)

add_library(tce_lccsd STATIC ${LCCSD_SOURCE_FILES})
#add_dependencies(tce_lccsd ${NWC_DEPENDENCIES})
target_include_directories(tce_lccsd PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_lccsd PUBLIC DEBUG_PRINT)
target_compile_options(tce_lccsd PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_lccsd PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
