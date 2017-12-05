
set(CISD_SOURCE_FILES
cisd_c1.F
cisd_e.F
cisd_c2.F
)

add_library(tce_cisd STATIC ${CISD_SOURCE_FILES})
#add_dependencies(tce_cisd ${NWC_DEPENDENCIES})
target_include_directories(tce_cisd PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_cisd PUBLIC DEBUG_PRINT)
target_compile_options(tce_cisd PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_cisd PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
