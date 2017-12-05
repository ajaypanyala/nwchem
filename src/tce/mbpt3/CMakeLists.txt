
set(MBPT3_SOURCE_FILES
mbpt3_t1.F
mbpt3_t1_sd.F
mbpt3_e.F
mbpt3_t2.F
mbpt3_t2_sd.F
)

add_library(tce_mbpt3 STATIC ${MBPT3_SOURCE_FILES})
#add_dependencies(tce_mbpt3 ${NWC_DEPENDENCIES})
target_include_directories(tce_mbpt3 PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_mbpt3 PUBLIC DEBUG_PRINT)

target_compile_options(tce_mbpt3 PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_mbpt3 PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
