
set(MBPT2_SOURCE_FILES
mbpt2_t1.F
mbpt2_t2.F
mbpt2_t1_sd.F
mbpt2_t2_b.F
mbpt2_e.F
mbpt2_t4.F
mbpt2_t2_sd.F
mbpt2_t3.F
)

add_library(tce_mbpt2 STATIC ${MBPT2_SOURCE_FILES})
#add_dependencies(tce_mbpt2 ${NWC_DEPENDENCIES})
target_include_directories(tce_mbpt2 PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_mbpt2 PUBLIC DEBUG_PRINT)

target_compile_options(tce_mbpt2 PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_mbpt2 PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
