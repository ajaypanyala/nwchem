
set(MBPT1_SOURCE_FILES
mbpt1_t2.F
mbpt1_t1.F
)

add_library(tce_mbpt1 STATIC ${MBPT1_SOURCE_FILES})
#add_dependencies(tce_mbpt1 ${NWC_DEPENDENCIES})
target_include_directories(tce_mbpt1 PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_mbpt1 PUBLIC DEBUG_PRINT)

target_compile_options(tce_mbpt1 PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_mbpt1 PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
