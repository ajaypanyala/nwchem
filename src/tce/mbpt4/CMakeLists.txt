
set(MBPT4_SOURCE_FILES
mbpt4_e.F
)

add_library(tce_mbpt4 STATIC ${MBPT4_SOURCE_FILES})
#add_dependencies(tce_mbpt4 ${NWC_DEPENDENCIES})
target_include_directories(tce_mbpt4 PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_mbpt4 PUBLIC DEBUG_PRINT)

target_compile_options(tce_mbpt4 PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_mbpt4 PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
