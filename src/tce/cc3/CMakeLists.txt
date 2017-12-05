
set(CC3_SOURCE_FILES
eomccsdt3_x3.F
ccsdt2_t3.F
ccsdt1b_t3.F
eomccsdt1_x3.F
cc3_lambda3.F
ccsdt3_t3.F
cc3_t3.F
eomcc3_x3.F
)

add_library(tce_cc3 STATIC ${CC3_SOURCE_FILES})
#add_dependencies(tce_cc3 ${NWC_DEPENDENCIES})
target_include_directories(tce_cc3 PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_cc3 PUBLIC DEBUG_PRINT)

target_compile_options(tce_cc3 PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_cc3 PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
