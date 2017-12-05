
set(CCSD_ACT_SOURCE_FILES
ccsd_t2a.F
diag2_act2.F
ccsd_t1a.F
activity.F
ccsd_ea.F
diag1_act2.F
)

add_library(tce_ccsd_act STATIC ${CCSD_ACT_SOURCE_FILES})
#add_dependencies(tce_ccsd_act ${NWC_DEPENDENCIES})
target_include_directories(tce_ccsd_act PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_ccsd_act PUBLIC DEBUG_PRINT)

target_compile_options(tce_ccsd_act PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_ccsd_act PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
