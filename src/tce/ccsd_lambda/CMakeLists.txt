
set(CCSD_LAMBDA_SOURCE_FILES
ccsd_lambda1_2.F
ccsd_lambda1_mod.F
ccsd_lambda2.F
ccsd_density1.F
ccsd_lambda1_1.F
)

add_library(tce_ccsd_lambda STATIC ${CCSD_LAMBDA_SOURCE_FILES})
#add_dependencies(tce_ccsd_lambda ${NWC_DEPENDENCIES})
target_include_directories(tce_ccsd_lambda PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_ccsd_lambda PUBLIC DEBUG_PRINT)

target_compile_options(tce_ccsd_lambda PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_ccsd_lambda PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
