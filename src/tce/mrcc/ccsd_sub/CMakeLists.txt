
set(CCSD_SUB_SOURCE_FILES
ccsd_dcl_sub.F
ccsd_st2.F
ccsd_st1.F
ccsd_se.F
)

add_library(tce_mrcc_ccsd_sub STATIC ${CCSD_SUB_SOURCE_FILES})
#add_dependencies(tce_mrcc_ccsd_sub ${NWC_DEPENDENCIES})
target_include_directories(tce_mrcc_ccsd_sub PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})

target_compile_definitions(tce_mrcc_ccsd_sub PUBLIC MRCC_METHODS)
target_compile_options(tce_mrcc_ccsd_sub PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_mrcc_ccsd_sub PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
