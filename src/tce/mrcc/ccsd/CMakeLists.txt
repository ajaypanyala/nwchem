
set(CCSD_SOURCE_FILES
T2_C2.F
T2_C2_1.F
T2_pg_C2.F
T2_pg_C2_1.F
mr_ccsd_t1.F
mr_ccsd_t2.F
mr_ccsd_dcl.F
mr_ccsd_e.F
)

add_library(tce_mrcc_ccsd STATIC ${CCSD_SOURCE_FILES})
#add_dependencies(tce_mrcc_ccsd ${NWC_DEPENDENCIES})
target_include_directories(tce_mrcc_ccsd PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
if(MRCC_METHODS)
  target_compile_definitions(tce_mrcc_ccsd PUBLIC MRCC_METHODS)
endif()
target_compile_options(tce_mrcc_ccsd PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_mrcc_ccsd PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
