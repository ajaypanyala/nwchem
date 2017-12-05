
set(USS_SOURCE_FILES
tce_uss_offdiagonal_4.F
tce_uss_T2_R2.F
ccsd_r3.F
ccsd_r4.F
tce_uss_offdiagonal_1.F
tce_uss_offdiagonal_3.F
tce_uss_T1_R1.F
tce_uss_aux.F
tce_uss_offdiagonal_2.F
)

add_library(tce_mrcc_uss STATIC ${USS_SOURCE_FILES})
#add_dependencies(tce_mrcc_uss ${NWC_DEPENDENCIES})
target_include_directories(tce_mrcc_uss PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
if(MRCC_METHODS)
  target_compile_definitions(tce_mrcc_uss PUBLIC MRCC_METHODS)
endif()
target_compile_options(tce_mrcc_uss PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_mrcc_uss PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
