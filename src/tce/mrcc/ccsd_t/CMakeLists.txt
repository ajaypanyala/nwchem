
set(CCSD_T_SOURCE_FILES
tce_mrcc_q3fnt2.F
tce_mrcc_q3fnt2_subg.F
tce_mrcc_q3vnt1.F
tce_mrcc_q3vnt2.F
tce_mrcc_q3vnt1_subg.F
tce_mrcc_ccsdpt_subg.F
tce_mrcc_ccsdpt.F
)

add_library(tce_mrcc_ccsd_t STATIC ${CCSD_T_SOURCE_FILES})
#add_dependencies(tce_mrcc_ccsd_t ${NWC_DEPENDENCIES})
target_include_directories(tce_mrcc_ccsd_t PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})

target_compile_options(tce_mrcc_ccsd_t PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_mrcc_ccsd_t PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
