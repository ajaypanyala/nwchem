
set(CCSD_SOURCE_FILES
icsd_t2.F
cc2_t2.F
ccsd_kernels.F
ccsd_1prdm_hh.F
ccsd_e.F
cc2_t1.F
ccsd_1prdm_pp.F
ccsd_t2_8.F
ccsd_1prdm_ph.F
icsd_t1.F
tce_1b_dens_print.F
ccsd_1prdm_hp.F
#sd_t2_8_loops.F
ccsd_t2.F
ccsd_t1.F
ccsd_1prdm.F
)

add_library(tce_ccsd STATIC ${CCSD_SOURCE_FILES})
#add_dependencies(tce_ccsd ${NWC_DEPENDENCIES})
target_include_directories(tce_ccsd PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_ccsd PUBLIC DEBUG_PRINT)


target_compile_options(tce_ccsd PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_ccsd PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
