
set(CCSD_YR_SOURCE_FILES
ccsd_yr1.F
ccsd_lr_d1_hp.F
ccsd_lr_d1_ph.F
ccsd_y_tr1.F
ccsd_yr2.F
ccsd_y_o2.F
ccsd_lr_d1_pp.F
ccsd_lr_d1_hh.F
ccsd_y_o1.F
ccsd_y_tr2.F
alpha_2_yr.F
)

add_library(tce_ccsd_yr STATIC ${CCSD_YR_SOURCE_FILES})
#add_dependencies(tce_ccsd_yr ${NWC_DEPENDENCIES})
target_include_directories(tce_ccsd_yr PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_ccsd_yr PUBLIC DEBUG_PRINT)

target_compile_options(tce_ccsd_yr PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_ccsd_yr PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
