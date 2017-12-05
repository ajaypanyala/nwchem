
set(EOMCCSD_SOURCE_FILES
eomccsd_density1.F
eomccsd_1prdm_hp.F
eomccsd_x2.F
cc2_x2.F
eomicsd_x2.F
eomccsd_y2.F
eomicsd_x1.F
eomccsd_1prdm_pp.F
eomccsd_1prdm_hh.F
cc2_x1.F
eomccsd_denominator.F
eomccsd_x1.F
eomccsd_1prdm.F
eomccsd_x2_old.F
eomccsd_1prdm_ph.F
eomccsd_x1_old.F
eomccsd_y1.F
)

add_library(tce_eomccsd STATIC ${EOMCCSD_SOURCE_FILES})
#add_dependencies(tce_eomccsd ${NWC_DEPENDENCIES})
target_include_directories(tce_eomccsd PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_eomccsd PUBLIC DEBUG_PRINT)
target_compile_options(tce_eomccsd PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_eomccsd PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
