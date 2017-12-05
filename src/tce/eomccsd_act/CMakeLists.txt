
set(EOMCCSD_ACT_SOURCE_FILES
eomccsd_x2a.F
q2na_hc1.F
eomccsd_x1a.F
)

add_library(tce_eomccsd_act STATIC ${EOMCCSD_ACT_SOURCE_FILES})
#add_dependencies(tce_eomccsd_act ${NWC_DEPENDENCIES})
target_include_directories(tce_eomccsd_act PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_eomccsd_act PUBLIC DEBUG_PRINT)

target_compile_options(tce_eomccsd_act PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_eomccsd_act PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
