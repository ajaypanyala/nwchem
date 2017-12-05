
set(CCSD_QR_SOURCE_FILES
ccsd_quad_xy1.F
ccsd_quad_ot1.F
ccsd_quad_ot2.F
ccsd_quad_xy2.F
)

add_library(tce_ccsd_qr STATIC ${CCSD_QR_SOURCE_FILES})
#add_dependencies(tce_ccsd_qr ${NWC_DEPENDENCIES})
target_include_directories(tce_ccsd_qr PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_ccsd_qr PUBLIC DEBUG_PRINT)

target_compile_options(tce_ccsd_qr PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_ccsd_qr PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
