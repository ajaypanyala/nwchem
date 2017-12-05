
set(QCISD_SOURCE_FILES
qcisd_t1.F
qcisd_t2.F
qcisd_e.F
)

add_library(tce_qcisd STATIC ${QCISD_SOURCE_FILES})
#add_dependencies(tce_qcisd ${NWC_DEPENDENCIES})
target_include_directories(tce_qcisd PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_qcisd PUBLIC DEBUG_PRINT)

target_compile_options(tce_qcisd PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_qcisd PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
