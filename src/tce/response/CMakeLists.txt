
set(RESPONSE_SOURCE_FILES
ccsd_o1.F
ccsd_o2.F
eomccsd_o2.F
ccsdtq_o1.F
eomccsd_o1.F
ccsdtq_o3.F
ccsdt_o1.F
ccsdt_o3.F
ccsdt_o2.F
ccsdtq_o4.F
ccsdtq_o2.F
)

add_library(tce_response STATIC ${RESPONSE_SOURCE_FILES})
#add_dependencies(tce_response ${NWC_DEPENDENCIES})
target_include_directories(tce_response PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_response PUBLIC DEBUG_PRINT)

target_compile_options(tce_response PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_response PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
