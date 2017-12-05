
set(CISDT_SOURCE_FILES
cisdt_c2.F
cisdt_c1.F
cisdt_c3.F
cisdt_e.F
)

add_library(tce_cisdt STATIC ${CISDT_SOURCE_FILES})
#add_dependencies(tce_cisdt ${NWC_DEPENDENCIES})
target_include_directories(tce_cisdt PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_cisdt PUBLIC DEBUG_PRINT)

target_compile_options(tce_cisdt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_cisdt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
