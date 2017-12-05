
set(LRH_SOURCE_FILES
tce_loc_j2.F
tce_loc_j1.F
)

add_library(tce_lrh STATIC ${LRH_SOURCE_FILES})
#add_dependencies(tce_lrh ${NWC_DEPENDENCIES})
target_include_directories(tce_lrh PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_lrh PUBLIC DEBUG_PRINT)

target_compile_options(tce_lrh PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_lrh PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
