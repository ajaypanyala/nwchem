
set(CONS_SOURCE_FILES
./cons_springs.F
./cons_utils.F
./cons_rtdb.F
./cons_data.F
./cons_input.F
./cons.F
)

add_library(cons STATIC ${CONS_SOURCE_FILES})
#add_dependencies(cons ${NWC_DEPENDENCIES})
target_include_directories(cons PUBLIC ${NWC_INCLUDES})
target_compile_options(cons PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(cons PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
