
set(NCC_SOURCE_FILES
./ncc_singles.F
./ncc_doubles.F
./ncc_util.F
./ncc_fourindex.F
)


add_library(ncc_un OBJECT ./ncc_driver.F ./ncc_input.F)
target_include_directories(ncc_un PUBLIC ${NWC_INCLUDES})
target_compile_definitions(ncc_un PUBLIC DETAILED_FREE)
target_compile_options(ncc_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})


add_library(ncc STATIC ${NCC_SOURCE_FILES} $<TARGET_OBJECTS:ncc_un>)
#add_dependencies(ncc ${NWC_DEPENDENCIES})
target_include_directories(ncc PUBLIC ${NWC_INCLUDES})
target_compile_definitions(ncc PUBLIC DETAILED_FREE)
target_compile_options(ncc PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(ncc PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
