
set(MM_SOURCE_FILES
./mm_bond.F
./mm_bq.F
./mm_neb.F
./mm_vdw_coords.F
./mm_cons.F
./mm_input.F
./mm_geom.F
./mm_vdw.F
./mm_bond_coords.F
./mm_link.F
./mm_coords.F
./mm_data.F
)

add_library(mm STATIC ${MM_SOURCE_FILES})
#add_dependencies(mm ${NWC_DEPENDENCIES})
target_include_directories(mm PUBLIC ${NWC_INCLUDES})
target_compile_options(mm PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(mm PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
