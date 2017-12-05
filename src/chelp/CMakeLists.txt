
set(CHELP_SOURCE_FILES
./task_chelp.F
./grid/chelp_grid_data.F
./grid/chelp_grid.F
./grid/chelp_atom_rad.F
)

add_library(chelp STATIC ${CHELP_SOURCE_FILES})
#add_dependencies(chelp ${NWC_DEPENDENCIES})
target_include_directories(chelp PUBLIC ${NWC_INCLUDES})
target_compile_options(chelp PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(chelp PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

