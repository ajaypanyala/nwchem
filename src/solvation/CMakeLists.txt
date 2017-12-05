
set(SOLVATION_SOURCE_FILES
./mnsol_interface.F
./mnsol_mem.F
./mnsol.F
./solv_data.F
./cosmo_smd_output.F
./cosmo_vem_output.F
)

set(SOLVATION_OPT_SOURCE_FILES
./cosmo_initialize.F
./cosmo_input.F
./cosmo_tidy.F
./cosmo_def_radii.F
./cosmo_charges.F
./hnd_cosmo_lib.F
./grad_hnd_cos.F
)

add_library(solvation_opt OBJECT ${SOLVATION_OPT_SOURCE_FILES})
target_include_directories(solvation_opt PUBLIC ${NWC_INCLUDES})
target_compile_options(solvation_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(solvation_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})


add_library(solvation STATIC ${SOLVATION_SOURCE_FILES} $<TARGET_OBJECTS:solvation_opt>)
#add_dependencies(solvation ${NWC_DEPENDENCIES})
target_include_directories(solvation PUBLIC ${NWC_INCLUDES})
target_compile_options(solvation PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
