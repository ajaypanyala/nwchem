
set(ETRANS_SOURCE_FILES
./et_vectors_inp.F
./et_omega2_mult.F
./et_dgefa.F
./et_fmo_util.F
./et_2e.F
./et_det.F
./et_calc_fmo.F
./et_movecs_read.F
./et_calc.F
./et_input.F
./et_fock.F
./et_dens_prod.F
)

add_library(etrans STATIC ${ETRANS_SOURCE_FILES})
#add_dependencies(etrans ${NWC_DEPENDENCIES})
target_include_directories(etrans PUBLIC ${NWC_INCLUDES})
target_compile_options(etrans PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(etrans PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
