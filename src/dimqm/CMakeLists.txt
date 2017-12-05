set(DIMQM_INCLUDES
${PROJECT_SOURCE_DIR}/src/perfm
${PROJECT_SOURCE_DIR}/src/nwdft/xc
${PROJECT_SOURCE_DIR}/src/nwdft/grid
${PROJECT_SOURCE_DIR}/src/nwdft/include
${PROJECT_SOURCE_DIR}/src/nwdft/scf_dft
${PROJECT_SOURCE_DIR}/src/NWints/int
${PROJECT_SOURCE_DIR}/src/NWints/hondo
)


set(DIMQM_SOURCE_FILES
./matvecComplex.F
./dimqm_broadcast.F
./dimqm_util.F
./dimqm_addop_uhf.F
./dimqm_EqmE.F
./seeded_qmr_real_augment.F
./qmr_complex.F
./qmr_seed_complex.F
./dim_grid_quadv0b.F
./dim_getv.F
./dimqm_prep.F
./dimqm_screening.F
./matvecReal.F
./dimqm_seeding.F
./CubicSpline.F
./pphilo.F
./dim_grid_quadv0.F
./dimqm_printDipole.F
./dimqm_rohf.F
./dimqm_polarizability.F
./dim_grid_quadv0a.F
./dimqm_iterative.F
./dim_dervxyz.F
./dimqm_fldpot.F
./dim_elfcon.F
./dim_tabcd.F
./dim_eval_fnl.F
./t2.F
./dim_elfder.F
./seeded_qmr_complex.F
./dimpar_input.F
./dim_fock_xc.F
./dimqm_indDipoles.F
./seeded_qmr_real.F
./dimqm_input.F
./dim_sxyz2.F
./dimqm_f2d.F
./dim_1eelec.F
./dimqm_Dens2Dipoles.F
#./dimqm_addop_cmplx.F
./qmr_seed_real.F
./dimqm_printAtomicDipoles.F
./dimqm_energy.F
./dimqm_main.F
./qmr_real.F
./dimqm_EqmN.F
./dimqm_tidy.F
./dim_sxyz.F
./dimqm_printAtoms.F
./dimqm_addop.F
./dimqm_lclfld.F
)

add_library(dimqm STATIC ${DIMQM_SOURCE_FILES})
#add_dependencies(dimqm ${NWC_DEPENDENCIES})
target_include_directories(dimqm PUBLIC ${NWC_INCLUDES} ${DIMQM_INCLUDES})

# If you add -DJEFF_OUTPUT, the dipole polarizability will be printed in the same format as the TCE.
# This is very useful (1) if you wrote scripts to parse that particular formatting and (2) if you want the anisotrophy.
target_compile_definitions(dimqm PUBLIC DBLEPREC_RYS)
target_compile_options(dimqm PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(dimqm PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

