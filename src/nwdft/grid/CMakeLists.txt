
set(GRID_OPT_SOURCE_FILES
./grid_read_ck.F
./grid_utils.F
./grid_bucket.F
./grid_quadv0b.F
./grid_sswold.F
./grid_quadvc.F
./grid_quadv0a.F
./grid_ylm.F
./grid_quadv0.F
./grid_quadvw1.F
./grid_sphere.F
./grid_file.F
./grid_qpts.F
./grid_pscreen.F
./grid_ssw.F
./grid_quadvw.F
./grid_eaf.F
./grid_sympts.F
./grid_beckew.F
./grid_signf.F
)

set(GRID_SOURCE_FILES
./grid_sic.F
./grid_checks.F
./grid_balance.F
./grid_rtrunc.F
./grid_atom_type_info.F
./grid_setquad.F
./grid_setrad.F
./grid_setang.F
./grid_list.F
./grid_invariant.F
./grid_gen_pts.F
./grid_acc_def.F
./grid_setspac_params.F
./grid_memory.F
./grid_update_rtdb.F
./grid_tasks.F
./grid_setspac.F
./grid_neighb.F
./grid_cleanup.F
./grid_gaussq.F
./grid_finpts.F
./grid_loop.F
./grid_init.F
./dftg_loop.F
./grid_xcfit.F
)


add_library(grid_opt OBJECT ${GRID_OPT_SOURCE_FILES})
target_include_directories(grid_opt PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(grid_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(grid_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

add_library(grid_un OBJECT ${GRID_SOURCE_FILES})
target_include_directories(grid_un PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(grid_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

