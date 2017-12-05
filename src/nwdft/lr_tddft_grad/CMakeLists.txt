set(LR_TDDFT_GRAD_SOURCE_FILES
tddft_grad_make_2pdm.F
tddft_grad_compute_dao.F
tddft_grad_compute_g.F
tddft_grad_print_info.F
tddft_grad_dfxc.F
tddft_grad_solve_precond.F
tddft_grad_store_gradient.F
tddft_grad_root_e.F
tddft_grad_ints.F
tddft_grad_gxc.F
tddft_grad_trans_mo2ao.F
tddft_grad_trans_ao2mo.F
tddft_grad_init.F
tddft_grad_compute_w.F
#tddft_grad_init_func.F
tddft_grad_mat_to_vec.F
tddft_grad_tp.F
tddft_grad_solve_save_x.F
tddft_grad_grad_1e.F
tddft_grad_print_result.F
tddft_grad_mo.F
tddft_grad_load_mo.F
tddft_grad_solve_load_x.F
tddft_gradient.F
tddft_grad_xy.F
tddft_grad_print_array.F
tddft_grad_compute_energy.F
tddft_grad_epsilon.F
tddft_grad_load_xy.F
tddft_grad_get_dens.F
tddft_grad_finalize.F
tddft_grad_r.F
tddft_grad_solve_product.F
tddft_grad_e.F
tddft_grad_util_print.F
tddft_grad_dvxc.F
tddft_grad_anti_symmetrize.F
tddft_grad_compute_p.F
tddft_grad_compute_r.F
tddft_grad_quadv0b.F
tddft_grad_solve.F
tddft_grad_compute_t.F
tddft_grad_input.F
tddft_grad_g.F
tddft_grad_dxcdmat.F
tddft_grad_irreg_acc.F
tddft_grad_fix_exc.F
tddft_grad_mpole.F
tddft_grad_fxc_l.F
tddft_grad_print_header.F
tddft_grad_grad_2e.F
)

add_library(lr_tddft_grad OBJECT ${LR_TDDFT_GRAD_SOURCE_FILES})
target_compile_definitions(lr_tddft_grad PUBLIC DEBUG_PRINT)
target_include_directories(lr_tddft_grad PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(lr_tddft_grad PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(lr_tddft_grad PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
