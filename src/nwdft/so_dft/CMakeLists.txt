set(SO_DFT_SOURCE_FILES
dft_mxovlso.F
getv_coul_cd.F
dft_densmso.F
dft_frac_mo_so.F
sym_mo_ap_op_so.F
movecs_anal_so.F
vectors_so.F
getv_coul.F
sym_op_cname_so.F
dft_scf_so.F
grad_dens_so.F
getv_exchY.F
getv_exch0.F
dft_scf_utils_so.F
scf_sym_adap_so.F
ga_orth_mos_so.F
diis_hamwgt_so.F
xc_exso.F
grad_dens_so2.F
sym_op_clsfy_so.F
diis_bsolve_so.F
grad_force_so.F
sym_char_tab_so.F
diis_bld12_so.F
getv_exchX.F
diis_driver_so.F
diis_core_so.F
getv_exchZ.F
getv_dftxc.F
diis_refvec_so.F
dft_gradients_so.F
sym_mo_adapt_so.F
grad1_so.F
ga_orthog_so.F
)

add_library(so_dft OBJECT ${SO_DFT_SOURCE_FILES})
target_include_directories(so_dft PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(so_dft PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(so_dft PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
