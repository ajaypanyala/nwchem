set(SCF_DFT_CG_SOURCE_FILES
dft_get_conv_info.F
dft_roks_hessv_xx.F
dft_roks_energy.F
dft_roks_search_precond.F
dft_uks_nr_solve.F
dft_uks_hessv.F
dft_uks_energy.F
dft_cg_guess.F
dft_cg_answer.F
dft_uks_search_precond.F
dft_cg_solve.F
dft_roks_nr_solve.F
dft_cg_info.F
dft_roks_hessv.F
dft_uks_step_energy.F
dft_uks_hessv_2e.F
dft_roks_step_energy.F
dft_roks_fock.F
dft_roks_hessv_2e.F
dft_cscf.F
)

add_library(scf_dft_cg OBJECT ${SCF_DFT_CG_SOURCE_FILES})
target_include_directories(scf_dft_cg PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(scf_dft_cg PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
