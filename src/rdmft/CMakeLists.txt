
set(RDMFT_SOURCE_FILES
./wfn1.F
./wfn1_calc_occ.F
./wfn1_lindep.F
./wfn1_f0df0f1.F
./wfn1_exp.F
./wfn1_dnrm2.F
./wfn1_configuration.F
./wfn1_print.F
./wfn1_2el_wfn1.F
./wfn1_allocate_mat.F
./wfn1_f1f2f3.F
./wfn1_2el_hf.F
./wfn1_zero_vv.F
./wfn1_size.F
./wfn1_test_rdms.F
./wfn1_trace.F
./wfn1_idmin.F
./wfn1_closest_three.F
./wfn1_guess_corr.F
./wfn1_fix_xm.F
./wfn1_df1df2df3.F
./wfn1_e_g.F
./wfn1_entropy.F
./wfn1_linesearch.F
./wfn1_insert.F
./wfn1_gen_trial_wfn.F
./wfn1_gen_mat.F
./wfn1_guess_movecs.F
./wfn1_allocate_drv.F
./wfn1_print_occ.F
./wfn1_e_g_2el.F
./wfn1_e_g_2el_hf.F
./wfn1_e_g_2el_s.F
./wfn1_print_analysis.F
./wfn1_fix_phase.F
./wfn1_allocate_wfn.F
./wfn1_ortho.F
)

add_library(rdmft STATIC ${RDMFT_SOURCE_FILES})
#add_dependencies(rdmft ${NWC_DEPENDENCIES})
target_include_directories(rdmft PUBLIC ${NWC_INCLUDES})
target_compile_options(rdmft PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(rdmft PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

