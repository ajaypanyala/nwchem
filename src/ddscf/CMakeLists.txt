
set(DDSCF_INCLUDES
${PROJECT_SOURCE_DIR}/src/nwdft/xc
${PROJECT_SOURCE_DIR}/src/nwdft/scf_dft
${PROJECT_SOURCE_DIR}/src/NWints/int
${PROJECT_SOURCE_DIR}/src/NWints/rel
${PROJECT_SOURCE_DIR}/src/NWints/auxil
)

set(DDSCF_UNOPT_SOURCE_FILES
./riscf_dovecs.F
./fock_2e_file.F
./fock_2e_cam.F
./fock_xc.F
./rohf_scat.F
./movecs_lock.F
./movecs_pr_anal.F
./fock_j_fit.F
./rohf_hessv3_ext.F
./rohf_dens.F
./int_dip_ga.F
./fock_get_blk.F
./uhf_precond_ext.F
./ga_atom_blk.F
./rohf.F
./scf_input.F
./rohf_wrap.F
./rohf_jkfac.F
./schwarz_pr.F
./scf_lindep.F
./print_integrals.F
./fock_2e.F
./rohf_k2cf.F
./riscf_io.F
./fock_2e_et.F
./fock_2e_a.F
./rohf_canon.F
./rohf_shift.F
./uhf_hessv2.F
./rohf_precon.F
./rhf_fock_2e.F
./rohf_nr_sol.F
./rohf_hessv2.F
./scf.F
./riscf_list.F
./comp4_insert.F
./rohf_h2e.F
./do_riscf.F
./comp4_bins.F
./rohf_diis.F
./fock_2e_rep.F
./uhf_jkfac.F
./scf_vec_guess.F
./schwarz_init.F
./scf_movecs.F
./vectors_inp.F
./rohf_lagr.F
./rohf_rotate.F
./riscf_fock.F
#./rohf_h1e3.F
./fock_upd_blk.F
./uhf_hessv2_ext.F
./print_integ_list.F
./uhf.F
./rohf_hxx.F
./matrix_exp.F
./int2e_balance.F
./schwarz_tidy.F
./vectors.F
./rohf_fock.F
./schwarz_max.F
./scf_pstat.F
./print_aoints.F
./multipole.F
./rohf_h1e.F
./ga_get2eri.F
./ga_cont_blk.F
./riscf_square.F
./riscf_init.F
./scf_dump_evp_input.F
./scf_sym_scr.F
./ga_orth_mos.F
./scf_frozemb.F
./ga_mat_to_vec.F
./scf_sym_adap.F
./movecs_frag.F
./rhf_dens_mo.F
./print_1eints.F
./fock_2e_cdf.F
./fock_2e_b.F
./mull_pop.F
./int_1e_ga.F
./rohf_fock2g.F
./riscf_tidy.F
./uhf_hessv.F
./riscf_trans.F
./localize.F
./rhf_fock.F
./movecs_phase.F
./scf_get_info.F
./movecs_rotate.F
./rohf_hessv2_ext.F
./rohf_step.F
./movecs_proj.F
./rohf_hessv3.F
./rohf_hscale.F
./ao_fock_2e.F
)

set(DDSCF_OPT_SOURCE_FILES
./fock_2e_slab.F
./fock_2e_lab.F
./fock_2e_task.F
./fock_2e_c.F
./fock_scale.F
./fock_labscal.F
./fock_densc.F
./schwarz_sh.F
./schwarz_atom.F
./comp4_ext.F
./int2e_file.F
./ao_replicated.F

./fast/newtestfmm.F
./fast/cart_trans.F
./fast/testsolver.F
./fast/cheby.F
./fast/anlfit.F
#./fast/fmm.F
./fast/xlm_poles.F
./fast/fastj.F
./fast/newfmm.F
./fast/pot_shell.F
./fast/fastj_direct.F
./fast/solver.F
./fast/potential.F
#./fast/testfmm.F
./fast/fastj_range.f
./fast/gauleg.f
#./fast/numbers.f
./fast/sint.f
#./fast/testcarttrans.f
#./fast/testcheby.f
./fast/vsint.f
)

add_library(ddscf_unopt OBJECT ${DDSCF_UNOPT_SOURCE_FILES})
target_include_directories(ddscf_unopt PUBLIC ${NWC_INCLUDES} ${DDSCF_INCLUDES})
target_compile_definitions(ddscf_unopt PUBLIC RISCF_USE_DRA)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(ddscf_unopt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS} -ffixed-line-length-none)
else()
  target_compile_options(ddscf_unopt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
endif()

add_library(ddscf STATIC ${DDSCF_OPT_SOURCE_FILES} $<TARGET_OBJECTS:ddscf_unopt>)
target_include_directories(ddscf PUBLIC ${NWC_INCLUDES} ${DDSCF_INCLUDES})
target_compile_definitions(ddscf PUBLIC FMM_LMAX=30 RISCF_USE_DRA)
target_compile_options(ddscf PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(ddscf PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
