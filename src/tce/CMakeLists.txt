set(TCE_INCLUDES 
${PROJECT_SOURCE_DIR}/src/perfm
${PROJECT_SOURCE_DIR}/src/tce
${PROJECT_SOURCE_DIR}/src/tce/ccsdt_yr
${PROJECT_SOURCE_DIR}/src/tce/ccsdtq
${PROJECT_SOURCE_DIR}/src/tce/lrh
${PROJECT_SOURCE_DIR}/src/tce/ccsdt
${PROJECT_SOURCE_DIR}/src/tce/eomccsd_act
${PROJECT_SOURCE_DIR}/src/tce/response
${PROJECT_SOURCE_DIR}/src/tce/mrcc
${PROJECT_SOURCE_DIR}/src/tce/mrcc/ccsd
${PROJECT_SOURCE_DIR}/src/tce/mrcc/include
${PROJECT_SOURCE_DIR}/src/tce/mrcc/ccsd_sub
${PROJECT_SOURCE_DIR}/src/tce/mrcc/ccsd_t
${PROJECT_SOURCE_DIR}/src/tce/eaccsd
${PROJECT_SOURCE_DIR}/src/tce/ccsd
${PROJECT_SOURCE_DIR}/src/tce/mbpt4
${PROJECT_SOURCE_DIR}/src/tce/emb
${PROJECT_SOURCE_DIR}/src/tce/cisdt
${PROJECT_SOURCE_DIR}/src/tce/ccsdtq_lambda
${PROJECT_SOURCE_DIR}/src/tce/ccsd_alpha
${PROJECT_SOURCE_DIR}/src/tce/cisd
${PROJECT_SOURCE_DIR}/src/tce/ccsd_yr
${PROJECT_SOURCE_DIR}/src/tce/gradients
${PROJECT_SOURCE_DIR}/src/tce/mbpt2
${PROJECT_SOURCE_DIR}/src/tce/ccd
${PROJECT_SOURCE_DIR}/src/tce/diis3
${PROJECT_SOURCE_DIR}/src/tce/cis
${PROJECT_SOURCE_DIR}/src/tce/ccsdt_alpha_new
${PROJECT_SOURCE_DIR}/src/tce/ccsd_lambda
${PROJECT_SOURCE_DIR}/src/tce/cisdtq
${PROJECT_SOURCE_DIR}/src/tce/crop
${PROJECT_SOURCE_DIR}/src/tce/ccsd_qr
${PROJECT_SOURCE_DIR}/src/tce/sort
${PROJECT_SOURCE_DIR}/src/tce/include
${PROJECT_SOURCE_DIR}/src/tce/eomccsdtq
${PROJECT_SOURCE_DIR}/src/tce/ccsdtq_alpha_offdiag
${PROJECT_SOURCE_DIR}/src/tce/mbpt3
${PROJECT_SOURCE_DIR}/src/tce/eomccsd
${PROJECT_SOURCE_DIR}/src/tce/cr_eomccsd_t_act
${PROJECT_SOURCE_DIR}/src/tce/lccd
${PROJECT_SOURCE_DIR}/src/tce/cc3
${PROJECT_SOURCE_DIR}/src/tce/qcisd
${PROJECT_SOURCE_DIR}/src/tce/ccsdt_beta
${PROJECT_SOURCE_DIR}/src/tce/diis
${PROJECT_SOURCE_DIR}/src/tce/ccsdt_lambda
${PROJECT_SOURCE_DIR}/src/tce/mbpt1
${PROJECT_SOURCE_DIR}/src/tce/eomccsdt
${PROJECT_SOURCE_DIR}/src/tce/ipccsd
${PROJECT_SOURCE_DIR}/src/tce/diis2
${PROJECT_SOURCE_DIR}/src/tce/cr-eomccsd_t
${PROJECT_SOURCE_DIR}/src/tce/lccsd
${PROJECT_SOURCE_DIR}/src/tce/ccsd_t
${PROJECT_SOURCE_DIR}/src/tce/ccsd_beta
${PROJECT_SOURCE_DIR}/src/tce/ccsd_act
${PROJECT_SOURCE_DIR}/src/tce/ccsdt_alpha_offdiag
)

if(TCE_CUDA AND NOT CCSD_T_GPU_OLD)  
  add_subdirectory(ttlg)
endif()

add_subdirectory(lrh)
add_subdirectory(eomccsd_act)
add_subdirectory(ccsd)
add_subdirectory(mbpt4)
add_subdirectory(emb)
add_subdirectory(cisdt)
add_subdirectory(cisd)
add_subdirectory(gradients)
add_subdirectory(mbpt2)
add_subdirectory(ccd)
add_subdirectory(diis3)
add_subdirectory(cis)
add_subdirectory(ccsd_lambda)
add_subdirectory(cisdtq)
add_subdirectory(ccsd_qr)
add_subdirectory(sort)
add_subdirectory(mbpt3)
add_subdirectory(eomccsd)
add_subdirectory(cr_eomccsd_t_act)
add_subdirectory(lccd)
add_subdirectory(cc3)
add_subdirectory(qcisd)
add_subdirectory(diis)
add_subdirectory(mbpt1)
add_subdirectory(diis2)
add_subdirectory(cr-eomccsd_t)
add_subdirectory(lccsd)
add_subdirectory(ccsd_t)
add_subdirectory(ccsd_act)


# ifeq ($(NWCHEM_TARGET),$(findstring $(NWCHEM_TARGET),LAPI64 IBM64))
#       OBJ += tce_eom_xguess.o tce_eom_yguess.o
# else
#       OBJ_OPTIMIZE += tce_eom_xguess.o tce_eom_yguess.o
# endif



SET(TCE_COMPILE_DEFS "")
SET(LIBTCE_OPT "")

if(MA_TRANS_BLOCKED)
    list(APPEND TCE_COMPILE_DEFS MA_TRANS_BLOCKED)
endif()

if(EACCSD)
    add_subdirectory(eaccsd)
    set(LIBTCE_OPT ${LIBTCE_OPT} tce_eaccsd)
    list(APPEND TCE_COMPILE_DEFS EACCSD)
endif()

if(IPCCSD)
    add_subdirectory(ipccsd)
    set(LIBTCE_OPT ${LIBTCE_OPT} tce_ipccsd)
    list(APPEND TCE_COMPILE_DEFS IPCCSD)
endif()

#Set to ON by default for backward compatibility
option(CCSDT "CCSDT" ON)
option(CCSDLR "CCSDTLR" ON)

if(CCSDTLR)
    set(LIBTCE_OPT ${LIBTCE_OPT} tce_ccsdt_yr tce_ccsdt_beta tce_ccsdt_alpha_new tce_ccsdt_alpha_offdiag)
    add_subdirectory(ccsdt_yr)
    add_subdirectory(ccsdt_beta)
    add_subdirectory(ccsdt_alpha_new)
    add_subdirectory(ccsdt_alpha_offdiag)
    list(APPEND TCE_COMPILE_DEFS CCSDTLR)
endif()

if(CCSDT)
   add_subdirectory(ccsdt)
   add_subdirectory(eomccsdt)
   add_subdirectory(ccsdt_lambda)
   list(APPEND TCE_COMPILE_DEFS CCSDT)
endif()
 

if(CCSDLR) 
    set(LIBTCE_OPT ${LIBTCE_OPT} tce_response tce_ccsd_yr tce_ccsd_alpha tce_ccsd_beta)
    add_subdirectory(response)
    add_subdirectory(ccsd_yr)
    add_subdirectory(ccsd_alpha)
    add_subdirectory(ccsd_beta)
    list(APPEND TCE_COMPILE_DEFS CCSDLR)
endif()

if(CCSDTQ) 
    set(LIBTCE_OPT ${LIBTCE_OPT} tce_response tce_ccsdtq tce_eomccsdtq tce_ccsdtq_lambda tce_ccsdtq_alpha_offdiag)
    add_subdirectory(response)
    add_subdirectory(ccsdtq)
    add_subdirectory(eomccsdtq)
    add_subdirectory(ccsdtq_lambda)
    add_subdirectory(ccsdtq_alpha_offdiag)
    list(APPEND TCE_COMPILE_DEFS CCSDTQ)
endif()

if(MRCC_METHODS)  
    set(LIBTCE_OPT ${LIBTCE_OPT} tce_mrcc)
    add_subdirectory(mrcc)
    list(APPEND TCE_COMPILE_DEFS MRCC_METHODS)
endif()

if(TCE_CUDA)  
    list(APPEND TCE_COMPILE_DEFS TCE_CUDA)
endif()

if(USE_PSTAT) 
    list(APPEND TCE_COMPILE_DEFS USE_PSTAT)
endif()

set(LIBTCE_OPT ${LIBTCE_OPT} ${TCE_GPU_LIBS})
set(TCE_CUDA_LIBS ${TCE_CUDA_LIBS} PARENT_SCOPE)

#Set parent scope variable just once, cannot append to it
SET(LIBTCE_OPT ${LIBTCE_OPT} PARENT_SCOPE)

set(TCE_SOURCE_FILES
tce_x1_offset.F
tce_tile.F
tce_property.F
tce_ao2e_cholesky_2d.F
ma_pack.F
add_block.F
get_mo_rdm_pp.F
sf_print.F
tce_create_x1.F
tce_clone_y4.F
tce_create_t3.F
tce_dgeev.F
deletefile.F
tce_mo2e_zones_4a_disk_4s_new_patch_m3.F
tce_multipole_create.F
tce_x2_offset_new.F
tce_ccsdtq_driver.F
tce_t1_offset_new.F
tce_guess_t2.F
reconcilefile.F
tce_eom_ydiagon.F
tce_reference_check.F
tce_create_y1.F
tce_mo1e.F
tce_y1_offset_new.F
tce_mo2e_hybrid_2eorb_split.F
tce_mo2e_zones_4a_disk_2s_new.F
tce_ccsd_lambda_perturbative.F
tce_y2_offset.F
createfile.F
tce_mo2e_zones_4a_disk_ga_act.F
tce_ccsdt_driver.F
tce_mo2e_zones_4a_disk_4s_new_patch.F
tce_clone_x2.F
ccsd_energy_loc.F
tce_mo2e_2emet4.F
tce_clone_y2.F
tce_t4_offset.F
tce_eom_xguess.F
tce_y4_offset.F
tce_mo2e_zones_4a_disk_ga_N5.F
tce_aod1.F
tce_mo2e_zones_4a_disk_ga_chop_act.F
tce_fock_destroy.F
tce_mo2e_offset_intorb.F
tce_mo1e_offset.F
#ccsd_imaginary.F
tce_restart.F
tce_mo2e_zones_4a_disk_ga_chop.F
tce_ao2e_cholesky_4d.F
tce_mo2e_zones_4a_disk_4s_new_patch_m2.F
tce_diagnose_t1.F
tce_ccsd_hyperpolar.F
tce_create_x4.F
get_mo_rdm_hp.F
tce_mod1.F
ccsd_lambda.F
get_mo_rdm_hh.F
tce_print_x2.F
tce_eom_init.F
tce_dens_ph_offset.F
tce_create_y4.F
tce_t4_offset_new.F
put_block.F
tce_ccsd_asym_polar.F
tce_restricted.F
tce_create_t4.F
ma_print_compact.F
tce_ao2e.F
tce_x3a_offset.F
tce_mo2e_offset_intorb_hybrid.F
tce_t3a_offset.F
tce_hash.F
tce_ccsdt_lambda_driver.F
ga_test.F
tce_mo2e_new.F
tce_xy1_offset.F
tce_dens_pp_offset.F
tce_mo2e_trans.F
get_block_ind.F
tce_ao2e_disk.F
atpart_nalength.F
tce_nud1.F
tce_filename.F
tce_ccsd_lambda_driver.F
tce_fourindex_cleanup.F
daxfile.F
tce_y4_offset_new.F
tce_mo2e_offset_size.F
tce_fourindex_driver.F
get_mo_rdm_ph.F
tce_x4_offset_new.F
tce_e_offset.F
tce_mo2e_incore_2eorb_split.F
tce_y3_offset.F
tce_create_x3.F
tce_xy2_offset.F
tce_zero.F
tce_t2_offset.F
tce_ccsd_response_driver.F
tce_residual_t1.F
tce_ccsd_lambda_response_driver.F
tce_mo2e.F
tce_mod1_new.F
tce_ao2e_cholesky.F
tce_create_y2.F
tce_guess_t1.F
daxpyfile.F
#tce_create_movecs.F
tce_mo2e_2emet4_n5.F
tce_fock_create.F
tce_ao1e.F
tce_multipole_destroy.F
tce_prop_start.F
tce_create_x2.F
hbar.F
tce_ao2e_new.F
tce_x3_offset_new.F
tce_mo2e_2emet3_n5.F
tce_mo2e_2emet3.F
tce_clone_y3.F
tce_create_y3.F
tce_clone_y1.F
tce_mod1_offset.F
tce_eom_tidy.F
tce_x2a_offset.F
get_block.F
get_hash_block_ma.F
ddotfile.F
tce_t1_offset.F
tce_mbpt2.F
sf_test.F
tce_t2a_offset.F
tce_clone_x3.F
gatoeaf.F
tce_residual_t2.F
dratoga.F
tce_t3_offset.F
tce_print_x1.F
copyfile.F
dscalfile.F
tce_hash_ind.F
tce_t2_offset_new.F
tce_t3_offset_new.F
ccsd_qr_beta.F
ao_rdm_write.F
tce_y1_offset.F
tce_clone_x4.F
tce_init_new.F
tce_mo2e_incore_2eorb.F
tensor_read_write.F
tce_residual_t3a.F
tce_schmidt.F
tce_restart_quadruples.F
get_block_ind_sf.F
tce_create_t2.F
tce_ccsd_perturbative.F
add_hash_block.F
tce_eom_xdiagon.F
tce_x3_offset.F
tce_mo2e_disk_2eorb.F
tce_ccsdtq_lambda_driver.F
tce_y2_offset_new.F
tce_x2_offset.F
tce_dens_hp_offset.F
tce_ccsd_symm_polar.F
tce_x1_offset_new.F
tce_clone_x1.F
tce_ccsd_driver.F
tce_tidy.F
tce_residual_t4.F
tce_cholesky.F
tce_residual_t3.F
#tce_ao2e_cholesky_orb.F
e2_offset_2_offset.F
tce_restart_triples.F
tce_ccsdtq_response_driver.F
ccsd_lr_alpha.F
get_hash_block.F
tce_mo2e_zones_4a_disk_ga_chop_N5.F
tce_mo2e_disk.F
tce_diis_act.F
tce_create_t1.F
ma_scale.F
tce_x4_offset.F
put_hash_block.F
tce_ccsdt_symm_polar.F
tce_tile_new.F
tce_ccsdt_response_driver.F
new_ga4ind_N5.F
new4ind_N5_osu.F
tce_ccsd_hyperpolar_static.F
tce_eom_yguess.F
tce_ccsdtq_symm_polar.F
tce_y3_offset_new.F
tce_mod1_offset_new.F
gatodra.F
sym_abelian_axis.F
cc2_energy.F
tce_mo2e_offset.F
tce_dens_hh_offset.F
tce_mo2e_zones_4a_disk_ga.F
)

if(MA_TRANS_BLOCKED)
    list(APPEND TCE_SOURCE_FILES ./ma_transpose_blocked.F)
else()
    list(APPEND TCE_SOURCE_FILES ./ma_transpose.F)
endif()


add_library(tce_un OBJECT tce_input.F tce_energy.F tce_energy_fragment.F tce_init.F)
target_include_directories(tce_un PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_un PUBLIC ${TCE_COMPILE_DEFS})
target_compile_options(tce_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(tce STATIC ${TCE_SOURCE_FILES} $<TARGET_OBJECTS:tce_un>)
#add_dependencies(tce ${NWC_DEPENDENCIES})
target_include_directories(tce PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce PUBLIC ${TCE_COMPILE_DEFS})
target_compile_options(tce PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

