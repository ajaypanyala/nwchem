set(NWINTS_API_SOURCE_FILES
int_lgen1e.F
int_l1e3ov.F
intp_2e3ct.F
int_projpole.F
int_l1eke.F
intpd_2e2c.F
intbd_init4c.F
intdd_init.F
intd_2e4c.F
intdd_2e3c.F
int_mem.F
int_term.F
intp_1e1cpe.F
intdd_1eecp1.F
int_func_c.F
int_dipvel.F
intdd_2e4c.F
int_chk_init.F
int_1cg.F
int_giaol10.F
intd_1epot.F
int_giao_2e.F
int_l2e3c.F
intbdd_init4c.F
intd_1eso.F
intp_mpole.F
int_giaos10.F
intpd_1e1cpe.F
intpd_2e3c.F
intd_1eke.F
int_vstat1e.F
exactd_mem.F
int_1epvpe.F
intpd_mpolel.F
intdd_1eov.F
int_giaoh11.F
intd_mpolel.F
int_1e3ov.F
int_dso.F
int_veloc.F
cando_sp.F
intd_1eov.F
int_l1eh1.F
intbdd_2e4c.F
intp_mpolel.F
int_2e2c.F
intpd_2e3ct.F
int_giaoh01.F
intp_1eke.F
int_mpole.F
int_l1epe.F
intp_txyz.F
int_1estv.F
intdd_2e2c.F
intp_2e3c.F
intpd_1eov.F
int_1er.F
int_giaos100.F
intd_1e3ov.F
int_l1eov.F
int_pso.F
int_nint.F
intpd_1eke.F
int_pgen1e.F
intp_1epe.F
int_init.F
intd_2e3c.F
intp_1eov.F
intd_init.F
int_l1eall.F
int_angmom.F
intp_2e2c.F
int_2e3c.F
int_l2e2c.F
int_1epxvpe.F
intso_1e.F
intd_2e2c.F
intdd_1eh1.F
int_giaotv10.F
int_chk_sh.F
int_canon.F
intbd_2e4c.F
intd_1eh1.F
int_1eelec.F
int_mpolel.F
int_giaobq10.F
intd_1epot_cosmo.F
int_acc.F
intpd_1epe.F
exact_mem.F
)

set(NWINTS_API_OPT_SOURCE_FILES
int_hf1sp.F
int_l2e4c.F 
intb_2e4c.F
intb_init4c.F
int_2e4c.F
)

if(USE_TEXAS)
	set(API_COMPILE_DEFS USE_TEXAS)
endif()



add_library(nwints_api_un OBJECT ${NWINTS_API_SOURCE_FILES})
target_compile_definitions(nwints_api_un PUBLIC ${API_COMPILE_DEFS})
target_include_directories(nwints_api_un PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_api_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(nwints_api_opt OBJECT ${NWINTS_API_OPT_SOURCE_FILES})
target_compile_definitions(nwints_api_opt PUBLIC ${API_COMPILE_DEFS})
target_include_directories(nwints_api_opt PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_api_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwints_api_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
