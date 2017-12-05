set(NWINTS_ECP_SOURCE_FILES
ecp_3j_prod.f
ecp_angint.f
ecp_angrad.f
ecp_bessel.f
ecp_bess_pow.f
ecp_c2slower.f
ecp_c2srecur.f
ecp_contract.f
ecp_down_m.f
ecp_gaussq.f
ecp_ind_var.f
ecp_s2ccoef.f
ecp_t2_init0.f
ecp_t2_init1.f
ecp_t2_init3.f
ecp_t2_init4.f
ecp_up_k.f
ecp_s2cnorm.F
ecp_bess_asy.F
ecp_int1.F
ecp_t2_p3pow.F
ecp_hessian.F
ecp_c2snorm.F
ecp_t1_Q2int.F
ecp_cart_xpd.F
ecp_t1_Q1pow.F
ecp_local2.F
ecp_t1_Q3pow.F
ecp_skipint.F
ecp_int2_k.F
#test_drive.F
ecp_t1_Q1asy.F
ecp_t1_Q3int.F
ecp_t2_init2.F
ecp_t1_Q0asy.F
ecp_init_con.F
ecp_local1.F
ecp_radint0.F
ecp_gradient.F
ecp_t1_Q3asy.F
ecp_t1_Q0int.F
ecp_init_c2s.F
ecp_grad_a.F
ecp_t1_Q4int.F
ecp_local0.F
ecp_local2_k.F
ecp_t2_ghq.F
ecp_t2_p1pow.F
ecp_locint1.F
ecp_so_ang.F
ecp_grad_b.F
ecp_radint2.F
ecp_int0.F
ecp_matpi.F
ecp_matpr.F
ecp_t1_Q0pow.F
ecp_integral.F
ecp_radint1.F
ecp_t1_Q1int.F
ecp_t2_p4pow.F
ecp_int2.F
ecp_t2_p0pow.F
ecp_cstrans.F
)


add_library(nwints_ecp_un OBJECT ecp_sph_tens.f int_ecp.F)
target_include_directories(nwints_ecp_un PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_ecp_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(nwints_ecp_opt OBJECT ${NWINTS_ECP_SOURCE_FILES})
target_include_directories(nwints_ecp_opt PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_ecp_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwints_ecp_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
