set(XC_SOURCE_FILES
xc_hfexch.F
xc_ssbD_1.F
xc_cm11.F
xc_spbe96.F
xc_fd.F
xc_signf.F
xc_op.F
xc_eval_fnl.F
xc_vwn_2.F
ts_eval_fnl.F
xc_sic.F
steric_energy.F
xc_oep.F
ts_vw.F
eumac.F
xc_3rd_deriv.F
xc_murakn.F
a_dist.F
xc_kop.F
xc_bnl.F
xc_becke88.F
xc_s12x.F
xc_xpkzb99.F
xc_lb94.F
xc_tarad.F
xc_xsogga.F
xc_hcth.F
xc_cm06.F
xc_cams12x.F
xc_dirac.F
xc_optc.F
xc_cvs98.F
xc_cm05.F
qdist.F
xc_xvs98.F
xc_camxpbe96.F
xc_xn12.F
xc_cpkzb99.F
xc_xm11.F
xc_vwn_4.F
xc_xmpw91.F
xc_vwn_3.F
xc_ft97.F
xc_xpbe96.F
xc_lyp88.F
xc_xdm.F
xc_xmvs15.F
xc_cMpbe96.F
xc_exso.F
xc_gill96.F
xc_att_xc.F
xc_bc95.F
xc_vwn_1_rpa.F
xc_tabcd.F
xc_camxlsd.F
xc_xpw86.F
xc_vwn_5.F
xc_perdew86.F
xc_csogga.F
xc_xpw6.F
xc_pbe96.F
xc_xtpss03.F
scat_mat.F
xc_xb86b.F
xc_perdew91.F
xc_ctpss03.F
xc_vdw_util.F
xc_pw91ldag.F
xc_pw91lda.F
localize_sic.F
xc_lebedev.F
xc_replicated.F
xc_vwn_1.F
xc_eval_basis.F
xc_rhogen.F
xc_vdw_main.F
xc_optx.F
xc_xwpbe.F
xc_xm05.F
xc_dirac_screened.F
ts_tf.F
xc_fitv.F
xc_atomblock.F
xc_kt1.F
xc_sic_pert.F
xc_camb88.F
setACmat.F
xc_xm06.F
xc_cs00.F
hse08fx.F
xc_xpw91.F
)

set(XC_UNOPT_SOURCE_FILES
xc_getv.F 
xc_getvxc.F
xc_exact_pot.F
xc_sicutil.F
chk_dens.F
xc_util.F
xc_chkdispauto.F
xc_vdw_init.F
)

add_library(nwdft_xc_c OBJECT 2ndDerivB97.c)
target_include_directories(nwdft_xc_c PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(nwdft_xc_c PUBLIC ${NWCHEM_C_COMPILE_FLAGS} ${NWCHEM_C_DEBUG_FLAGS})
#target_compile_options(nwdft_xc_c PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})

add_library(nwdft_xc_un OBJECT ${XC_UNOPT_SOURCE_FILES})
target_include_directories(nwdft_xc_un PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(nwdft_xc_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(nwdft_xc OBJECT ${XC_SOURCE_FILES})
target_include_directories(nwdft_xc PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(nwdft_xc PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwdft_xc PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
