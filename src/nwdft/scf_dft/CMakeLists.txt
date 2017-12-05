set(SCF_DFT_SOURCE_FILES
dft_pstat.F
dft_diagn.F
dft_fockbld.F
dft_frac_mo.F
dft_nwmolden.F
dft_symbr.F
dft_canorg.F
dft_dump_info.F
dft_swapab.F
dft_scfutil.F
fukui.F
print_irrep_occp.F
dft_scf.F
dft_fant.F
dft_main0d.F
dft_printinfo.F
dft_fock_efield.F
)

set(SCF_DFT_OPT_SOURCE_FILES
dft_dencvg.F
dft_densm.F
dft_invdiag.F
dft_mxovl.F
dft_mxspin_ovlp.F
dft_s2_value.F
dft_scfcvg.F
diis_bld12.F
diis_bsolve.F
diis_core.F
diis_driver.F
diis_hamwgt.F
dft_fdist_init.F
cdft_util.F
)

add_library(scf_dft_un OBJECT ${SCF_DFT_SOURCE_FILES})
target_include_directories(scf_dft_un PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(scf_dft_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(scf_dft_opt OBJECT ${SCF_DFT_OPT_SOURCE_FILES})
target_include_directories(scf_dft_opt PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(scf_dft_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(scf_dft_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
