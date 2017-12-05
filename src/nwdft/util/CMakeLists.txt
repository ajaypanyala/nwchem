set(DFT_UTIL_SOURCE_FILES
dft_scatter.F
dft_mpole.F
nbf_to_mbf.F
movecs_log.F
dft_tstamp.F
dft_genutils.F
dft_pstrace.F
sca_lab_mat.F
dft_frozemb.F
dft_gather.F
dft_header.F
dft_utilmap.F
dft_frozemb_fde.F
gacrea.F
get_col.F
int_1e_fde.F
)

add_library(dft_util_un OBJECT ga_chol.F ga_inverse.F ga_chol_seq.F ga_inv_seq.F 
			       dft_invio.F dft_print_mats.F dft_save_mats.F)
target_include_directories(dft_util_un PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(dft_util_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(dft_util OBJECT ${DFT_UTIL_SOURCE_FILES})
target_include_directories(dft_util PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(dft_util PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(dft_util PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
