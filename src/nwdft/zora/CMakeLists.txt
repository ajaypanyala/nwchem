set(ZORA_SOURCE_FILES
calc_EFG.F
calc_zora_EFGZ4_SO.F
get_ints_zora.F
dft_zora_EPR.F
zora_tabcd.F
calc_zora_sf.F
zora_getv_EPR.F
calc_zora_facs.F
#grid_zoracorr_FA.F
merf.F
calc_zora_NMRCS_SR.F
zora_getv_NMRCS_SR.F
dft_zora_rhos.F
zora_tabcd_EFG.F
dft_zora_NMR.F
dft_zora_utils.F
dft_zora_utils_so.F
gridHartreePotential.F
grid_zoracorr.F
dft_zora_Hyperfine.F
calc_zora_EFGZ4_SR.F
zora_getv.F
zora_getv_sf.F
dft_zora_EPR-NMR_tools.F
zora_getv_so.F
zora_getv_EFGZ4_SO.F
calc_zora_EPR.F
calc_zora_so.F
gridNuclearPotential.F
calc_zora_HFine.F
zora_getv_EFGZ4_SR.F
zora_recalc.F
dft_zora_rhos_so.F
gridQpqPotential.F
zora_getv_HFine.F
)

add_library(nwdft_zora OBJECT ${ZORA_SOURCE_FILES})
target_include_directories(nwdft_zora PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(nwdft_zora PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwdft_zora PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})