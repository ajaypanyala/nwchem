set(NWDFT_COULOMB_SOURCE_FILES
dft_fitcd.F
dft_get2eri.F
dft_getspm.F
dft_fitvc.F
dft_n3c.F
dft_mem3c.F
dft_3cincor.F
dft_getvc.F
dft_genspm.F
)

add_library(nwdft_coulomb OBJECT ${NWDFT_COULOMB_SOURCE_FILES})
target_include_directories(nwdft_coulomb PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(nwdft_coulomb PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwdft_coulomb PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})