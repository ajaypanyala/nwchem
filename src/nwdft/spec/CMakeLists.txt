set(NWDFT_SPEC_SOURCE_FILES
dft_vspec.F
dft_xps.F
dft_vspec_so.F
)

add_library(nwdft_spec OBJECT ${NWDFT_SPEC_SOURCE_FILES})
target_include_directories(nwdft_spec PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(nwdft_spec PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwdft_spec PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})