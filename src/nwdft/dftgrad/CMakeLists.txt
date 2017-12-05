set(DFT_GRAD_SOURCE_FILES
dft_grad_top.F
dftg_gridv0.F
dft_gradients.F
dftg_getxc.F
grad_dft.F
dftg_gridv0b.F
dftg_getxden.F
dftg_gridv0a.F
dftg_gridpart.F
dftg_griddo.F
dftg_getcden.F
)

set(DFT_GRAD_OPT_SOURCE_FILES
dftg_cdfit.F
dftg_quadv0b.F
dftg_grad2.F
dint_block_j.F
j_nucder_gen.F
dftg_beckew.F
dftg_wderiv.F
dftg_signf.F
dftg_gridssw.F
dftg_gridsswold.F
)

add_library(dft_grad_un OBJECT ${DFT_GRAD_SOURCE_FILES})
target_include_directories(dft_grad_un PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(dft_grad_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(dft_grad_opt OBJECT ${DFT_GRAD_OPT_SOURCE_FILES})
target_include_directories(dft_grad_opt PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(dft_grad_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(dft_grad_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
