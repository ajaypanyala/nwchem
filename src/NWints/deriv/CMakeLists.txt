set(NWINTS_DERIV_SOURCE_FILES
hf1PRtoAB.f
hf2PRtoAB.f
hf2QStoCD.f
hfd2oi.f
hfdnai.f
hfefi.f
hf2d_2b.F
hf2d_2a.F
hf1d.F
hf2d_3a.F
hfdkei.F
hf2d_2a2b.F
hfderi.F
hfd3oi.F
hf2d.F
hf1d_cosmo.F
hf2d_all.F
hf2d_3a2b.F
)

add_library(nwints_deriv OBJECT ${NWINTS_DERIV_SOURCE_FILES})
target_include_directories(nwints_deriv PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_deriv PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwints_deriv PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})