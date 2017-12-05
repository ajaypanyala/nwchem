set(NWINTS_DK_SOURCE_FILES
int_1edkso_ga.F
bas_create_uncont.F
dk_input.F
int_1edk_ga.F
int_pxvp_ga.F
dkdiv_energy.F
)

add_library(nwints_dk OBJECT ${NWINTS_DK_SOURCE_FILES})
target_include_directories(nwints_dk PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_dk PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwints_dk PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})