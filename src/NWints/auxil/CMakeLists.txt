set(NWINTS_AUXIL_OPT_SOURCE_FILES
nmcoef.F
nmcoef2c.F
#xlm_norm.org.F
int_chk_sum.F
int_norm_2c.F
int_spcart.F
spherical.F
int_norm.F
)

set(NWINTS_AUXIL_SOURCE_FILES
intb_chk_sum.F
do_chk_sum.F
xlmcoeff.F
xlm_norm.F
int_spcart2.F
)

add_library(nwints_auxil_un OBJECT ${NWINTS_AUXIL_SOURCE_FILES})
target_include_directories(nwints_auxil_un PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_auxil_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(nwints_auxil_opt OBJECT ${NWINTS_AUXIL_OPT_SOURCE_FILES})
target_include_directories(nwints_auxil_opt PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_auxil_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwints_auxil_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
