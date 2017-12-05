set(NWINTS_SIMINT_OPT_SOURCE_FILES
source/nwcsim_facef90.F 
source/nwcsim_hf2.F
)

set(NWINTS_REL_SOURCE_FILES
int_rel_init.F
rel_input.F
modelpotential_input.F
set_modelpotential_params.F
get_modelpotential_params.F
calc_modelpotential.F
)

add_library(nwints_rel_un OBJECT ${NWINTS_REL_SOURCE_FILES})
target_include_directories(nwints_rel_un PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_rel_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(nwints_rel_opt OBJECT ${NWINTS_REL_OPT_SOURCE_FILES})
target_include_directories(nwints_rel_opt PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_rel_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwints_rel_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
