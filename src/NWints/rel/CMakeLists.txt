set(NWINTS_REL_OPT_SOURCE_FILES
rel_init_con.F
rel_LLSSd.F
rel_SSLLd.F
rel_oneld.F
rel_pot2.F
rel_LLSS.F
rel_onel.F
rel_SSSS.F
rel_SSLL.F
rel_mpole.F
rel_2e4c_sf.F
rel_pvp.F
zora_input.F
dmd_input.F
rel_SSSSd.F
rel_pot.F
rel_oneld_cosmo.F
rel_2e4cd_sf.F
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
