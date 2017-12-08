set(NWINTS_SIMINT_OPT_SOURCE_FILES
source/nwcsim_facef90.F 
source/nwcsim_hf2.F
module/nwcsim_module.F
)

add_library(nwints_simint OBJECT ${NWINTS_SIMINT_OPT_SOURCE_FILES})
target_include_directories(nwints_simint PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_simint PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwints_simint PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
