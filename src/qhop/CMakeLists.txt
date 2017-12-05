
set(QHOP_SOURCE_FILES
./qhop_asite.F
./qhop_pindex.F
./qhop_setup_good.F
#./qhop_newpc.F
./qhop_setup_tj.F
./qhop_init.F
./qhop_dsite.F
./qhop_setup.F
./qhop_final.F
./qhop_prob.F
)

add_library(qhop STATIC ${QHOP_SOURCE_FILES})
#add_dependencies(qhop ${NWC_DEPENDENCIES})
target_include_directories(qhop PUBLIC ${NWC_INCLUDES})
target_compile_options(qhop PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(qhop PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
