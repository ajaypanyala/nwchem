
set(PSTAT_SOURCE_FILES
./pstat_term.F
./pstat_pr_det.F
./pstat_off.F
./pstat_free.F
./pstat_alloc.F
./pstat_pr_han.F
./pstat_init.F
./pstat_pr_all.F
./pstat_on.F
)

add_library(pstat STATIC ${PSTAT_SOURCE_FILES})
#add_dependencies(pstat ${NWC_DEPENDENCIES})
target_include_directories(pstat PUBLIC ${NWC_INCLUDES})
target_compile_options(pstat PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(pstat PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
