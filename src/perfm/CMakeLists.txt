
set(PERFM_SOURCE_FILES
./perfm_end.F
./perfm_start.F
./perfm_flop.F
#./testdist.F
)

add_library(perfmC OBJECT ./dist.c)
target_include_directories(perfmC PUBLIC ${NWC_INCLUDES})
target_compile_options(perfmC PUBLIC ${NWCHEM_C_COMPILE_FLAGS})
target_compile_options(perfmC PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})

add_library(perfm STATIC ${PERFM_SOURCE_FILES} $<TARGET_OBJECTS:perfmC>)
#add_dependencies(perfm ${NWC_DEPENDENCIES})
target_include_directories(perfm PUBLIC ${NWC_INCLUDES})
target_compile_options(perfm PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(perfm PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
