
set(MOINTS_SOURCE_FILES
./moints_stats.F
./moints_1e.F
./mointdriver.F
./direct_mp2.F
./moints_full.F
./atwtest.F
./moints_mem.F
./moints_ccsd.F
./moints_dbg.F
./sifs1.F
./sifs0.F
./moints_mp2.F
./get_JKelt.F
./moints_gsh1.F
)

add_library(moints_opt OBJECT ./moints_aodisk.F ./moints_2x.F ./moints_6x.F ./moints_aux2.F ./moints_semi.F)
target_include_directories(moints_opt PUBLIC ${NWC_INCLUDES})
target_compile_options(moints_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(moints_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

add_library(mointsC OBJECT ./moints_cutil.c) #./crayio.c
target_include_directories(mointsC PUBLIC ${NWC_INCLUDES})
target_compile_options(mointsC PUBLIC ${NWCHEM_C_COMPILE_FLAGS})
target_compile_options(mointsC PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})

add_library(moints STATIC ${MOINTS_SOURCE_FILES} $<TARGET_OBJECTS:mointsC> $<TARGET_OBJECTS:moints_opt>)
#add_dependencies(moints ${NWC_DEPENDENCIES})
target_include_directories(moints PUBLIC ${NWC_INCLUDES})
target_compile_options(moints PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
