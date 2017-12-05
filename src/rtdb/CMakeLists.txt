
set(RTDB_SOURCE_FILES
./test.F
./davetest.F
./testgr.F
./cntx.F
)

set(RTDB_C_FILES
./context.c
./context_f2c.c
./rtdb_f2c.c
./rtdb.c
./rtdb_seq.c
./hdbm/hdbm.c
)

add_library(rtdbC OBJECT ${RTDB_C_FILES})
target_include_directories(rtdbC PUBLIC ${NWC_INCLUDES})
target_compile_definitions(rtdbC PUBLIC USE_HDBM)
target_compile_options(rtdbC PUBLIC ${NWCHEM_C_COMPILE_FLAGS})
target_compile_options(rtdbC PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})

add_library(rtdb STATIC ${RTDB_SOURCE_FILES} $<TARGET_OBJECTS:rtdbC>)
#add_dependencies(rtdb ${NWC_DEPENDENCIES})
target_include_directories(rtdb PUBLIC ${NWC_INCLUDES})
target_compile_definitions(rtdb PUBLIC USE_HDBM)

target_compile_options(rtdb PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(rtdb PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
