
set(INP_SOURCE_FILES
#./test_ilist.F
#./test.F
./inp.F
./inp_irange.F
./inp_ilist.F
)

add_library(nwc_un OBJECT ./hnd_rdfree.F)
target_include_directories(nwc_un PUBLIC ${NWC_INCLUDES})
target_compile_options(nwc_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(nwc_cinp OBJECT ./c_inp.c)
target_include_directories(nwc_cinp PUBLIC ${NWC_INCLUDES})
target_compile_options(nwc_cinp PUBLIC ${NWCHEM_C_COMPILE_FLAGS} ${NWCHEM_C_DEBUG_FLAGS})
#target_compile_options(nwc_cinp PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})

add_library(nwc_inp STATIC ${INP_SOURCE_FILES} $<TARGET_OBJECTS:nwc_cinp> $<TARGET_OBJECTS:nwc_un>)
target_include_directories(nwc_inp PUBLIC ${NWC_INCLUDES})
target_compile_options(nwc_inp PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwc_inp PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
