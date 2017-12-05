set(NWC_DEVELOP_INCLUDES
${PROJECT_SOURCE_DIR}/src/NWints/rel
${PROJECT_SOURCE_DIR}/src/NWints/ecp
)

set(DEVELOP_SOURCE_FILES
./rak25.F
./rak26.F
./rak18.F
./rak30.F
./rak24.F
./rak27.F
./rjhtest.F
./raktest.F
./rak20.F
./rak22.F
./rak21.F
./rak29.F
./cct_UHF.F
./raktest_semi.F
./rak19.F
./gaussq_netlib.F
./rak28.F
./rak31.F
./debtest.F
./rjhtrans.F
./rak23.F
./kgdtest.F
./jvltest.F
)


add_library(develop_opt OBJECT ./jantest.F ./uccsdtest.F)
target_include_directories(develop_opt PUBLIC ${NWC_INCLUDES} ${NWC_DEVELOP_INCLUDES})
target_compile_options(develop_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(develop_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

add_library(develop STATIC ${DEVELOP_SOURCE_FILES} $<TARGET_OBJECTS:develop_opt>)
#add_dependencies(develop ${NWC_DEPENDENCIES})
target_include_directories(develop PUBLIC ${NWC_INCLUDES} ${NWC_DEVELOP_INCLUDES})
target_compile_options(develop PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})


