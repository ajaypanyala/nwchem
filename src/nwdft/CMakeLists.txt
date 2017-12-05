set(NWDFT_INCLUDES
${PROJECT_SOURCE_DIR}/src/perfm
${PROJECT_SOURCE_DIR}/src/NWints/rel
${PROJECT_SOURCE_DIR}/src/NWints/int
${PROJECT_SOURCE_DIR}/src/nwdft/input_dft
${PROJECT_SOURCE_DIR}/src/nwdft/include
${PROJECT_SOURCE_DIR}/src/nwdft/spec
${PROJECT_SOURCE_DIR}/src/nwdft/grid
${PROJECT_SOURCE_DIR}/src/nwdft/util
${PROJECT_SOURCE_DIR}/src/nwdft/xc
${PROJECT_SOURCE_DIR}/src/nwdft/scf_dft
)

set(NWDFT_SOURCE_FILES
./nwdft.F
./input_dft/cdft_inp.F
./input_dft/conv_inp.F
./input_dft/occup_input.F
./input_dft/grid_inp.F
./input_dft/tol_inp.F
./input_dft/dft_input.F
./input_dft/dft_rdinput.F
./input_dft/xc_inp.F
./input_dft/dft_inpana.F
)

add_subdirectory(xc)
add_subdirectory(grid)
add_subdirectory(zora)
add_subdirectory(util)
add_subdirectory(spec)
add_subdirectory(coulomb)
add_subdirectory(dftgrad)
add_subdirectory(so_dft)
add_subdirectory(scf_dft)
add_subdirectory(scf_dft_cg)
add_subdirectory(rt_tddft)
add_subdirectory(lr_tddft)
add_subdirectory(lr_tddft_grad)

add_library(nwdft STATIC ${NWDFT_SOURCE_FILES} $<TARGET_OBJECTS:nwdft_xc> $<TARGET_OBJECTS:grid_un>
 $<TARGET_OBJECTS:nwdft_xc_un> $<TARGET_OBJECTS:nwdft_xc_c> $<TARGET_OBJECTS:rt_tddft>
 $<TARGET_OBJECTS:grid_opt> $<TARGET_OBJECTS:lr_tddft> $<TARGET_OBJECTS:lr_tddft_grad>
 $<TARGET_OBJECTS:nwdft_zora> $<TARGET_OBJECTS:so_dft> $<TARGET_OBJECTS:dft_util> 
 $<TARGET_OBJECTS:dft_util_un> $<TARGET_OBJECTS:scf_dft_opt> $<TARGET_OBJECTS:scf_dft_un>
 $<TARGET_OBJECTS:dft_grad_un> $<TARGET_OBJECTS:dft_grad_opt> $<TARGET_OBJECTS:nwdft_spec>
 $<TARGET_OBJECTS:nwdft_coulomb> $<TARGET_OBJECTS:scf_dft_cg>)

target_include_directories(nwdft PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(nwdft PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})


