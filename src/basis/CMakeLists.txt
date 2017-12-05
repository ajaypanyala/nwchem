set(BASIS_INCLUDES
${PROJECT_SOURCE_DIR}/src/NWints/rel
)

set(BASIS_SOURCE_FILES
./bas_vec_info.F
./testbasis.F
./libcheck.F
./bas_contrib.F
./bas_input.F
./bas_checksum.F
./basisP.F
./basis.F
./newbasis.F
./compare.F
./bas_blas.F
./bas_blasP.F
)

add_library(basis STATIC ${BASIS_SOURCE_FILES})
#add_dependencies(basis ${NWC_DEPENDENCIES})
target_include_directories(basis PUBLIC ${NWC_INCLUDES} ${BASIS_INCLUDES})
target_compile_options(basis PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(basis PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

