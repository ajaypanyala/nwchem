
set(CPHF_SOURCE_FILES
./cphf_solve.F
./cphf_solve2.F
./cphf_solve3.F
./cphf_int_dbg.F
./cphf_solve4.F
./cphf_poliz.F
./cphf.F
#./mjdtest.F
)

add_library(cphf STATIC ${CPHF_SOURCE_FILES})
#add_dependencies(cphf ${NWC_DEPENDENCIES})
target_include_directories(cphf PUBLIC ${NWC_INCLUDES})

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(cphf PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS} -ffixed-line-length-none)
else()
  target_compile_options(cphf PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
endif()
#target_compile_options(cphf PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
