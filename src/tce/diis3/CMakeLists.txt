
set(DIIS3_SOURCE_FILES
tce_jacobi2_ir1.F
tce_jacobi_ir4.F
tce_jacobi_ir2.F
tce_jacobi_ir3.F
tce_jacobi2_ir4.F
tce_jacobi2_ir2.F
tce_diis3.F
tce_jacobi_ir1.F
tce_diis3c.F
tce_jacobi2_ir3.F
)

add_library(tce_diis3 STATIC ${DIIS3_SOURCE_FILES})
#add_dependencies(tce_diis3 ${NWC_DEPENDENCIES})
target_include_directories(tce_diis3 PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_diis3 PUBLIC DEBUG_PRINT)

target_compile_options(tce_diis3 PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_diis3 PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
