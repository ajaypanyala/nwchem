
set(DIIS_SOURCE_FILES
tce_jacobi_t3a.F
tce_diis.F
tce_jacobi_l1.F
tce_jacobi_t2.F
tce_jacobi_y3.F
tce_jacobi_x3.F
tce_jacobi_l3.F
tce_jacobi_y2.F
tce_diis_tidy.F
tce_jacobi_l4.F
tce_diis_switch.F
tce_jacobi_t4.F
tce_diis_init.F
tce_jacobi_x4.F
tce_jacobi_x2.F
tce_jacobi_y4.F
tce_jacobi_t3.F
tce_jacobi_t1.F
tce_jacobi_x1.F
tce_jacobi_y1.F
tce_jacobi_l2.F
tce_jacobi_x3a.F
)

add_library(tce_diis STATIC ${DIIS_SOURCE_FILES})
#add_dependencies(tce_diis ${NWC_DEPENDENCIES})
target_include_directories(tce_diis PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_diis PUBLIC DEBUG_PRINT)

target_compile_options(tce_diis PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_diis PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
