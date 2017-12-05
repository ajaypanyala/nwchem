
set(DIIS2_SOURCE_FILES
tce_jacobi_lr4.F
tce_jacobi_lr2.F
tce_jacobi_tr4.F
tce_jacobi_tr1.F
tce_diis2.F
tce_diis2_init.F
tce_jacobi_tr3a.F
tce_jacobi_lr3.F
tce_jacobi_tr2.F
tce_jacobi_tr3.F
tce_jacobi_lr1.F
tce_diis2_tidy.F
)

add_library(tce_diis2 STATIC ${DIIS2_SOURCE_FILES})
#add_dependencies(tce_diis2 ${NWC_DEPENDENCIES})
target_include_directories(tce_diis2 PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_diis2 PUBLIC DEBUG_PRINT)

target_compile_options(tce_diis2 PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_diis2 PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
