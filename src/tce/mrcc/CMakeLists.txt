add_definitions(-DMRCC_METHODS)

add_subdirectory(ccsd)
add_subdirectory(uss)
add_subdirectory(ccsd_sub)
add_subdirectory(ccsd_t)

set(TCE_MRCC_INCLUDES
${PROJECT_SOURCE_DIR}/src/NWints/rel
)

set(MRCC_SOURCE_FILES
tce_mrcc_diis_new.F
pgcreatefile.F
tce_mrcc_refs_sym.F
create_mr_evl_sorted.F
tce_mrcc_diis_init.F
tce_mrcc_uhf.F
tce_mrcc_mo2e.F
tce_mrcc_diis.F
atpart_nalength_mrcc.F
tce_mrcc_macpy.F
tce_mrcc_c1_offset.F
tce_t4ac_offset.F
tce_mrcc_tile.F
tce_mrcc_c2.F
tce_jacobi_sub_t2.F
tce_jacobi_sub_t1.F
tce_mrcc_calcexcit.F
tce_mrcc_create_cas.F
tce_mrcc_energy.F
tce_mrcc_c1.F
tce_mrcc_diis_tidy.F
tce_mrcc_jacobi_t1.F
tce_mrcc_c2_offset.F
tce_mrcc_create_cas1.F
tce_mrcc_print.F
eaf_mrcc.F
tce_mrcc_ao2e.F
tce_mrcc_gacp.F
tce_mrcc_mul_1.F
tce_mrcc_fillindexarray.F
tce_mrcc_diis_save.F
tce_mrcc_heff.F
tce_mrcc_readref.F
tce_t3ac_offset.F
tce_mrcc_jacobi_t2.F
tce_subnxtask.F
tce_mrcc_daxfile.F
new_ga4ind_N5uhf_mrcc_lm2.F
tce_mrcc_force_orben.F
tce_mrcc_mul_2.F
tce_mrcc_sortint.F
)

add_library(tcemrcc_un OBJECT tce_mrcc_input.F)
target_include_directories(tcemrcc_un PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES} ${TCE_MRCC_INCLUDES})
target_compile_options(tcemrcc_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(tcemrccC OBJECT tce_mrcc_sublock.c)
target_include_directories(tcemrccC PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES} ${TCE_MRCC_INCLUDES})
target_compile_options(tcemrccC PUBLIC ${NWCHEM_C_COMPILE_FLAGS})
target_compile_options(tcemrccC PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})

add_library(tce_mrcc STATIC ${MRCC_SOURCE_FILES} $<TARGET_OBJECTS:tcemrccC> $<TARGET_OBJECTS:tcemrcc_un>)
#add_dependencies(tce_mrcc ${NWC_DEPENDENCIES})
target_include_directories(tce_mrcc PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES} ${TCE_MRCC_INCLUDES})
if(MRCC_METHODS)
  target_compile_definitions(tcemrccC PUBLIC MRCC_METHODS)
  target_compile_definitions(tcemrcc_un PUBLIC MRCC_METHODS)
  target_compile_definitions(tce_mrcc PUBLIC MRCC_METHODS)
endif()
target_compile_options(tce_mrcc PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_mrcc PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

