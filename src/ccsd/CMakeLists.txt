
set(CCSD_SOURCE_FILES
./ccsd_itm_omp.F
./ccsd.F
./ccsd_iterdrv2.F
#./ccsd_tengy_bgp2.F
./ccsd_tengy.F
./ccsd_input.F
./qstat.F
./ccsd_tzinta.F
./aoccsd2.F
./ccsd_itm.F
./ccsd_z2pm.F
./ccsd_t2pm.F
./ccsd_get_param.F
./rddcom.F
./ccsd_trpdrv.F
./moints_trp.F
./convert_single_double.F
./ccden_driver.F
./ccsd_driver.F
./cputime.F
./ccsd_idx2.F
./ccsd_tzintb.F
./ccsd_aosig.F
./ccsd_pzamp.F
./ccsd_trpdrv_nb.F
./ccden_interm4.F
./ccsd_pdiis.F
./ccsd_idx34.F
./ccsd_idx1.F
./ccsd_mem.F
./ccsd_graddrv.F
./ccden_interm2.F
./ccsd_pstat.F
./ccsd_pampt3.F
./ccden_1pdm.F
./ccden_2pdma.F
./ccsd_trpmo.F
./ccsd_zsig2.F
./geaxpy.F
./ccsd_rdtrpo.F
./ccsd_mktau_ga.F
./ccsd_trpdrv_bgp2.F
./ccsd_fsig2.F
./ccsd_aux.F
./ccsd_zvecs.F
./ccsd_tengy_omp.F
./ccden_2pdmb.F
./ccsd_trpdrv_omp_mixprec.F
./ccsd_mktask.F
./ccsd_trpdrv_omp.F
./ccsd_pampt2.F
./ccsd_mkt2.F
./ccsd_zsig1.F
./ccsd_tengy_bgp.F
./ccsd_sortmo.F
./ccsd_fsig1.F
)

add_library(ccsd STATIC ${CCSD_SOURCE_FILES})
#add_dependencies(ccsd ${NWC_DEPENDENCIES})
target_include_directories(ccsd PUBLIC ${NWC_INCLUDES})
target_compile_options(ccsd PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(ccsd PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "XL")
  target_compile_definitions(ccsd PUBLIC TRPMIXP_OFF)
endif()