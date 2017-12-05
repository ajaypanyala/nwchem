set(NWINTS_HONDO_SOURCE_FILES
hnd_sxyz.F
hnd_tvd_ij.F
hnd_giaol10.F
hnd_root5.F
hnd_rt123.F
hnd_wleg.F
hnd_dso.F
hnd_pso.F
hnd_dipxyz.F
hnd_giaoh11.F
hnd_deriaj.F
hnd_dervxyz.F
hnd_initdd.F
hnd_init_dip.F
hnd_init_dso.F
hnd_init_1emtpl.F
hnd_giaspd.F
hnd_giaotv10.F
hnd_veloc.F
hnd_ddtxyz.F
hnd_giaos100.F
hnd_oechrg.F
hnd_hlfd_ij.F
hnd_1emem.F
hnd_dd1emem.F
hnd_spdgia.F
hnd_initd.F
hnd_txyz.F
hnd_giaobq10.F
hnd_giaos10.F
hnd_jkwrys.F
hnd_stvintd.F
hnd_init_1eelec.F
hnd_stvintdd.F
hnd_giaoh01.F
hnd_init.F
hnd_dipvel.F
hnd_dderij.F
hnd_giaxyz.F
hnd_terminate.F
hnd_giao_init.F
hnd_giagnm.F
hnd_giao_1emem.F
hnd_elfder.F
hnd_giahnd.F
hnd_init_pso.F
hnd_angmom.F
hnd_giao_2emem.F
hnd_jkbcdf.F
hnd_root4.F
hnd_stvint.F
hnd_d1emem.F
hnd_set_acc.F
hnd_denxyz.F
hnd_xyzgia.F
hnd_herxyz.F
hnd_sd_ij.F
hnd_whermt.F
hnd_hlf_ij.F
hnd_droot.F
#hnd_dipder.F
)


add_library(nwints_hondo_un OBJECT hnd_vstat.F)
target_compile_definitions(nwints_hondo_un PUBLIC DBLEPREC_RYS)
target_include_directories(nwints_hondo_un PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_hondo_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(nwints_hondo_opt OBJECT ${NWINTS_HONDO_SOURCE_FILES})
target_compile_definitions(nwints_hondo_opt PUBLIC DBLEPREC_RYS)
target_include_directories(nwints_hondo_opt PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_hondo_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwints_hondo_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
