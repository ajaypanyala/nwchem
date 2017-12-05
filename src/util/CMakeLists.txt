set(UTIL_INCLUDES
${PROJECT_SOURCE_DIR}/src/nwdft/xc
${PROJECT_SOURCE_DIR}/src/hessian/analytic
)


set(UTIL_OPT_SOURCE_FILES
ga_matpow.F 
util_pack.F 
dabssum.F 
dabsmax.F 
dfill.F 
ifill.f  
mabyte_fill.F 
zfill.f 
dgami.f                 
ga_it_lsolve.F 
ga_it_orth.F 
ga_orthog.F 
idamin.f 
ga_copy_dz.F 
util_eri_cholesky.F
util_jacobi.F 
util_rsg.f 
stpr_sjacobi.F  
nxtask.F 
util_mirror.F 
util_sillysort.F 
util_test_cholesky.F
)
                 

set(UTIL_SOURCE_FILES
./ao_1prdm_read.F
./ao_1prdm_write.F
./util_getenv.F
./ga_iter_diag.F
./util_intrsc2008.F
./ga_local_mdot.F
./indint.F
./ga_lkain_2cpl3_ext.F
./ma_print.F
./dgewr.F
./ga_it_proj.F
./ga_mix.F
./util_erf.F
./ga_maxelt.F
./util_ascii_plot.F
./input_echo.F
./nwmolden.F
./corr_mk_ref.F
./util_getarg.F
./util_xyz.F
./ga_lkain_2cpl3.F
./ga_mat2col.F
./bgj.F
./util_print.F
./atoi.F
#./error.F
./ga_lkain_ext.F
./util_nwchemrc.F
./ga_pcg_min.F
./testecce.F
./seq_output.F
./two_ind_trn.F
./ga_extra.F
./util_rtdb_state.F
./ga_normf.F
./swrite.F
./mk_fit_xf.F
./util_flush.F
./ga_list.F
./util_nwchem_srcdir.F
./util_pname.F
./ga_antisymmetr.F
./linux_shift.F
./get_density.F
./util_sgroup.F
./util_test_eig.F
./util_test.F
./ga_screen.F
./util_dra_test.F
./util_perf_test.F
./ga_orth_vec.F
./util_dgeev.F
./int_2c_ga.F
./ga_complex_diag.F
./util_nwsort.F
./ga_get_diag.F
./ga_it2.F
./movecs_ecce.F
./util_file_name.F
./sread.F
./print_center.F
./util_eaf_test.F
./util_patch_test.F
./util_rtdb_speak.F
./line_search.F
./util_ch_brd.F
./util_nwchem_paper.F
./util_io_unit.F
./util_ma.F
./util_v_stub.F
./util_transpose.F
./util_ga_test.F
./nw_inp_from_file.F
./util_cpusec.F
./freeze_input.F
#./util_labelconv.F
./ga_ran_fill.F
./util_md.F
./util_units.F
./util_legal.F
./util_allocga.F
./util_io.F
./util_job.F
./ga_conjg.F
./testsolve.F
./ga_trace_dg.F
./errquit.F
./util_ndim_test.F
./util_sf_test.F
./moeig_read.F
./vtune.F
./util_random.F
./ga_zgop.F
./util_test_lu.F
./ffflush.F
./util_wallsec.F
./nga_idamax.F
./util_length_scale.F
./nxtval_ga.F
#./util_module_avail.F
#./md5/ftest.F
banner.f
dgefa.f
dsum.f
gather.f
icopy.f
lcopy.f
output.f
scatter.f
zsum.f
)

set(UTIL_C_OPT_SOURCES 
util_memcpy.c 
ga_accback.c 
ga_asymmetr.c 
util_gnxtval.c 
)

set(UTIL_C_SOURCES
ecce_print.c
erfc.c
fortchar.c
ga_accback.c
#ibm_cputime.c
#ieeetrap.c
linux_cpu.c
#linux_gnu_trap.c
#linux_printaff.c
linux_random.c
linux_setfpucw.c
#macx_trapfpe.c
#sgi_flush6.c
#speechclient.c
#speechserver.c
#udp.c
util.c
util_cflush.c
#utilc_nwchem_srcdir.c
#util_cuda_support.c
util_date.c
util_debug.c
util_fadvise.c
util_file_copy.c
util_file_unlink.c
util_getppn.c
util_host.c
util_md_c.c
util_md_sockets.c
util_memcpy.c
#util_mic_support.c
#util_mingw.c
util_mpinap.c
util_mygabcast.c
util_sleep.c
util_speak.c
util_system.c
util_wall_remain.c
#vtune_bind.c
#win32_cpu.c
md5/md5.c
md5/md5wrap.c
#md5/md5wrap_dox.c
#md5/mddriver.c
)


add_library(nwcutil_copt OBJECT ${UTIL_C_OPT_SOURCES})
add_library(nwcutil_fopt OBJECT ${UTIL_OPT_SOURCE_FILES})
target_include_directories(nwcutil_copt PUBLIC ${NWC_INCLUDES} ${UTIL_INCLUDES})
target_include_directories(nwcutil_fopt PUBLIC ${NWC_INCLUDES} ${UTIL_INCLUDES})

target_compile_options(nwcutil_copt PUBLIC ${NWCHEM_C_COMPILE_FLAGS})
target_compile_options(nwcutil_copt PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})
target_compile_options(nwcutil_fopt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} -ffixed-line-length-none)
target_compile_options(nwcutil_fopt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})


if (${CMAKE_HOST_SYSTEM_NAME} STREQUAL Linux)
    set(UTIL_C_SOURCES ${UTIL_C_SOURCES} ./linux_printaff.c)
endif ()

add_library(nwcutilC OBJECT ${UTIL_C_SOURCES})
target_include_directories(nwcutilC PUBLIC ${NWC_INCLUDES} ${UTIL_INCLUDES})
target_compile_options(nwcutilC PUBLIC ${NWCHEM_C_COMPILE_FLAGS} ${NWCHEM_C_DEBUG_FLAGS})

add_library(nwcutil STATIC ${UTIL_SOURCE_FILES} $<TARGET_OBJECTS:nwcutilC> $<TARGET_OBJECTS:nwcutil_copt> $<TARGET_OBJECTS:nwcutil_fopt>)
target_include_directories(nwcutil PUBLIC ${NWC_INCLUDES} ${UTIL_INCLUDES})

if (${CMAKE_HOST_SYSTEM_NAME} STREQUAL Linux OR ${CMAKE_HOST_SYSTEM_NAME} STREQUAL Darwin)
    target_compile_definitions(nwcutil PUBLIC NEED_LOC BYTE_SHORT_PACKING NWCHEM_SRCDIR="${PROJECT_SOURCE_DIR}/src")
    target_compile_definitions(nwcutilC PUBLIC NEED_LOC BYTE_SHORT_PACKING NWCHEM_SRCDIR="${PROJECT_SOURCE_DIR}/src")
    target_compile_definitions(nwcutil_copt PUBLIC NEED_LOC BYTE_SHORT_PACKING NWCHEM_SRCDIR="${PROJECT_SOURCE_DIR}/src")
    target_compile_definitions(nwcutil_fopt PUBLIC NEED_LOC BYTE_SHORT_PACKING NWCHEM_SRCDIR="${PROJECT_SOURCE_DIR}/src")
endif ()

target_compile_options(nwcutil PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS} -ffixed-line-length-none)



