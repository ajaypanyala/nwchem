
set(RIMP2_SOURCE_FILES
./iamax.f
./lnblnk.f  
#./panel_diag.f
./pr_corr_map.F
./close_3ci.F
./open_3ci.F
./max_v_diag.F
./put_3cints.F
./rimp2_xf.F
./get_3cints.F
./xf3ci_step3.F
./rimp2_parm.F
./dgecop.F
./rimp2_timers.F
./ga_by_shells.F
./mk_v.F
./rimp2_v_e2.F
./occ_ranges.F
./xf3ci_mem12.F
./xf3ci_step12.F
./rimp2.F
./mk_1d_shdist.F
./rimp2_e_mem.F
./prog_report.F
./xform_3cints.F
)


add_library(rimp2_opt OBJECT rimp2_denom.F sum_pairs.F xf3ci_s12a.F xf3ci_s12b.F eval_3ci.F)
target_include_directories(rimp2_opt PUBLIC ${NWC_INCLUDES})
target_compile_options(rimp2_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(rimp2_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})


add_library(rimp2 STATIC ${RIMP2_SOURCE_FILES} $<TARGET_OBJECTS:rimp2_opt>)
#add_dependencies(rimp2 ${NWC_DEPENDENCIES})
target_include_directories(rimp2 PUBLIC ${NWC_INCLUDES})
target_compile_options(rimp2 PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})


