
set(STEPPER_SOURCE_FILES
./stpr_clnstp.F
./stpr_filenames.F
./stepper_inp.F
./stpr_sumstc.F
./stpr_walk.F
./stpr_mkstep.F
./stpr_gh_at.F
./stpr_wrtsq.F
./stpr_readsq.F
./stpr_hesupd.F
./stpr_stepcor.F
./stpr_face.F
./stpr_output.F
./stpr_chkhess.F
./stpr_fndmde.F
./stpr_gen_aux.F
./stpr_save.F
./stpr_place.F
#./read_genr.F
./stpr_wlkchk.F
./stpr_formup.F
./stpr_driver.F
./distang.F
./stpr_inihes.F
./stpr_cheker.F
./stpr_cneed.f
./stpr_cross.f
./stpr_datestp.f
./stpr_gmass.f
./stpr_hmass.f
./stpr_locate.f
./stpr_matout.f
./stpr_mgeom.f
./stpr_mktri.f
./stpr_pmat.f
./stpr_prntpd.f
./stpr_projec.f
./stpr_putbck.f
./stpr_setmass.f
./stpr_stepper.f
./stpr_tprntd.f
./stpr_uprcase.f
)


add_library(stepper_opt OBJECT ./stpr_partit.F)
target_include_directories(stepper_opt PUBLIC ${NWC_INCLUDES})
target_compile_options(stepper_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(stepper_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})


add_library(stepper STATIC ${STEPPER_SOURCE_FILES} $<TARGET_OBJECTS:stepper_opt>)
#add_dependencies(stepper ${NWC_DEPENDENCIES})
target_include_directories(stepper PUBLIC ${NWC_INCLUDES})
target_compile_options(stepper PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})


