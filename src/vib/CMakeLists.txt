
set(VIB_SOURCE_FILES
./vib_chkhess.F
./vib_nmass.F
./vib_wrtfreq.F
./vib_prdata.F
./vib_phonon.F
./vib_rdinp.F
./vib_cebdab.F
./vib_setup.F
./vib_vib.F
./vib_hmass.F
./vib_bend.F
./vib_istep.F
./vib_dtrngl.F
./vib_bdiha5.F
./vib_cross.F
./vib_wmatrix.F
./vib_zone_input.F
./vib_thermo.F
./vib_transp.F
./vib_zmatin.F
./vib_clean.F
./vib_prntpd.F
./vib_bcalc.F
./vib_dsquar.F
./vib_eckart.F
./vib_intense.F
./vib_input.F
./vib_bstr.F
./vib_tors.F
./vib_libe.F
./vib_gtors.F
)

add_library(vib_opt OBJECT vib_callg.F)
target_include_directories(vib_opt PUBLIC ${NWC_INCLUDES})
target_compile_options(vib_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(vib_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

add_library(vib STATIC ${VIB_SOURCE_FILES} $<TARGET_OBJECTS:vib_opt>)
#add_dependencies(vib ${NWC_DEPENDENCIES})
target_include_directories(vib PUBLIC ${NWC_INCLUDES})
target_compile_options(vib PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
