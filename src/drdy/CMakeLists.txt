
set(DRDY_SOURCE_FILES
./pot1.F
./setpth.F
./potsp.F
./react.F
./drdy_synch.F
./saddle.F
./io.F
./rgeom.F
./pythag.F
./projct.F
./Input_for_drdy.F
./drdy_center.F
./react2.F
./drdy.F
./rstat.F
./rotf.F
./cfloat.F
./restor.F
./rgener.F
./read_drdy_rtdb.F
./drdy_nwc_get.F
./nwc_gen.F
./moiten.F
./fdiag.F
./rword.F
./diagf.F
./tranf.F
./trans.F
./nwc_s2.F
./potset_nwchem.F
./mxlneq.F
./rphwrt.F
./drdy_sort.F
./rpath.F
./atsym.F
./euler.F
./rotg.F
./pmcube.F
./sumout.F
./rgsad.F
./pmlqsp.F
./getstr.F
./pmclqa.F
./pmlqa.F
./dsdt.F
./icint.F
./pot2.F
./rotx.F
./rsp.F
./nwc_s1.F
./drdy_input.F
#./rline.F
./path.F
./rotmat.F
./wfiles.F
./ratoms.F
./path2.F
)

add_library(drdy STATIC ${DRDY_SOURCE_FILES})
#add_dependencies(drdy ${NWC_DEPENDENCIES})
target_include_directories(drdy PUBLIC ${NWC_INCLUDES})
target_compile_options(drdy PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(drdy PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

