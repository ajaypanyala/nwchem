
set(ANALYZ_SOURCE_FILES
./ana_wrtpdb.F
./ana_rdref.F
./ana_bond.F
./ana_xaver.F
./ana_hbond.F
./ana_esppt.F
./ana_groups.F
./ana_rdfram.F
./ana_pltgrd.F
./ana_tag.F
./ana_rmsdev.F
./ana_index.F
./ana_super.F
./ana_input.F
./ana_rdtop.F
./ana_rdhdr.F
./ana_chgdst.F
./ana_task.F
./ana_wthdr.F
./ana_initx.F
./ana_diel.F
./ana_sizref.F
./ana_finish.F
./analyz.F
./ana_inirms.F
./ana_torsion.F
./ana_data.F
./ana_rama.F
./ana_wtfram.F
./ana_init.F
./super.F
./ana_angle.F
./ana_covar.F
./ana_rtdbin.F
./ana_select.F
./ana_rdf.F
./ana_wrtxyz.F
./ana_define.F
./ana_report.F
./ana_atrad.F
)

add_library(analyz STATIC ${ANALYZ_SOURCE_FILES})
#add_dependencies(analyz ${NWC_DEPENDENCIES})
target_include_directories(analyz PUBLIC ${NWC_INCLUDES})

target_compile_options(analyz PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(analyz PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

