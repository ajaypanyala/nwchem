
set(DIANA_SOURCE_FILES
./dia_initx.F
./dia_rtdbin.F
./dia_report.F
./dia_rdref.F
./dia_wtfram.F
./dia_rmsdev.F
./dia_rdf.F
./dia_angle.F
./dia_diel.F
./dia_rama.F
./dia_rdhdr.F
./dia_wrtpdb.F
./dia_super.F
./dia_define.F
./dia_input.F
./dia_gethdr.F
./dia_esppt.F
./dia_tag.F
./dia_init.F
./dia_average.F
./dia_getref.F
./dia_select.F
./dia_task.F
./dia_pltgrd.F
./diana.F
./dia_inirms.F
./dia_bond.F
./dia_alloc.F
./dia_index.F
./dia_hbond.F
./dia_wthdr.F
./dia_rdfram.F
./dia_torsion.F
./dia_gettop.F
./dia_distance.F
./dia_wrtxyz.F
#./dia_saver.F
./dia_xaver.F
./dia_chgdst.F
./dia_rdtop.F
./dia_groups.F
./dia_evaluate.F
./dia_covar.F
./dia_sizref.F
./dia_getsel.F
./dia_gettrj.F
./dia_data.F
./dia_atrad.F
./dia_finish.F
./dia_espplt.F
)

add_library(diana STATIC ${DIANA_SOURCE_FILES})
#add_dependencies(diana ${NWC_DEPENDENCIES})
target_include_directories(diana PUBLIC ${NWC_INCLUDES})
target_compile_options(diana PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(diana PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

