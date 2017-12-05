set(NWINTS_TEXAS_SOURCE_FILES
lab_reorder.F
zeroint.F
spec_block.F
getmem.F
screened_op.F
uniq_pairs.F
texas_face.F
texas_hf.F
spec_calcint.F
memoha.F
shells.F
spec_destiny.F
precalc2.F
rt123.F
d_datlog.F
derivat.F
transfor.F
prepint.F
assemblx.F
big_service.F
trobsa.F
conbuf2_der.F
amshift.F
convert.f
cshneg.f
destdul.f
gencon.f
missing.f
num_diff.f
test_calc2.f
time.f
texas_stubs.F
)


add_library(nwints_texas_un OBJECT abcd.F)
target_compile_definitions(nwints_texas_un PUBLIC USE_TEXAS)
target_include_directories(nwints_texas_un PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_texas_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(nwints_texas_opt OBJECT ${NWINTS_TEXAS_SOURCE_FILES})
target_compile_definitions(nwints_texas_opt PUBLIC USE_TEXAS)
target_include_directories(nwints_texas_opt PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_texas_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwints_texas_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
