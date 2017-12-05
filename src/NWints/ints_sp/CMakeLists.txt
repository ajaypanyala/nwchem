set(NWINTS_SP_SOURCE_FILES
tabgen.f
sinfo.F
qout70.F
sp0101.F
init70.F
sp0011.F
genr70.F
sp0111.F
sp0000.F
sp0001.F
)

add_library(nwints_sp_un OBJECT sp1111.F)
target_include_directories(nwints_sp_un PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_sp_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(nwints_sp OBJECT ${NWINTS_SP_SOURCE_FILES})
target_include_directories(nwints_sp PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_sp PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwints_sp PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
