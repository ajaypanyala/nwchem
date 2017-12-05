
set(EMB_SOURCE_FILES
c1v.F
emb1.F
emb2.F
c2vt1_s.F
)

add_library(tce_emb STATIC ${EMB_SOURCE_FILES})
#add_dependencies(tce_emb ${NWC_DEPENDENCIES})
target_include_directories(tce_emb PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_emb PUBLIC DEBUG_PRINT)
target_compile_options(tce_emb PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_emb PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
