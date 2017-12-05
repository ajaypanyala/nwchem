
set(CISDTQ_SOURCE_FILES
cisdtq_c3.F
cisdtq_c4.F
cisdtq_e.F
cisdtq_c2.F
cisdtq_c1.F
)

add_library(tce_cisdtq STATIC ${CISDTQ_SOURCE_FILES})
#add_dependencies(tce_cisdtq ${NWC_DEPENDENCIES})
target_include_directories(tce_cisdtq PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
if(CCSDTQ)
  target_compile_definitions(tce_cisdtq PUBLIC DEBUG_PRINT CCSDTQ)
else()
  target_compile_definitions(tce_cisdtq PUBLIC DEBUG_PRINT)
endif()
target_compile_options(tce_cisdtq PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_cisdtq PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
