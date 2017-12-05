
set(BQ_SOURCE_FILES
./bq_main.F
./bq_input.F
./bq_force.F
./bq_checksum.F
./bq_rtdb.F
./bq_data.F
)

add_library(bq STATIC ${BQ_SOURCE_FILES})
#add_dependencies(bq ${NWC_DEPENDENCIES})
target_include_directories(bq PUBLIC ${NWC_INCLUDES})
target_compile_options(bq PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(bq PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

