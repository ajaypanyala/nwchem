

set(CAFE_SOURCE_FILES
./cf_list.F
./cf_grid.F
./cf_forcep.F
./pme_grid.F
./pme_init.F
./cf_rdtop.F
./cf_d.F
./cf_uhop.F
./cf_shake.F
./cf_finish.F
./pme_fft.F
./cf_forces.F
./cf_main.F
./cf_rdgeom.F
./cf_hopping.F
./cf_init.F
./cf_f.F
./cf_induce.F
./cf_pairs.F
./cf_fp.F
./cf_pbc.F
)

add_library(cafe STATIC ${CAFE_SOURCE_FILES})
#add_dependencies(cafe ${NWC_DEPENDENCIES})
target_include_directories(cafe PUBLIC ${NWC_INCLUDES})
if (${CMAKE_HOST_SYSTEM_NAME} STREQUAL Linux OR ${CMAKE_HOST_SYSTEM_NAME} STREQUAL Darwin)
    target_compile_definitions(cafe PUBLIC NEED_LOC)
endif ()
target_compile_options(cafe PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(cafe PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

