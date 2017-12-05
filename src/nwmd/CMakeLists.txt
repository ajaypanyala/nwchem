
set(NWMD_SOURCE_FILES
./md_input.F
./prp_init.F
./prp_main.F
./nwmd.F
./md_start.F
./md_main.F
./md_finish.F
)

add_library(nwmd STATIC ${NWMD_SOURCE_FILES})
#add_dependencies(nwmd ${NWC_DEPENDENCIES})
target_include_directories(nwmd PUBLIC ${NWC_INCLUDES})
if (${CMAKE_HOST_SYSTEM_NAME} STREQUAL Linux OR ${CMAKE_HOST_SYSTEM_NAME} STREQUAL Darwin)
    target_compile_definitions(nwmd PUBLIC NEED_LOC)
endif ()
target_compile_options(nwmd PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwmd PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
