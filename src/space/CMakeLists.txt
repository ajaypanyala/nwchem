
set(SPACE_SOURCE_FILES
./sp_main.F
./sp_init.F
)

add_library(space STATIC ${SPACE_SOURCE_FILES})
#add_dependencies(space ${NWC_DEPENDENCIES})
target_include_directories(space PUBLIC ${NWC_INCLUDES})
target_compile_options(space PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

