
set(DNTMC_SOURCE_FILES
./dntmc_setgeom.F
./dntmc_input.F
)

add_library(dntmc_opt OBJECT ./gibbs.F)
target_include_directories(dntmc_opt PUBLIC ${NWC_INCLUDES})
target_compile_options(dntmc_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(dntmc_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})


add_library(dntmc STATIC ${DNTMC_SOURCE_FILES} $<TARGET_OBJECTS:dntmc_opt>)
#add_dependencies(dntmc ${NWC_DEPENDENCIES})
target_include_directories(dntmc PUBLIC ${NWC_INCLUDES})
target_compile_options(dntmc PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
