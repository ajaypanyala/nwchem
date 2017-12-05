
set(LEPS_SOURCE_FILES
crossing.F
gasdev.F
mc_init.F
mc_main.F
md_data.F
md_driver.F
md_set.F
mymd_input.F
pred-corr.F
urand.F
)

add_library(leps_opt OBJECT leps_energy.F leps_gradient.F mc_driver.F mymc_input.F mc_data.F)
target_include_directories(leps_opt PUBLIC ${NWC_INCLUDES})
target_compile_options(leps_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(leps_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

add_library(leps STATIC ${LEPS_SOURCE_FILES} $<TARGET_OBJECTS:leps_opt>)
#add_dependencies(leps ${NWC_DEPENDENCIES})
target_include_directories(leps PUBLIC ${NWC_INCLUDES})
target_compile_options(leps PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
