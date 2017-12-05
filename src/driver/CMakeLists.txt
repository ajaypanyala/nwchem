
set(DRIVER_SOURCE_FILES
./lst_drv.F
./lbfgs.F
./lbfgs-b_driver.F
./oniom.F
./lbfgs_utils.F
./opt_drv.F
./lbfgs_driver.F
./lbfgs-b.F
)


add_library(driver_opt OBJECT ./driver_input.F)
target_include_directories(driver_opt PUBLIC ${NWC_INCLUDES})
target_compile_options(driver_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(driver_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})


add_library(driver STATIC ${DRIVER_SOURCE_FILES} $<TARGET_OBJECTS:driver_opt>)
#add_dependencies(driver ${NWC_DEPENDENCIES})
target_include_directories(driver PUBLIC ${NWC_INCLUDES})
target_compile_options(driver PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
