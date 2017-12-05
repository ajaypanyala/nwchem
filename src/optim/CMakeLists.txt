
set(OPTIM_SOURCE_FILES
./string/string_input.F
./string/string.F
./mepgs/mepgs_input.F
./tropt/tropt_input.F
./deloc/deloc.F
./deloc/deloc_driver.F
./deloc/deloc_utils.F
./neb/neb_input.F
./neb/neb_drv.F
./neb/neb_utils.F
./neb/bead_list.F
)


add_library(optim_un OBJECT ./mepgs/mepgs_drv.F ./tropt/tropt_drv.F)
target_include_directories(optim_un PUBLIC ${NWC_INCLUDES})
target_compile_options(optim_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(optim STATIC ${OPTIM_SOURCE_FILES} $<TARGET_OBJECTS:optim_un>)
#add_dependencies(optim ${NWC_DEPENDENCIES})
target_include_directories(optim PUBLIC ${NWC_INCLUDES})
target_compile_options(optim PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(optim PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
