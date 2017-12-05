
set(GRADIENTS_INCLUDES
${PROJECT_SOURCE_DIR}/src/nwdft/xc
${PROJECT_SOURCE_DIR}/src/NWints/int
)

set(GRADIENTS_SOURCE_FILES
./gradients.F
#./grad_store.F
./grad_force.F
./scf_gradient.F
./grad_inp.F
./grad1.F
./ga_reorder.F
./grad_dens.F
)

add_library(gradients_opt OBJECT ./grad2.F ./grad_getdens.F)
target_include_directories(gradients_opt PUBLIC ${NWC_INCLUDES} ${GRADIENTS_INCLUDES})
target_compile_options(gradients_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(gradients_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

add_library(gradients STATIC ${GRADIENTS_SOURCE_FILES} $<TARGET_OBJECTS:gradients_opt>)
#add_dependencies(gradients ${NWC_DEPENDENCIES})
target_include_directories(gradients PUBLIC ${NWC_INCLUDES} ${GRADIENTS_INCLUDES})
target_compile_options(gradients PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
