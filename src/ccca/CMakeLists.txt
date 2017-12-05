
set(CCCA_SOURCE_FILES
./ccca_utils.F
./task_ccca.F
)

add_library(ccca STATIC ${CCCA_SOURCE_FILES})
#add_dependencies(ccca ${NWC_DEPENDENCIES})
target_include_directories(ccca PUBLIC ${NWC_INCLUDES})
target_compile_options(ccca PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})


