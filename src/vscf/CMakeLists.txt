
set(VSCF_SOURCE_FILES
./vscf_input.F
./vscf.F
)

add_library(vscf STATIC ${VSCF_SOURCE_FILES})
#add_dependencies(vscf ${NWC_DEPENDENCIES})
target_include_directories(vscf PUBLIC ${NWC_INCLUDES})
target_compile_options(vscf PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(vscf PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
