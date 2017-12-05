
set(GENINTERFACE_SOURCE_FILES
./geninterface.F
./NWChemWrap.F
)

add_library(geninterface STATIC ${GENINTERFACE_SOURCE_FILES})
target_include_directories(geninterface PUBLIC ${NWC_INCLUDES})
target_compile_options(geninterface PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
