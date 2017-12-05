
set(DANGCHANG_SOURCE_FILES
./dc_initialize.F
./dc_setgeom.F
)

add_library(dcopt OBJECT ./dc_energy.F)
target_include_directories(dcopt PUBLIC ${NWC_INCLUDES})
target_compile_options(dcopt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(dcopt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

add_library(dangchang STATIC ${DANGCHANG_SOURCE_FILES} $<TARGET_OBJECTS:dcopt>)
target_include_directories(dangchang PUBLIC ${NWC_INCLUDES})
target_compile_options(dangchang PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

