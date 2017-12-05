
set(CROP_SOURCE_FILES
tce_crop.F
)

add_library(tce_crop STATIC ${CROP_SOURCE_FILES})
#add_dependencies(tce_crop ${NWC_DEPENDENCIES})
target_include_directories(tce_crop PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_crop PUBLIC DEBUG_PRINT)

target_compile_options(tce_crop PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_crop PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
