set(NWPW_LIBPS_SOURCE_FILES
nwpw_set_library.F
nwpw_libfile.F
nwpw_library.F
)


add_library(nwpw_libps OBJECT ${NWPW_LIBPS_SOURCE_FILES})
target_include_directories(nwpw_libps PUBLIC ${NWC_INCLUDES} ${NWPW_INCLUDES})
target_compile_options(nwpw_libps PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwpw_libps PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})