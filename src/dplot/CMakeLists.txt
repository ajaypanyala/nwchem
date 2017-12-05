
set(DPLOT_SOURCE_FILES
./dplot_input.F
./get_grid.F
./get_dens.F
./get_charges.F
./dmat.F
./ao_rdm_read.F
./dplot_dos.F
./get_transden.F
./get_diffden.F
./dplot.F
./errmem.F
./dplot_defgrid.F
./create_contour.F
./dplot_dump.F
)

add_library(dplot STATIC ${DPLOT_SOURCE_FILES})
#add_dependencies(dplot ${NWC_DEPENDENCIES})
target_include_directories(dplot PUBLIC ${NWC_INCLUDES})
target_compile_options(dplot PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(dplot PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
