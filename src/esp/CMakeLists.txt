
set(ESP_SOURCE_FILES
./esp_esp.F
./esp_jacobi.F
./esp_wrtinp.F
./esp_dalloc.F
./esp_fit.F
./esp_getg.F
./esp_getden.F
./esp_out.F
./esp_elemom.F
./esp_rms.F
./esp_matinv.F
./esp_grid.F
./esp_input.F
./esp_fastesp.F
./esp_plt.F
./esp_rtdbin.F
./esp.F
./esp_denmat.F
./esp_grdcnt.F
./esp_getx.F
./esp_malloc.F
)

add_library(esp STATIC ${ESP_SOURCE_FILES})
#add_dependencies(esp ${NWC_DEPENDENCIES})
target_include_directories(esp PUBLIC ${NWC_INCLUDES})
target_compile_options(esp PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(esp PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
