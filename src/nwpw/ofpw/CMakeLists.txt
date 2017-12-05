set(OFPW_INCLUDES
${PROJECT_SOURCE_DIR}/src/nwpw/ofpw
${PROJECT_SOURCE_DIR}/src/nwpw/ofpw/minimizer
${PROJECT_SOURCE_DIR}/src/nwpw/ofpw/cpsd
${PROJECT_SOURCE_DIR}/src/nwpw/ofpw/lib
${PROJECT_SOURCE_DIR}/src/nwpw/ofpw/lib/electron
${PROJECT_SOURCE_DIR}/src/nwpw/ofpw/lib/WGC
${PROJECT_SOURCE_DIR}/src/nwpw/ofpw/lib/chi
${PROJECT_SOURCE_DIR}/src/nwpw/ofpw/lib/ke
)

set(NWPW_OFPW_SOURCE_FILES
./ofpw_energy.F
./task_ofpw.F
./minimizer/ofpw_cgsd_energy.F
./minimizer/ofpw_minimizer.F
./minimizer/ofpw_cgminimize.F
./cpsd/inner_loop_ofpw.F
./cpsd/ofpwsd.F
./cpsd/chi_lmbda.f
./lib/electron/chi_H.F
./lib/electron/chi_Hv4.F
./lib/WGC/wgc.F
./lib/chi/chi_read.F
./lib/chi/chi_new.F
./lib/chi/chi_expander.F
./lib/chi/chi_write.F
)


add_library(ofpw_un OBJECT ./lib/ke/v_thomasfermi.F)
target_include_directories(ofpw_un PUBLIC ${NWC_INCLUDES} ${NWPW_INCLUDES} ${OFPW_INCLUDES})
target_compile_options(ofpw_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(ofpw_opt OBJECT ${NWPW_OFPW_SOURCE_FILES})
target_include_directories(ofpw_opt PUBLIC ${NWC_INCLUDES} ${NWPW_INCLUDES} ${OFPW_INCLUDES})
target_compile_options(ofpw_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(ofpw_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
