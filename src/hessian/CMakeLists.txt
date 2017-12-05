set(HESSIAN_INCLUDES
${PROJECT_SOURCE_DIR}/src/nwdft/xc
${PROJECT_SOURCE_DIR}/src/nwdft/include
${PROJECT_SOURCE_DIR}/src/NWints/int
${PROJECT_SOURCE_DIR}/src/NWints/hondo
)

set(HESSIAN_SOURCE_FILES
#./rhf_hessian.F
./scf_hessian.F
./analytic/hess_check.F
./analytic/hess_anal.F
./analytic/twodd_cont.F
./analytic/uhf_fock_setup.F
./analytic/hess_sxsy.F
./analytic/hess_init.F
./analytic/nucdd_cont.F
./analytic/uhf_fock_finish.F
./analytic/hess_wrt.F
./analytic/print_dens_bit.F
./analytic/twodd_coul_ex.F
./analytic/shell_fock_build.F
./analytic/sx_dra.F
./analytic/hess_restart.F
./analytic/hess_energyinit.F
./analytic/hess_hssout.F
./analytic/shell_fock_build_cam.F
./analytic/get_wdens.F
./analytic/hess_utils.F
./analytic/hess_pxfy.F
./analytic/oned_calc.F
./analytic/hess_dipder.F
./analytic/hess_tidy.F
./analytic/onedd_calc.F
./analytic/hess_print.F
./analytic/nucdd_calc.F
./analytic/shell_fock_build_ext.F
./analytic/hess_cphf.F
./analytic/onedd_cont.F
./analytic/dft/xc_init_index.F
./analytic/dft/xc_cpks_coeff.F
./analytic/dft/xc_exit_index.F
./analytic/dft/xc_nucder_gen.F
./analytic/dft/dft_3dacc.F
)

add_library(hessian_un OBJECT ./analytic/dft/xc_d2expl.F ./hess_blas.F ./hess_input.F)
target_include_directories(hessian_un PUBLIC ${NWC_INCLUDES} ${HESSIAN_INCLUDES})
target_compile_options(hessian_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(hessian STATIC ${HESSIAN_SOURCE_FILES} $<TARGET_OBJECTS:hessian_un>)
#add_dependencies(hessian ${NWC_DEPENDENCIES})
target_include_directories(hessian PUBLIC ${NWC_INCLUDES} ${HESSIAN_INCLUDES})

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(hessian PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} -ffixed-line-length-none)
else()
  target_compile_options(hessian PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
endif()
target_compile_options(hessian PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

