set(BAND_INCLUDES
${PROJECT_SOURCE_DIR}/src/nwxc
${PROJECT_SOURCE_DIR}/src/nwpw

${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/include
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/utilities

${PROJECT_SOURCE_DIR}/src/nwpw/band
${PROJECT_SOURCE_DIR}/src/nwpw/band/structure
${PROJECT_SOURCE_DIR}/src/nwpw/band/minimizer
${PROJECT_SOURCE_DIR}/src/nwpw/band/cell_optimize
${PROJECT_SOURCE_DIR}/src/nwpw/band/dplot
${PROJECT_SOURCE_DIR}/src/nwpw/band/cpsd
${PROJECT_SOURCE_DIR}/src/nwpw/band/analysis
${PROJECT_SOURCE_DIR}/src/nwpw/band/makepsi
${PROJECT_SOURCE_DIR}/src/nwpw/band/lib
${PROJECT_SOURCE_DIR}/src/nwpw/band/lib/hfx
${PROJECT_SOURCE_DIR}/src/nwpw/band/lib/borbs
${PROJECT_SOURCE_DIR}/src/nwpw/band/lib/cpsp
${PROJECT_SOURCE_DIR}/src/nwpw/band/lib/ewald
${PROJECT_SOURCE_DIR}/src/nwpw/band/lib/electron
${PROJECT_SOURCE_DIR}/src/nwpw/band/lib/psi
${PROJECT_SOURCE_DIR}/src/nwpw/band/lib/Grsm_list
${PROJECT_SOURCE_DIR}/src/nwpw/band/lib/coulomb
${PROJECT_SOURCE_DIR}/src/nwpw/band/lib/brillioun
${PROJECT_SOURCE_DIR}/src/nwpw/band/lib/ke
${PROJECT_SOURCE_DIR}/src/nwpw/band/kbpp
)

set(BAND_SOURCE_FILES
./kbpp/setup_kbpp_band.f
./band_stress.F
./band_energy.F
./task_band.F
./band_gradient.F
./band_noit_energy.F
./structure/band_fft_structure.F
./structure/band_kp_structure.F
./structure/band_interpolate_structure.F
./structure/band_structure.F
./minimizer/c_cgsd_noit_energy.F
./minimizer/c_bybminimize.F
./minimizer/c_cgsd_energy.F
./minimizer/c_cgsd_excited.F
./minimizer/band_minimizer.F
./minimizer/c_bfgsminimize.F
./minimizer/c_cgsd_energy_gradient.F
./minimizer/c_cgsd_energy_stress.F
./minimizer/c_nominimize.F
./minimizer/c_cgminimize.F
./minimizer/c_bybminimize2.F
./minimizer/c_sdminimize.F
./cell_optimize/band_update_unitcell.F
./cell_optimize/band_cell_optimize.F
./cell_optimize/band_stress_numerical.F
./cell_optimize/band_stress_fd.F
./dplot/band_dplot_input.F
./dplot/band_dplot_iteration.F
./dplot/band_dplot.F
./cpsd/band_inner_loop_md.F
./cpsd/band_inner_loop.F
./cpsd/band_cpmd.F
./cpsd/band_sd.F
./analysis/orb_pop_borb.F
./analysis/band_analysis.F
./makepsi/c_v_wvfnc_new.F
./makepsi/pspw_reformat_c_wvfnc.F
./makepsi/c_wvfnc_expander.F
./makepsi/c_wvfnc_new.F
./lib/hfx/band_hfx.F
./lib/borbs/borb_projector.F
./lib/borbs/borbs.F
#./lib/cpsp/cpsp_U.F
./lib/cpsp/c_semicore_force_a.F
./lib/cpsp/cpsp_stress.F
./lib/cpsp/cpsp.F
./lib/cpsp/c_semicore.F
./lib/cpsp/cpsp_projector.F
#./lib/cpsp/c_semicore_force_b.F
./lib/cpsp/c_semicore_force.F
./lib/cpsp/cpsp_relativistic.F
./lib/ewald/c_ewald.F
./lib/electron/c_electron.F
./lib/electron/cpsi_H.F
./lib/psi/Stiefel.F
./lib/psi/cpsi_write.F
./lib/psi/cpsi_read.F
./lib/psi/cpsi.F
./lib/psi/cpsi_spin2.F
./lib/psi/cpsi_new.F
./lib/psi/v_cpsi_new.F
./lib/psi/cpsi_lmbda.F
./lib/psi/v_cpsi_write.F
./lib/psi/BGrsm.F
./lib/psi/c_geodesic.F
./lib/psi/v_cpsi_read.F
./lib/psi/c_check.F
./lib/psi/cpsi_KS.F
./lib/Grsm_list/band_Grsm_list.F
./lib/Grsm_list/band_lmbfgs.F
./lib/coulomb/c_coulomb_screened.F
./lib/coulomb/c_coulomb.F
#./lib/brillioun/brillioun.F - not compiled
./lib/ke/cke.F
./kbpp/integrate_kbpp_e_band_nonlocal.F
./kbpp/integrate_kbpp_band_ray.F
./kbpp/cpsp_stress_formatter_auto.F
./kbpp/integrate_kbpp_e_band_ray.F
./kbpp/integrate_kbpp_e_band_stress_nonlocal.F
./kbpp/kbpp_band_ray.F
./kbpp/kbpp_e_band.F
./kbpp/integrate_kbpp_band_orb.F
./kbpp/integrate_kbpp_e_band_stress_nonlocal_new.F
./kbpp/kbpp_band_stress.F
./kbpp/integrate_kbpp_band_stress_nonlocal.F
./kbpp/Spherical_band.F
./kbpp/kbpp_e_band_stress.F
./kbpp/integrate_kbpp_band_local.F
./kbpp/cpsp_formatter_auto.F
./kbpp/kbpp_band_orb.F
./kbpp/integrate_kbpp_band_nonlocal.F
./kbpp/integrate_kbpp_e_band_nonlocal_new.F
./kbpp/kbpp_band.F
./kbpp/integrate_kbpp_band_local_new.F
./kbpp/HGH_band_nonlocal.F
./kbpp/integrate_kbpp_e_band_stress_ray.F
./kbpp/integrate_kbpp_band_stress_nonlocal_new.F
./kbpp/HGH_band_local.F
./kbpp/integrate_kbpp_band_stress_ray.F
./kbpp/Tesseral_band.F
./kbpp/integrate_kbpp_e_band_local.F
./kbpp/borbs_formatter.F
./kbpp/rkbpp_band.F
./kbpp/integrate_kbpp_band_nonlocal_new.F
./kbpp/integrate_kbpp_e_band_local_new.F
./kbpp/integrate_kbpp_band_stress_local_new.F
./kbpp/integrate_kbpp_e_band_stress_local_new.F
./kbpp/integrate_kbpp_e_band_stress_local.F
./kbpp/hghpp_band.F
./kbpp/integrate_kbpp_band_stress_local.F
)


add_library(band_un OBJECT ./band_input.F)
target_include_directories(band_un PUBLIC ${NWC_INCLUDES} ${BAND_INCLUDES})
target_compile_definitions(band_un PUBLIC ${NWPW_COMPILE_DEFS})
target_compile_options(band_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(band STATIC ${BAND_SOURCE_FILES} $<TARGET_OBJECTS:band_un>)
#add_dependencies(band ${NWC_DEPENDENCIES})
target_include_directories(band PUBLIC ${NWC_INCLUDES} ${BAND_INCLUDES})
target_compile_definitions(band PUBLIC ${NWPW_COMPILE_DEFS})
target_compile_options(band PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(band PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

