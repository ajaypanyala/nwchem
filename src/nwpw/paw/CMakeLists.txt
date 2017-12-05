set(PAW_INCLUDES
${PROJECT_SOURCE_DIR}/src/nwxc
${PROJECT_SOURCE_DIR}/src/nwpw

${PROJECT_SOURCE_DIR}/src/nwpw/paw
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_density
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_proj
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_mult
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_core
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_vloc
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_matrix
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_utils
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_basis
${PROJECT_SOURCE_DIR}/src/nwpw/paw/cell_optimize
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_kinetic
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_gaunt
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_ovlp
${PROJECT_SOURCE_DIR}/src/nwpw/paw/cpsd
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_psi
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_ma
${PROJECT_SOURCE_DIR}/src/nwpw/paw/cgsd
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_special_functions
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_comp_charge
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_geometry
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_nonlocal
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_ion_potential
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_hartree
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_comp_pot
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_Grsm_list
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_force
${PROJECT_SOURCE_DIR}/src/nwpw/paw/paw_xc

${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/include
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/fractional
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/utilities
)

set(PAW_SOURCE_FILES
./paw_density/paw_density.F
./paw_proj/paw_proj.F
./paw_mult/paw_mult_interaction.F
./paw_mult/paw_mult.F
./paw_mult/paw_mult_energy.F
./paw_core/paw_matrix_core.F
./paw_core/paw_energy_core_atom.F
./paw_vloc/paw_vloc.F
./paw_vloc/paw_matrix_vloc.F
./paw_vloc/paw_energy_vloc_atom.F
./paw_matrix/paw_comp_charge_matrix.F
./paw_utils/paw_utils.F
./paw_utils/io_utils.F
./paw_utils/paw_spher_grid.F
./paw_utils/chkerr.F
./paw_utils/paw_gauss_weights.F
./paw_utils/bessel_transform.F
./paw_basis/paw_basis.F
./paw_basis/paw_small_basis.F
./paw_basis/paw_basis_data.F
./cell_optimize/paw_cell_optimize.F
./cell_optimize/paw_stress_fd.F
./cell_optimize/paw_update_unitcell.F
./cell_optimize/paw_stress_numerical.F
./paw_kinetic/paw_energy_kin_atom.F
./paw_kinetic/paw_matrix_kin.F
./paw_gaunt/paw_gaunt.F
./paw_ovlp/paw_ovlp.F
./paw_ovlp/paw_overlap_matrix.F
./cpsd/paw_inner_loop.F
./cpsd/paw_inner_loop_md.F
./cpsd/psi_marat_write.F
./cpsd/paw_output.F
./cpsd/paw_cpmd.F
./cpsd/paw_sd.F
./paw_psi/paw_geodesic.F
./paw_psi/paw_psi.F
./paw_psi/paw_electron.F
./paw_psi/paw_psi2.F
./paw_psi/paw_psi_H.F
./paw_psi/paw_psi_new.F
./paw_ma/paw_ma.F
./cgsd/paw_gradient.F
./cgsd/paw_cgsd_energy.F
./cgsd/paw_bfgsminimize.F
./cgsd/paw_cgsd_noit_energy.F
./cgsd/paw_noit_energy.F
./cgsd/paw_energy.F
./cgsd/paw_nominimize.F
./cgsd/paw_sdminimize.F
./cgsd/paw_minimizer.F
./cgsd/paw_stress.F
./cgsd/paw_cgminimize.F
./paw_special_functions/paw_spher_func.F
./paw_special_functions/paw_double_factorial.F
./paw_special_functions/paw_gamma_function.F
./paw_special_functions/paw_hypergeom1F1.F
./paw_comp_charge/paw_comp_charge.F
./paw_nonlocal/paw_nonlocal.F
./paw_ion_potential/paw_matrix_ion.F
./paw_ion_potential/paw_energy_ion_atom.F
./paw_hartree/paw_matrix_hartree.F
./paw_hartree/paw_energy_hartree_atom.F
./paw_hartree/paw_pot_hartree.F
./paw_comp_pot/paw_pot_comp.F
./paw_comp_pot/paw_matrix_comp_pot.F
./paw_Grsm_list/paw_lmbfgs.F
./paw_Grsm_list/paw_Grsm_list.F
./paw_force/paw_force.F
./paw_xc/paw_gga.F
./paw_xc/paw_xc.F
./paw_xc/paw_vosko.F
./paw_xc/paw_dirac_x.F
)



add_library(paw_un OBJECT ./task_paw.F ./paw_utils/integrate.F)
target_include_directories(paw_un PUBLIC ${NWC_INCLUDES} ${PAW_INCLUDES})
target_compile_definitions(paw_un PUBLIC ${NWPW_COMPILE_DEFS})
target_compile_options(paw_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})


add_library(paw STATIC ${PAW_SOURCE_FILES} $<TARGET_OBJECTS:paw_un>)
#add_dependencies(paw ${NWC_DEPENDENCIES})
target_include_directories(paw PUBLIC ${NWC_INCLUDES} ${PAW_INCLUDES})
target_compile_definitions(paw PUBLIC ${NWPW_COMPILE_DEFS})
target_compile_options(paw PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(paw PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
