set(NWPWLIB_INCLUDES
${PROJECT_SOURCE_DIR}/src/nwxc
${PROJECT_SOURCE_DIR}/src/nwpw
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/D3dB
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/TAMD
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/control
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/K1dB
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/brillouin_zone
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/linesearch
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/utilities
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/utilities/paw_utilities
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/io
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/unitcell
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/md_xs
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/md_xs/FEFF6L
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/ewald
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/FMM
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/Pneb
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/pseudopotential
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/pseudopotential/paw_atom
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/lattice
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/ilist
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/Parallel
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/Gaunt
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/MetaDynamics
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/include
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/C3dB
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/thermostats
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/DMatrix
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/ion
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/ion/shake
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/mapping
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/CMatrix
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/fractional
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/timing
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/D1dB

)

set(NWPW_C_SOURCES
./D3dB/hcurve.c
./D3dB/hilbert.c
./D3dB/olist.c
./utilities/nwpw_emachine.c
./utilities/nwpw_talker.c
./ion/shake/pspw_atom.c
./ion/shake/pspw_bond.c
./ion/shake/pspw_molecule.c
)

set(NWPWLIB_SOURCE_FILES
./D3dB/trans.F
./D3dB/D3dB-new.F
./D3dB/D3dB_pfft.F
./D3dB/generate_map_indexes.F
./TAMD/tamd.F
./TAMD/tamd_input.F
./control/control.F
./control/psi_get_header.F
./control/psi_get_version.F
./control/cpsi_get_header.F
./control/psi_get_ne.F
./control/get_psp_type.F
./K1dB/K1dB.F
./brillouin_zone/brillouin_zone_input.F
#./brillouin_zone/Pbrillioun.F
./brillouin_zone/brillioun.F
./TMatrix/dgamm_omp.F
./TMatrix/ffm_multiply_omp_group.F
./TMatrix/dscal_omp.F
./TMatrix/daxpy_omp.F
./TMatrix/dgemm_omp.F
./linesearch/Order_Values.F
./linesearch/linesearch.F
./utilities/matrix.f
./utilities/MOTION_nframes.F
./utilities/nwpw_list.F
./utilities/butter.F
./utilities/BA_memory.F
./utilities/analysis_expansion_coef.F
./utilities/pair_dist.F
./utilities/nwpw_kbpp_ray.F
./utilities/cellgeometry.F
./utilities/densityofstates.F
./utilities/simp.F
./utilities/nwpw_expression.F
./utilities/generate_unfolded_xyz.F
./utilities/spdf_name.F
./utilities/nwpw_matrix_invert.F
./utilities/psi_filefind.F
./utilities/cpmd_properties.F
./utilities/nwpw_scratch.F
./utilities/cell.F
./utilities/log_integrate.F
./utilities/kerker_G.F
./utilities/psi_data.F
./utilities/nwpw_kain.F
./utilities/nwpw_double_factorial.F
#./utilities/nwpw_4index.F
./utilities/nwpw_parse.F
./utilities/task_paulsmatrix.F
./utilities/beom.F
./utilities/pair_distab.F
./utilities/nwpw_spline.F
./utilities/is_sORd.F
./utilities/nwpw_diis.F
./utilities/nwpw_interp.F
./utilities/nwpw_filefind.F
./utilities/print_elapsed_time.F
./utilities/getfilling.F
./utilities/pair_dist2.F
./utilities/cpsi_data.F
./utilities/nwpw_scf_mixing.F
./utilities/CIF_write.F
./utilities/v_psi_filefind.F
./utilities/auto_corr.F
./utilities/nwpw_rdf.F
./utilities/btdb.F
./utilities/paw_utilities/nwpw_hypergeom1F1.F
./utilities/paw_utilities/nwpw_xc.F
./utilities/paw_utilities/nwpw_hartree.F
./utilities/paw_utilities/nwpw_gga.F
./utilities/paw_utilities/nwpw_compcharge.F
./utilities/paw_utilities/nwpw_gamma_function.F
./utilities/paw_utilities/nwpw_SpecialKummer.F
./utilities/paw_utilities/nwpw_gaunt.F
./utilities/paw_utilities/nwpw_WGaussian.F
./utilities/paw_utilities/nwpw_GaussBessel.F
./utilities/paw_utilities/nwpw_vosko.F
./utilities/paw_utilities/nwpw_gauss_weights.F
#./utilities/paw_utilities/nwpw_compcharge2.F
./unitcell/simulation_cell_input.F
./md_xs/md_xs_input.F
./md_xs/md_xs.F
./ewald/ewald.F
./ewald/ewald_strfac.F
./FMM/fmm.F
./Pneb/Pneb.F

./ion/center.f
./ion/fcoord.f
./ion/incell1.f
./ion/incell2.f
./ion/incell3.f
./ion/unfold.f
./lattice/G_init.F
./lattice/lattice.f
./lattice/c_G_init.F
./lattice/cstrfac.F
./lattice/balance.F
./lattice/mask_init.F
./lattice/cram.F
./lattice/pack.F
./lattice/cloak.F
./lattice/strfac.F
./lattice/c_balance.F
./ilist/nwpw_ilist.F
./Parallel/Parallel.F
#./Gaunt/nwpw_gaunt.F
./MetaDynamics/metadynamics_input.F
./MetaDynamics/metadynamics_aux.F
./MetaDynamics/metadynamics.F
./C3dB/C3dB-new.F
./C3dB/C3dB_pfft.F
./C3dB/c_trans.F
./thermostats/nose-hoover.F
./DMatrix/DMatrix.F
./ion/ion_scaling_atoms.F
./ion/xyz.F
./ion/seperate_pointcharge.F
./ion/ion.F
./ion/ion_FixIon.F
./ion/katm_list.F
./ion/seperate_molpsp.F

./ion/shake/shake_bonddiff.F
./ion/shake/shake_chain.F
./ion/shake/shake_analysis.F
./ion/shake/molecular_constraint.F
./ion/shake/shake_chain3.F
./ion/shake/shake_chain2.F
./ion/shake/molecular.F
./mapping/nwpw_mapping.F
./CMatrix/CMatrix.F
./fractional/frac_occ_input.F
./fractional/fractional.F
./fractional/frac_occ.F
./timing/nwpw_timing.F
./D1dB/Dne.F
./D1dB/Dne_omp.F
./D1dB/Dnex.F
./D1dB/skew.F
./D1dB/D1dB.F
)

#SET(NWPWLIB_COMPILE_DEFS "")

if (USE_MPIF)
  #list(APPEND NWPWLIB_COMPILE_DEFS "MPI")
  list(APPEND NWPWLIB_SOURCE_FILES ./CMatrix/CMatrix-mpi.F ./Parallel/Parallel-mpi.F
   ./D1dB/D1dB-mpi.F ./D3dB/D3dB-mpi.F
   ./K1dB/K1dB-mpi.F ./DMatrix/DMatrix-mpi.F ./C3dB/C3dB-mpi.F)
else()
  #add_definitions(-DTCGMSG)
  #list(APPEND NWPWLIB_COMPILE_DEFS "TCGMSG")
  list(APPEND NWPWLIB_SOURCE_FILES 
  ./CMatrix/CMatrix-tcgmsg.F ./Parallel/Parallel-tcgmsg.F ./C3dB/C3dB-tcgmsg.F
   ./D1dB/D1dB-tcgmsg.F ./D3dB/D3dB-tcgmsg.F 
   ./DMatrix/DMatrix-tcgmsg.F ./K1dB/K1dB-tcgmsg.F)
endif()


add_subdirectory(pseudopotential)

add_library(nwpwio OBJECT ./io/second.F ./io/nwpw_message.F)
target_include_directories(nwpwio PUBLIC ${NWC_INCLUDES})
target_compile_options(nwpwio PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwpwio PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

add_library(nwpwioC OBJECT ./io/compressed_io.c)
target_include_directories(nwpwioC PUBLIC ${NWC_INCLUDES})
target_compile_options(nwpwioC PUBLIC ${NWCHEM_C_COMPILE_FLAGS})
target_compile_options(nwpwioC PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})


add_library(nwpwlibC OBJECT ${NWPW_C_SOURCES})
target_include_directories(nwpwlibC PUBLIC ${NWC_INCLUDES} ${NWPWLIB_INCLUDES})
target_compile_definitions(nwpwlibC PUBLIC ${NWPW_COMPILE_DEFS}) # ${NWPWLIB_COMPILE_DEFS})
target_compile_options(nwpwlibC PUBLIC ${NWCHEM_C_COMPILE_FLAGS})
target_compile_options(nwpwlibC PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})


add_library(nwpwlib STATIC ${NWPWLIB_SOURCE_FILES} 
  $<TARGET_OBJECTS:nwpwio> $<TARGET_OBJECTS:nwpwioC> $<TARGET_OBJECTS:nwpwlibC>
  $<TARGET_OBJECTS:nwpwlib_pp_F> $<TARGET_OBJECTS:nwpwlib_pp_un> $<TARGET_OBJECTS:nwpwlib_pp_opt>)
#add_dependencies(nwpwlib ${NWC_DEPENDENCIES})
target_include_directories(nwpwlib PUBLIC ${NWC_INCLUDES} ${NWPWLIB_INCLUDES})
target_compile_definitions(nwpwlib PUBLIC ${NWPW_COMPILE_DEFS}) # ${NWPWLIB_COMPILE_DEFS})
target_compile_options(nwpwlib PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwpwlib PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

