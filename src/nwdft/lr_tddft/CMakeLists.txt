set(LR_TDDFT_SOURCE_FILES
tddft_dftinit.F
tddft_s2.F
tddft_iguess.F
tddft_simpleso.F
tddft_fold.F
tddft_davidson.F
tddft_finalize.F
tddft_subspace.F
tddft_energy.F
tddft_oscstr.F
tddft_transfm.F
tddft_dump.F
tddft_restart.F
tddft_contract.F
tddft_civecs_phase.F
tddft_init.F
tddft_residual.F
tddft_patch.F
tddft_input.F
tddft_checkpoint.F
tddft_epsilon.F
tddft_bdtensor.F
tddft_memory.F
tddft_analysis.F
tddft_diagon.F
tddft_schmidt.F
tddft_contract2.F
tddft_lanczos.F 
tddft_lanczos_tda.F
tddft_dipole.F 
tddft_prod_v.F
)

add_library(lr_tddft OBJECT ${LR_TDDFT_SOURCE_FILES})
target_compile_definitions(lr_tddft PUBLIC DEBUG_PRINT)
target_include_directories(lr_tddft PUBLIC ${NWC_INCLUDES} ${NWDFT_INCLUDES})
target_compile_options(lr_tddft PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(lr_tddft PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
