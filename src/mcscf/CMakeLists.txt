
set(MSCF_INCLUDES
${PROJECT_SOURCE_DIR}/src/mcscf/detci
)

set(MCSCF_SOURCE_FILES
./detci/detci_driver.F
./detci/detci_aux.F
./detci/detci_intf.F
./mcscf_etr.F
./mcscf_fock.F
./mcscf_orb.F
./mcscf.F
./mcscf_debug.F
./mcscf_ohv1e.F
./mcscf_ohv2e.F
./mcscf_fcore.F
./mcscf_ohvao.F
./mcscf_hdiag.F
./mcscf_ohvmo.F
./mcscf_prcnd.F
./mcscf_aux.F
./mcscf_input.F
./hmat.F
./mcscf_ohv.F
./mcscf_ci.F
./mcscf_fdiff.F
./linesrch.F
./mcscf_b.F
)

add_library(mcscf_opt OBJECT ./detci/detci_dens.F ./detci/detci_diag.F ./detci/detci_sigma.F ./detci/detci_spin.F)
target_include_directories(mcscf_opt PUBLIC ${NWC_INCLUDES} ${MSCF_INCLUDES})
target_compile_definitions(mcscf_opt PUBLIC OLD_DAVIDSON)
target_compile_options(mcscf_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(mcscf_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})


add_library(mcscf STATIC ${MCSCF_SOURCE_FILES} $<TARGET_OBJECTS:mcscf_opt>)
#add_dependencies(mcscf ${NWC_DEPENDENCIES})
target_include_directories(mcscf PUBLIC ${NWC_INCLUDES} ${MSCF_INCLUDES})
target_compile_definitions(mcscf PUBLIC OLD_DAVIDSON)
target_compile_options(mcscf PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
