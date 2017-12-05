set(NWPW_PP_OPT_SOURCES
paw_atom/paw_atom.c
paw_atom/paw_atom_driver.c
paw_atom/paw_basis.c
paw_atom/paw_bisect.c
paw_atom/paw_comp_charge.c
paw_atom/paw_core.c
paw_atom/paw_dirac_exchange.c
paw_atom/paw_error_function.c
paw_atom/paw_get_inverse.c
paw_atom/paw_hartree.c
paw_atom/paw_ion.c
paw_atom/paw_kinetic_energy.c
paw_atom/paw_loggrid.c
paw_atom/paw_matrix_elements.c
paw_atom/paw_my_constants.c
paw_atom/paw_my_memory.c
paw_atom/paw_orbitals.c
paw_atom/paw_output.c
paw_atom/paw_pauli.c
paw_atom/paw_potential.c
paw_atom/paw_pred_cor.c
paw_atom/paw_scattering.c
paw_atom/paw_schrodin.c
paw_atom/paw_sdir.c
paw_atom/paw_utilities.c
paw_atom/paw_vosko.c
)

set(NWPW_PP_SOURCES
atom.c
becke_exchange.c
carter_parse.c
cpi_parse.c
debug.c
dft.c
dirac.c
dirac_exchange.c
gaussj.c
generate_rho_semicore.c
get_word.c
#gga_exchange.c
grids.c
hamann.c
hartree.c
hgh_parse.c
#integrate.c
loggrid.c
#lyp_correlation2.c
lyp_correlation.c
name.c
pauli.c
pbe_correlation.c
pbe_exchange.c
perdew_wang.c
perdew_zunger.c
pred_cor.c
psp1d.c
#psp.c
pspsolve.c
qmmm_parse.c
#rappe.c
revpbe_correlation.c
revpbe_exchange.c
rhamann.c
rpsp.c
#rpspsolve.c
rtroullier.c
schrodin.c
spline.c
teter_parse.c
troullier.c
# tt2.c
# tt.c
# ttt.c
vanderbilt.c
vosko.c
xpansion2.c
xpansion.c
zora.c
)


add_library(nwpwlib_pp_F OBJECT psp_generator_input.F)
target_include_directories(nwpwlib_pp_F PUBLIC ${NWC_INCLUDES} ${NWPWLIB_INCLUDES})
target_compile_definitions(nwpwlib_pp_F PUBLIC ${NWPW_COMPILE_DEFS}) 
target_compile_options(nwpwlib_pp_F PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})


add_library(nwpwlib_pp_un OBJECT ${NWPW_PP_SOURCES})
target_include_directories(nwpwlib_pp_un PUBLIC ${NWC_INCLUDES} ${NWPWLIB_INCLUDES})
target_compile_definitions(nwpwlib_pp_un PUBLIC ${NWPW_COMPILE_DEFS}) 
target_compile_options(nwpwlib_pp_un PUBLIC ${NWCHEM_C_COMPILE_FLAGS} ${NWCHEM_C_DEBUG_FLAGS})

add_library(nwpwlib_pp_opt OBJECT ${NWPW_PP_OPT_SOURCES})
target_include_directories(nwpwlib_pp_opt PUBLIC ${NWC_INCLUDES} ${NWPWLIB_INCLUDES})
target_compile_definitions(nwpwlib_pp_opt PUBLIC ${NWPW_COMPILE_DEFS}) # ${NWPWLIB_COMPILE_DEFS})
target_compile_options(nwpwlib_pp_opt PUBLIC ${NWCHEM_C_COMPILE_FLAGS})
target_compile_options(nwpwlib_pp_opt PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})
