
set(CR_EOMCCSD_T_ACT_SOURCE_FILES
q3rexpt2_act.F
cr_ccsd_t_E_act.F
cr_eomccsd_tasim.F
cr_eomt_simac.F
creomccsd_t_n2_act_mem.F
cr_eomccsd_t_act.F
ccsd_t_d_act.F
cr_ccsd_t_N_act.F
)

add_library(tce_cr_eomccsd_t_act STATIC ${CR_EOMCCSD_T_ACT_SOURCE_FILES})
#add_dependencies(tce_cr_eomccsd_t_act ${NWC_DEPENDENCIES})
target_include_directories(tce_cr_eomccsd_t_act PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_cr_eomccsd_t_act PUBLIC DEBUG_PRINT)

target_compile_options(tce_cr_eomccsd_t_act PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_cr_eomccsd_t_act PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
