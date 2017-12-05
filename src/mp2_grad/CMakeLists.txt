
set(MP2_GRAD_SOURCE_FILES
mp2_make_cphf_rhs.F
mp2_make_tot_dens.F
mp2_input.F
mp2_wai_uhf.F
mp2_make_nos.F
mp2_make_piy.F
mp2_grad_doit.F
mp2_grad.F
mp2_make_dip_mom.F
mp2_pstat.F
mp2_grad_utils.F
mp2_make_tot_wdens.F
mp2_energy.F
)

set(MP2_GRAD_OPT_SOURCE_FILES
mp2_back_transform.F 
mp2_read_tiajb.F
mp2_pijab_uhf.F 
mp2_wijab_uhf.F 
mp2_lai_uhf.F 
mp2_lai_fock_uhf.F 
mp2_memory.F 
mp2_make_tuhf.F
)


add_library(mp2_grad_opt OBJECT ${MP2_GRAD_OPT_SOURCE_FILES})
target_compile_definitions(mp2_grad_opt PUBLIC DEBUG_PRINT)
target_include_directories(mp2_grad_opt PUBLIC ${NWC_INCLUDES})
target_compile_options(mp2_grad_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(mp2_grad_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

add_library(mp2_grad STATIC ${MP2_GRAD_SOURCE_FILES} $<TARGET_OBJECTS:mp2_grad_opt>)
target_compile_definitions(mp2_grad PUBLIC DEBUG_PRINT)
#add_dependencies(mp2_grad ${NWC_DEPENDENCIES})
target_include_directories(mp2_grad PUBLIC ${NWC_INCLUDES})
target_compile_options(mp2_grad PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

