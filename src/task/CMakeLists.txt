
set(TASK_SOURCE_FILES
./task_optimize.F
./task_bsse.F
./task_dntmc.F
./task_rfcf.F
./task_property.F
./task_dynamics.F
./task_shell_input.F
./task_save_state.F
./task_gradient.F
./task_vscf.F
./task_freq.F
./task_thermo.F
./task_energy.F
./task_et.F
./task_hessian.F
./task_input.F
./task_ncc.F
./task_jefftce.F
./task_num_grad.F
./task.F
)

add_library(nwc_task STATIC ${TASK_SOURCE_FILES})
target_include_directories(nwc_task PUBLIC ${NWC_INCLUDES})
target_compile_options(nwc_task PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(nwc_task PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
