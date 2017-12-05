set(QMMM_INCLUDES
${PROJECT_SOURCE_DIR}/src/nwmd
)

set(QMMM_SOURCE_FILES
./qmmm_cons.F
./qmmm_bq.F
./task_qmmm_abc.F
./qmmm_bq_data.F
./qmmm_input.F
./qmmm_esp.F
./qmmm_data.F
./task_qmmm_fep.F
./qmmm_output.F
./task_qmmm_dynamics.F
./qmmm_main.F
./task_qmmm_energy.F
./mm_interface.F
./task_qmmm_optimize.F
./task_qmmm_gradient.F
./qmmm_utils.F
./qmmm_check_forces.F
./qmmm_link.F
./qmmm_lbfgsb.F
./qmmm_property.F
./qmmm_pes.F
./qmmm.F
./mm_utils.F
./qmmm_geom.F
)

add_library(qmmm STATIC ${QMMM_SOURCE_FILES})
#add_dependencies(qmmm ${NWC_DEPENDENCIES})
target_include_directories(qmmm PUBLIC ${NWC_INCLUDES} ${QMMM_INCLUDES})
if (${CMAKE_HOST_SYSTEM_NAME} STREQUAL Linux OR ${CMAKE_HOST_SYSTEM_NAME} STREQUAL Darwin)
    target_compile_definitions(qmmm PUBLIC NEED_LOC)
endif ()

target_compile_options(qmmm PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(qmmm PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
