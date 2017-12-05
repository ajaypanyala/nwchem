set(QMD_INCLUDES
${PROJECT_SOURCE_DIR}/src/nwmd
${PROJECT_SOURCE_DIR}/src/NWints/int
)

set(QMD_SOURCE_FILES
./qmd_thermostat.F
./qmd_main.F
./qmd_namd_init.F
./qmd_coup.F
./qmd_start.F
./qmd_freeze_inactive.F
./qmd_nucl_pos.F
./qmd_tdks_enmo.F
./qmd_randomize.F
./qmd_nucl_vel.F
./qmd_init.F
./task_qmd.F
./qmd_civecs_read.F
./qmd_write_rst.F
./qmd_elec_prop.F
./qmd_geom_store.F
./qmd_ran1.F
./qmd_nucl_dt.F
./qmd_move_civecs.F
./qmd_driver.F
./qmd_fix_momentum.F
./qmd_nucl_dt2.F
./qmd_read_rst.F
./qmd_set_den.F
./qmd_kinetic.F
./qmd_energy.F
./qmd_sumno.F
./qmd_gradient.F
./qmd_gasdev.F
./qmd_namd_input.F
./qmd_tddft_enmo.F
./qmd_coords_vels.F
./qmd_nucl_step.F
./qmd_gamdev.F
./qmd_write_trj.F
./qmd_input.F
./qmd_cavity.F
./qmd_property.F
)

add_library(qmd STATIC ${QMD_SOURCE_FILES})
#add_dependencies(qmd ${NWC_DEPENDENCIES})
target_include_directories(qmd PUBLIC ${NWC_INCLUDES} ${QMD_INCLUDES})
target_compile_options(qmd PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(qmd PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

