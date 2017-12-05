set(NWC_SYMMETRY_INCLUDES
${PROJECT_SOURCE_DIR}/src/NWints/auxil
)

set(SYMMETRY_SOURCE_FILES
./sym_op_clsfy.F
./sym_put_geom.F
./sym_pr_ops.F
./sym_irrepname.F
./sym_apply_op.F
./sym_geom_prj.F
./sym_char_tab.F
./sym_vec_sym.F
./cross_prod.F
./sym_pr_all.F
./sym_cent_map.F
#./sym_op_mtab.F
./sym_tr_bs_op.F
./graveyard.F
./sym_nwc.F
./sym_map.F
./sym_get_cart.F
./sym_op_cname.F
./spgen.F
./sym_ops_get.F
./sym_abelian.F
./deter3.F
./sym_bs_irrep.F
./sym_num_ops.F
./sym_grp_name.F
./dosymops.F
./gensym_extra.F
./sym_g_sym.F
./sym_op_type.F
./wrcell.F
./sym_ap_cart.F
./sym_pr_ctab.F
./sym_mol.F
./sym_mo_adapt.F
./ludcmp.F
./sym_inv_op.F
./dctr.f  
./gensym.f  
./mprint.f  
./opprint.f
)

set(SYMMETRY_OPT_SOURCE_FILES
./sym_sym.F ./sym_mo_ap_op.F ./sym_bas_op.F ./sym_bs_ir_id.F
)

if (${CMAKE_HOST_SYSTEM_NAME} STREQUAL Linux OR ${CMAKE_HOST_SYSTEM_NAME} STREQUAL Darwin)
    set(SYMMETRY_OPT_SOURCE_FILES ${SYMMETRY_OPT_SOURCE_FILES} ./sym_sh_pair.F)
endif ()

add_library(symmetry_opt OBJECT ${SYMMETRY_OPT_SOURCE_FILES})
target_include_directories(symmetry_opt PUBLIC ${NWC_INCLUDES} ${NWC_SYMMETRY_INCLUDES})
target_compile_options(symmetry_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(symmetry_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

add_library(symmetry STATIC ${SYMMETRY_SOURCE_FILES} $<TARGET_OBJECTS:symmetry_opt>)
#add_dependencies(symmetry ${NWC_DEPENDENCIES})
target_include_directories(symmetry PUBLIC ${NWC_INCLUDES} ${NWC_SYMMETRY_INCLUDES})
target_compile_options(symmetry PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})


