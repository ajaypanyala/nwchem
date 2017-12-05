
#file(COPY ${PROJECT_SOURCE_DIR}/src/nwxc/nwxc.fh DESTINATION ${PROJECT_SOURCE_DIR}/src/include/)

set(NWXC_SOURCE_FILES
./nwxc_vdw3a.F
./nwxc_vdw3b.F
./nwxc_vdw3c.F
./nwxc_vdw3d.F
#./nwxc_print.F
#./nwxc.F 
./nwxc_eval.F
./nwxc_add.F 
./nwxc_block.F 
./nwxc_query.F
./nwxc_nwchem.F 
#./nwxcm_eval.F
# ./nwxc_vdw_main.F 
# ./nwxc_vdw_e.F 
# ./nwxc_vdw_der.F 
# ./nwxc_vdw_hess.F 
# ./nwxc_vdw_util.F 
# ./nwxcm_heaviside.F
)

add_library(nwxc STATIC ${NWXC_SOURCE_FILES})
#add_dependencies(nwxc ${NWC_DEPENDENCIES})
target_include_directories(nwxc PUBLIC ${NWC_INCLUDES} ./nwad)
