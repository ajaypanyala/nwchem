set(NWINTS_INCLUDES
${PROJECT_SOURCE_DIR}/src/NWints/int
${PROJECT_SOURCE_DIR}/src/NWints/rel
${PROJECT_SOURCE_DIR}/src/NWints/ecp
${PROJECT_SOURCE_DIR}/src/NWints/auxil
${PROJECT_SOURCE_DIR}/src/NWints/hondo
${PROJECT_SOURCE_DIR}/src/NWints/ints_sp
)

set(NWINTS_SOURCE_FILES
./intgrl_input.F
)

set(USE_TEXAS ON)

add_subdirectory(dk)
add_subdirectory(api)
add_subdirectory(int)
add_subdirectory(ecp)
add_subdirectory(rel)
add_subdirectory(auxil)
add_subdirectory(deriv)
add_subdirectory(hondo)
add_subdirectory(texas)
add_subdirectory(ints_sp)
#add_subdirectory(simint)

add_library(NWints STATIC ${NWINTS_SOURCE_FILES} $<TARGET_OBJECTS:nwints_api_opt> 
$<TARGET_OBJECTS:nwints_api_un> $<TARGET_OBJECTS:nwints_hondo_opt> $<TARGET_OBJECTS:nwints_hondo_un>
$<TARGET_OBJECTS:nwints_texas_un> $<TARGET_OBJECTS:nwints_texas_opt> $<TARGET_OBJECTS:nwints_int>
$<TARGET_OBJECTS:nwints_ecp_un> $<TARGET_OBJECTS:nwints_ecp_opt> $<TARGET_OBJECTS:nwints_deriv>
$<TARGET_OBJECTS:nwints_sp_un> $<TARGET_OBJECTS:nwints_sp> $<TARGET_OBJECTS:nwints_dk>
$<TARGET_OBJECTS:nwints_auxil_un> $<TARGET_OBJECTS:nwints_auxil_opt>
$<TARGET_OBJECTS:nwints_rel_un> $<TARGET_OBJECTS:nwints_rel_opt> 
)

target_include_directories(NWints PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(NWints PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(NWints PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
