set(NWPW_INCLUDES
${PROJECT_SOURCE_DIR}/src/nwxc
${PROJECT_SOURCE_DIR}/src/nwpw/paw
${PROJECT_SOURCE_DIR}/src/nwpw/pspw
${PROJECT_SOURCE_DIR}/src/nwpw/band
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/include
${PROJECT_SOURCE_DIR}/src/nwpw/nwpwlib/utilities
)

set(NWPW_SOURCE_FILES
./nwpw_vectors_input.F
./nwpw_input.F
./nwpw_default.F
)


SET(NWPW_COMPILE_DEFS "")

if (USE_MPIF)
  list(APPEND NWPW_COMPILE_DEFS "MPI")
else()
  list(APPEND NWPW_COMPILE_DEFS "TCGMSG")
endif()

if (USE_MPIF4)
    list(APPEND NWPW_COMPILE_DEFS "MPI4")
endif()

if (USE_MLIB)
    list(APPEND NWPW_COMPILE_DEFS "MLIB")
endif()


add_subdirectory(paw)
add_subdirectory(band)
add_subdirectory(pspw)
add_subdirectory(ofpw)
add_subdirectory(nwpwlib)
add_subdirectory(libraryps)

add_library(nwpw STATIC ${NWPW_SOURCE_FILES} $<TARGET_OBJECTS:ofpw_opt> 
				 $<TARGET_OBJECTS:ofpw_un> $<TARGET_OBJECTS:nwpw_libps>)
#add_dependencies(nwpw ${NWC_DEPENDENCIES})
target_include_directories(nwpw PUBLIC ${NWC_INCLUDES} ${NWPW_INCLUDES})
target_compile_options(nwpw PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(nwpw PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

