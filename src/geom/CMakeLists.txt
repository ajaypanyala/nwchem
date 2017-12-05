
set(GEOM_SOURCE_FILES
./geom_print_ecce.F
./geom_fragment.F
./geom_input.F
./geom.F
./geom_getsym.F
./geom_freeze.F
./geom_1d.F
./geom_2d.F
./geom_checksum.F
./geom_3d.F
./geom_input2.F
./geom_numcore.F
)

add_library(geom_opt OBJECT ./geom_hnd.F)
target_include_directories(geom_opt PUBLIC ${NWC_INCLUDES})
target_compile_options(geom_opt PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(geom_opt PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

add_library(geom STATIC ${GEOM_SOURCE_FILES} $<TARGET_OBJECTS:geom_opt>)
#add_dependencies(geom ${NWC_DEPENDENCIES})
target_include_directories(geom PUBLIC ${NWC_INCLUDES})

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    if(USE_OPENMP)
        target_compile_options(geom PUBLIC -fno-recursive -fno-openmp) 
	target_compile_options(geom_opt PUBLIC -fno-recursive -fno-openmp) 
    endif()
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    target_compile_options(geom PUBLIC -save) 
    target_compile_options(geom_opt PUBLIC -save) 
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "XL")
    target_compile_options(geom PUBLIC -qsave) 
    target_compile_options(geom_opt PUBLIC -qsave) 
endif()

target_compile_options(geom PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

