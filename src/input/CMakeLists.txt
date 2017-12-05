set(INPUT_SOURCE_FILES
./input_parse.F
./get_inp_file.F
./input_mem.F
./memory_def.F
./speech_input.F
./memory_input.F
./input_title.F
./input_set.F
)

add_library(nwc_input STATIC ${INPUT_SOURCE_FILES})
target_include_directories(nwc_input PUBLIC ${NWC_INCLUDES})
if(MRCC_METHODS)
    target_compile_definitions(nwc_input PUBLIC MRCC_METHODS)
endif()

target_compile_options(nwc_input PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(nwc_input PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
