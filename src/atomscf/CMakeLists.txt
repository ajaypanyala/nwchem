set(ATOMSCF_INCLUDES
${PROJECT_SOURCE_DIR}/src/NWints/ecp
)

set(ATOMSCF_SOURCE_FILES
./trafsd.F
./shalfd.F
./tramad.F
./hamild.F
./denat.F
./tracd.F
./pdfded.F
./densid_z.F
./oeigd.F
./teigd.F
./orderd.F
./jacod.F
./square.F
./atomd.F
./atcond.F
./cmergd.F
./guess_mem.F
./starcd.F
./printbasis.F
#./osatod.F
./denmad.F
./writel.F
./creded.F
./atscf_atscf.F
./datoms.F
./atscf_ecp.F
./densid.F
./outpud.F
)

add_library(atomscf STATIC ${ATOMSCF_SOURCE_FILES})
#add_dependencies(atomscf ${NWC_DEPENDENCIES})
target_include_directories(atomscf PUBLIC ${NWC_INCLUDES} ${ATOMSCF_INCLUDES})

target_compile_options(atomscf PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(atomscf PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

