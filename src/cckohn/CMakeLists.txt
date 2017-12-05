
set(CCKOHN_SOURCE_FILES
./fforth2.F
./fforth1.F
./obftrans.F
./top.F
./denom.F
./cckohn.F
./bforthog.F
./cknew.F
./hbftrans.F
./cckohn_input.F
./bessel/gg.F
./bessel/ubesplg.F
./bessel/zero.F
./bessel/splinec.F
./bessel/sjymec.F
./bessel/outward.F
./bessel/cckohn_spline.F
)

add_library(cckohn STATIC ${CCKOHN_SOURCE_FILES})
#add_dependencies(cckohn ${NWC_DEPENDENCIES})
target_include_directories(cckohn PUBLIC ${NWC_INCLUDES})
target_compile_options(cckohn PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})
#target_compile_options(cckohn PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

