
set(64TO32BLAS_SOURCE_FILES
./iyamax.F
./yspsvx.F
./yhseqr.F
./ygeevx.F
./ylacpy.F
./yposv.F
./ygesvd.F
./yscal.F
./ylange.F
./ylarnv.F
./ynrm2.F
./ylascl.F
./ygeev.F
./xcopy.F
./ygebak.F
./xaxpy.F
./xgesvd.F
./yger.F
./xheev.F
./yorghr.F
./ygemm.F
./ylagtf.F
./ylagts.F
./ygemv.F
./xdotc.F
./ysyev.F
./yrot.F
./ixamax.F
./ysterf.F
./yaxpy.F
./ytrevc.F
./xgemm.F
./ycopy.F
./ysygv.F
./ygehrd.F
./ygesv.F
./xscal.F
./ifily.F
./yasum.F
./yswap.F
./ygebal.F
./ydot.F
./ylaset.F
)

add_library(64to32blas STATIC ${64TO32BLAS_SOURCE_FILES})
#add_dependencies(64to32blas ${NWC_DEPENDENCIES})
target_include_directories(64to32blas PUBLIC ${NWC_INCLUDES})

if (BLAS_SIZE)
  if(${BLAS_SIZE} STREQUAL 8)
     target_compile_definitions(64to32blas PUBLIC USE_INTEGER8)
  endif()
endif()

target_compile_options(64to32blas PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(64to32blas PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

