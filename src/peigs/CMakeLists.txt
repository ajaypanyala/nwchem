add_definitions(-DTCGMSG)

set(PEIGS_INCLUDES
${PROJECT_SOURCE_DIR}/src/peigs
${PROJECT_SOURCE_DIR}/src/peigs/h
${PROJECT_SOURCE_DIR}/src/peigs/src
${PROJECT_SOURCE_DIR}/src/peigs/src/c
${PROJECT_SOURCE_DIR}/src/peigs/src/f77
${PROJECT_SOURCE_DIR}/src/peigs/new_src
${PROJECT_SOURCE_DIR}/src/peigs/h
${PROJECT_SOURCE_DIR}/src/peigs/comm
${PROJECT_SOURCE_DIR}/src/peigs/ctof
)

set(PEIGS_C_SOURCE_FILES
./ctof/bortho_f.c
./ctof/choleski_f.c
./ctof/inverse_f.c
./ctof/mxm25_f.c
./ctof/mxm5x_f.c
./ctof/mxm88_f.c
./ctof/onenorm_f.c
./ctof/ortho_f.c
./ctof/pdspev_f.c
./ctof/pdspevx_f.c
./ctof/pdspgv_f.c
./ctof/pdspgvx_f.c
#./ctof/pdsptri_f.c
./ctof/resid_f.c
./ctof/residual_f.c
./ctof/sonenorm_f.c
./ctof/tresid_f.c
./src/c/b_ortho.c
./src/c/chol_bcst.c
./src/c/choleski9.c
./src/c/ci_entry.c
#./src/c/clustrf3.c
#./src/c/clustrf4.c
./src/c/clustrf5.c
./src/c/clustrfix.c
./src/c/conjug22.c
./src/c/de_sym.c
./src/c/exit2.c
./src/c/exit.c
./src/c/forLL.c
./src/c/forLU.c
./src/c/gmax.c
./src/c/inverse.c
#./src/c/inv_it3.c
#./src/c/inv_it44.c
#./src/c/inv_it4.c
./src/c/inv_it5.c
#./src/c/inv_it6.c
./src/c/lu_mxm2.c
./src/c/mapchk.c
./src/c/mapdif1.c
./src/c/mapdif.c
./src/c/matmax.c
./src/c/mdif2b.c
./src/c/mdiff1.c
./src/c/mdiff2.c
#./src/c/memreq2.c
./src/c/memreq.c
#./src/c/memreq_eval.c
./src/c/memreq_f.c
#./src/c/mgs4.c
./src/c/mgs5.c
./src/c/mgscs.c
./src/c/mxm25.c
./src/c/mxm5x.c
./src/c/mxm88.c
./src/c/mxm_ll1.c
./src/c/onenorm.c
#./src/c/ortho2.c
./src/c/ortho.c
#./src/c/pdcomplex.c
./src/c/pdiff.c
#./src/c/pdspeval.c
./src/c/pdspev_c.c
./src/c/pdspevx.c
./src/c/pdspgv_c.c
./src/c/pdspgvx.c
#./src/c/pdsptri2.c
./src/c/pdsptri.c
./src/c/peigs_dlasq1.c
./src/c/peigs_tldlfact.c
./src/c/pgexit.c
./src/c/pipe_fut.c
./src/c/pmmlsl2.c
./src/c/pmmLSL.c
./src/c/pmmLUL.c
./src/c/pscale.c
./src/c/pstebz10.c
./src/c/pstebz11.c
#./src/c/pstebz9.c
#./src/c/pstein3.c
#./src/c/pstein44.c
./src/c/pstein4.c
./src/c/pstein5.c
./src/c/pxerbla.c
./src/c/qsort.c
./src/c/reducemap.c
./src/c/resid.c
./src/c/residual.c
./src/c/r_ritz.c
./src/c/sfnorm.c
./src/c/shellsort.c
./src/c/soluf.c
./src/c/sonenorm.c
./src/c/soort.c
./src/c/sort.c
./src/c/sorteig.c
./src/c/tred22.c
./src/c/treesort3.c
./src/c/tresid.c
./src/c/util.c
./src/c/xerbla.c
./src/c/clustrxx5.c
./src/c/inv_it3.c
./src/c/clustrxx4.c
)

SET(PEIGS_F_SOURCE_FILES
./src/c/gdot.F
./src/c/peigs_cmod.f
./src/f77/neblw2.F
./src/f77/dgetavec.F
./src/f77/dgetavec2.f
./src/f77/dgetavec3.f
./src/f77/dlaebz2.F
./src/f77/inverm.F
./src/f77/dspevx2.F
./src/f77/dstebz3.F
./src/f77/damax.f
./src/f77/heapsort.f
./src/f77/types1.f
./src/f77/menode.f
./src/f77/pairup.f  
./src/f77/xerbl2.f  
./src/f77/xstop.f
./comm/mxbrod.F
./comm/mxcombv1.F
./comm/peigs_dgop.F
./comm/maxx.F
./comm/sum.F
./comm/mxsubs.F
./example/gendrv.F
##./src/f77/slaebz2.F
##./src/f77/sstebz3.F
#./src/f77/dgetavec4.F
##./src/f77/sneblw2.F
#./src/f77/slaebz.F
#./new_src/ga_pdspev.F
)

add_library(peigsc STATIC ${PEIGS_C_SOURCE_FILES})
#add_dependencies(peigsc ${NWC_DEPENDENCIES})
target_include_directories(peigsc PUBLIC ${NWC_INCLUDES} ${PEIGS_INCLUDES})
target_compile_definitions(peigsc PRIVATE IEEE)
target_compile_options(peigsc PUBLIC ${NWCHEM_C_COMPILE_FLAGS})
target_compile_options(peigsc PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})

add_library(peigsf STATIC ${PEIGS_F_SOURCE_FILES})
#add_dependencies(peigsf ${NWC_DEPENDENCIES})
target_include_directories(peigsf PUBLIC ${NWC_INCLUDES} ${PEIGS_INCLUDES})
target_compile_options(peigsf PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(peigsf PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

