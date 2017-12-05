cmake_minimum_required(VERSION 3.7.0 FATAL_ERROR)

# NWCHEM superbuild is invoked first. Do not change this option
option (NWCHEM_SUPERBUILD "This will build GA first" ON)

if (NWCHEM_SUPERBUILD)
  project (NWCHEM_SUPERBUILD C CXX Fortran)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  

  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
  include(set_options)
  include(globalarrays)

  #Begin Driver for NWCHEM cmake. 
  #Options,variables set in superbuild through toolchain are preserved. 
  #Variables not set in toolchain file, but set in cmakemodules are lost
  #when cmake is called with superbuild turned OFF. We pass those variables 
  #(ga_root_dir,nwc_install_dir) here again.
  # When using env vars instead of toolchain, all of the options are lost. To 
  # handle env vars case as well and be general, we pass all variables here. 
  # if an env var is not set, its value is empty ie OFF. 
  # NOTE: call set_options again to handle both cases - this does not seem like a clean
  # solution and also does not work as expected.
  include(ExternalProject)
    ExternalProject_Add (NWCHEM_SB
    DEPENDS GLOBALARRAYS
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}
    #Only this line is needed when building with toolchain - print all vars and double check.
    #CMAKE_ARGS -DNWCHEM_SUPERBUILD=OFF -DGA_ROOT_DIR=${GA_ROOT_DIR} -DCMAKE_INSTALL_PREFIX=${NWCHEM_INSTALL_DIR}
    #The following command is not needed if everyone agrees to use toolchain files instead of env vars
    CMAKE_ARGS -DNWCHEM_SUPERBUILD=OFF -DGA_ROOT_DIR=${GA_ROOT_DIR} -DCMAKE_INSTALL_PREFIX=${NWCHEM_INSTALL_DIR}
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DMPI_INCLUDE_PATH=${MPI_INCLUDE_PATH} 
    -DMPI_LIBRARY_PATH=${MPI_LIBRARY_PATH} -DMPI_LIBRARIES=${MPI_LIBRARIES}
    -DBLAS_LIBRARIES=${BLAS_LIBRARIES} -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES} 
    -DSCALAPACK_LIBRARIES=${SCALAPACK_LIBRARIES} 
    -DNWCHEM_MODULES=${NWCHEM_MODULES} -DUSE_OPENMP=${USE_OPENMP} -DUSE_MPI=${USE_MPI}
    -DUSE_OFFLOAD=${USE_OFFLOAD} -DUSE_MPIF=${USE_MPIF} -DUSE_MPIF4=${USE_MPIF4}
    -DUSE_MLIB=${USE_MLIB} -DMA_TRANS_BLOCKED=${MA_TRANS_BLOCKED} -DIPCCSD=${IPCCSD}
    -DEACCSD=${EACCSD} -DMRCC_METHODS=${MRCC_METHODS} -DCCSDTLR=${CCSDTLR} 
    -DCCSDTQ=${CCSDTQ} -DTCE_CUDA=${TCE_CUDA} -DCUDA_INCLUDE=${CUDA_INCLUDE}
    -DUSE_PSTAT=${USE_PSTAT} -DTCE_NEW_OPENMP=${TCE_NEW_OPENMP} -DOPENMP_OFFLOAD=${OPENMP_OFFLOAD}
    INSTALL_COMMAND make install
    )
   # End Driver for NWCHEM cmake

  return() # stop further processing
else()
  project (NWCHEM C CXX Fortran)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  
  include(set_compile_flags)

    execute_process(
      COMMAND ./ga-config --libs
      WORKING_DIRECTORY ${GA_ROOT_DIR}/bin
      OUTPUT_VARIABLE GA_EXP_LIBS 
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
endif()

#message(STATUS "GA_EXP_LIBS= " ${GA_EXP_LIBS})
set(GA_LIBS ${GA_ROOT_DIR}/lib/libga.a)

if(EXISTS ${GA_LIBS})
    message(STATUS "Found Global Arrays at: ${GA_ROOT_DIR}")
else()
    #This should never happen
    message(FATAL_ERROR "GA not found. Something went wrong with the build. Please run cmake again.")
endif()


set(NWC_INCLUDES
${MPI_INCLUDE_PATH}
${GA_ROOT_DIR}/include
${PROJECT_SOURCE_DIR}/src
${PROJECT_SOURCE_DIR}/src/include
${PROJECT_SOURCE_DIR}/src/util
#${PROJECT_SOURCE_DIR}/src/config
${PROJECT_SOURCE_DIR}/src/basis
${PROJECT_SOURCE_DIR}/src/geom
${PROJECT_SOURCE_DIR}/src/symmetry
${PROJECT_SOURCE_DIR}/src/pstat
${PROJECT_SOURCE_DIR}/src/property
${PROJECT_SOURCE_DIR}/src/ddscf
${PROJECT_SOURCE_DIR}/src/bq
${PROJECT_SOURCE_DIR}/src/inp
${PROJECT_SOURCE_DIR}/src/rtdb
${PROJECT_SOURCE_DIR}/src/NWints/api
)

set(CMAKE_C_FLAGS_RELEASE -DNDEBUG)
set(CMAKE_CXX_FLAGS_RELEASE -DNDEBUG)
set(CMAKE_Fortran_FLAGS_RELEASE -O2) 
#message("CXX RELEASE = ${CMAKE_C_FLAGS_RELEASE}") 
#message("FC RELEASE = ${CMAKE_Fortran_FLAGS_RELEASE}")

add_subdirectory(src)
add_definitions(-DCOMPILATION_DATE="04-24-2017")
add_definitions(-DCOMPILATION_DIR=${PROJECT_SOURCE_DIR})
add_definitions(-DNWCHEM_BRANCH="RELEASE")

set(LIBTCE_ALL tce ${LIBTCE_OPT} tce_lrh tce_ccsdt tce_eomccsd_act 
tce_ccsd tce_mbpt4 tce_emb tce_cisdt 
 tce_cisd  tce_gradients tce_mbpt2 
tce_ccd tce_diis3 tce_cis  tce_ccsd_lambda tce_cisdtq tce_ccsd_qr 
tce_mbpt3 tce_eomccsd tce_cr_eomccsd_t_act 
tce_lccd tce_cc3 tce_qcisd  tce_diis tce_ccsdt_lambda tce_mbpt1 
tce_eomccsdt tce_diis2 tce_cr-eomccsd_t tce_lccsd tce_ccsd_t 
tce_ccsd_act tce_ccsdt tce tce_sort
#These are compiled even if MRCC_METHODS is not set
tce_mrcc_ccsd tce_mrcc_uss tce_mrcc_ccsd_sub tce_mrcc_ccsd_t)


set(NWC_SmallQM_LIBS NWints atomscf ddscf gradients moints nwdft rimp2 hessian bq cons optim
  stepper driver vib qmd mcscf property solvation selci dplot mp2_grad esp etrans vscf nwmd dimqm)

set(NWC_NWPW_LIBS nwpw paw band nwpwlib pspw)
set(ALL_NWC_LIBS nwc_task nwc_input nwc_task ccsd mcscf 
atomscf selci mp2_grad moints stepper driver optim nwdft gradients cphf geom 
dimqm drdy driver qmd mm bq mcscf dplot property ccca ${NWC_NWPW_LIBS} esp ddscf dangchang hessian 
vib nwcutil symmetry rimp2 property solvation basis NWints prepar nwmd ${NWC_NWPW_LIBS} cafe 
space analyz qhop pfft dplot drdy vscf qmmm ${NWC_NWPW_LIBS} qmd etrans bq mm cons perfm dntmc ccca  
dimqm drdy stepper qmd prepar ccsd esp dplot ccca qmmm nwcutil symmetry rdmft)

set(OPT_LIBS ${NETLIB_BLAS} ncc chelp cckohn)

#link_directories(${GA_ROOT_DIR}/lib)
#link_directories(${MPI_LIBRARY_PATH})
add_executable(nwchem ${PROJECT_SOURCE_DIR}/src/nwchem.F ${PROJECT_SOURCE_DIR}/src/stubs_cmake.F) 
set_target_properties(nwchem PROPERTIES LINKER_LANGUAGE Fortran)
target_include_directories(nwchem PUBLIC ${NWC_INCLUDES})
add_dependencies(nwchem ${LIBTCE_ALL} ${ALL_NWC_LIBS} ${OPT_LIBS})

#Internal use only.
if (NETLIB_BLAS)
  set(TAMM_LIBRARIES_OPTIONAL -lmpi_cxx)
endif()
#set(TAMM_LIBRARIES -L(TAMM_LIB_PATH) -ltamm -lstdc++ (TAMM_LIBS) ${TAMM_LIBRARIES_OPTIONAL})

if(NOT "${CMAKE_HOST_SYSTEM_NAME}" STREQUAL "Darwin")
  set(NWCHEM_LINKER_FLAGS "${NWCHEM_LINKER_FLAGS} -lrt")
endif()

#message("GA_LIBS = ${GA_EXP_LIBS}")
target_link_libraries(nwchem ${NWCHEM_Fortran_COMPILE_FLAGS} 
-L${GA_ROOT_DIR}/lib -L${MPI_LIBRARY_PATH} ${GA_LIBS}
 ${ALL_NWC_LIBS} ga ${LIBTCE_ALL} geom armci ga peigsc peigsf 
cons bq nwc_input symmetry nwcutil nwc_inp rtdb pstat perfm nwcutil  
rtdb ga ${TAMM_LIBRARIES} ${LIBTCE_ALL} ${NWCHEM_LINKER_FLAGS}  ${TCE_CUDA_LIBS}
${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} 
${MPI_LIBRARIES} ${GA_EXP_LIBS})

install(TARGETS nwchem DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

#install(TARGETS ${ALL_NWC_LIBS} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)


