
message(STATUS "Set compiler flags")

set(CMAKE_C_FLAGS                "")
set(CMAKE_C_FLAGS_DEBUG          "")
set(CMAKE_C_FLAGS_MINSIZEREL     "")
set(CMAKE_C_FLAGS_RELEASE        "")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "")

set(CMAKE_CXX_FLAGS                "")
set(CMAKE_CXX_FLAGS_DEBUG          "")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "")
set(CMAKE_CXX_FLAGS_RELEASE        "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")

set(CMAKE_Fortran_FLAGS                "")
set(CMAKE_Fortran_FLAGS_DEBUG          "")
set(CMAKE_Fortran_FLAGS_MINSIZEREL     "")
set(CMAKE_Fortran_FLAGS_RELEASE        "")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "")


set(NWCHEM_COMPILEDEFS "-DEXT_INT -DLINUX -DLINUX64")
#set(NWCHEM_Fortran_COMPILE_FLAGS "")
#set(NWCHEM_Fortran_OPTIMIZE_FLAGS -O2)
#set(NWCHEM_C_COMPILE_FLAGS "")
#set(NWCHEM_C_OPTIMIZE_FLAGS "")
#set(NWCHEM_LINKER_FLAGS "")
set(NWCHEM_C_DEBUG_FLAGS -g)

set(CPU_X86_64 OFF)
set(CPU_ppc64le OFF)

if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
  set(CPU_X86_64 ON)
elseif("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le")
  set(CPU_ppc64le ON)
else()
  message(FATAL_ERROR "NWChem CMake build does not support ${CMAKE_SYSTEM_PROCESSOR} CPU currently")
endif()

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" AND "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  # Using GNU Compilers

  set(NWCHEM_Fortran_COMPILE_FLAGS ${NWCHEM_Fortran_COMPILE_FLAGS} -m64 -ffast-math -Warray-bounds -fdefault-integer-8)

  set(NWCHEM_Fortran_DEBUG_FLAGS -O2 -g -fno-aggressive-loop-optimizations)
  
  if(CPU_X86_64)
    set(NWCHEM_Fortran_COMPILE_FLAGS ${NWCHEM_Fortran_COMPILE_FLAGS} -finline-functions)
    set(NWCHEM_Fortran_DEBUG_FLAGS ${NWCHEM_Fortran_DEBUG_FLAGS} -g -O)
  endif()

  set(NWCHEM_Fortran_OPTIMIZE_FLAGS ${NWCHEM_Fortran_OPTIMIZE_FLAGS} -O2 -Wuninitialized -fno-aggressive-loop-optimizations)
  
  if(CPU_X86_64)
    set(NWCHEM_Fortran_OPTIMIZE_FLAGS ${NWCHEM_Fortran_OPTIMIZE_FLAGS} -O3 -mfpmath=sse -ffast-math -fprefetch-loop-arrays -ftree-vectorize -fopt-info-vec -mtune=native) 
    if (CMAKE_C_COMPILER_VERSION GREATER_EQUAL 6.0)
      set(NWCHEM_Fortran_OPTIMIZE_FLAGS ${NWCHEM_Fortran_OPTIMIZE_FLAGS} -fno-tree-dominator-opts)
    endif()
  endif()

  set(NWCHEM_C_COMPILE_FLAGS ${NWCHEM_C_COMPILE_FLAGS} -m64 -Wall)

  if(CPU_X86_64)  
    set(NWCHEM_C_COMPILE_FLAGS ${NWCHEM_C_COMPILE_FLAGS} -O3 -funroll-loops -ffast-math)
    set(NWCHEM_C_OPTIMIZE_FLAGS ${NWCHEM_C_OPTIMIZE_FLAGS} -O1)
  elseif(CPU_ppc64le)
    set(NWCHEM_C_COMPILE_FLAGS ${NWCHEM_C_COMPILE_FLAGS} -m64 -O)
  endif()

  set(NWCHEM_COMPILEDEFS "${NWCHEM_COMPILEDEFS} -DGFORTRAN -DCHKUNDFLW -DGCC4 -DGCC46 ")

  if(USE_OPENMP)
    set(NWCHEM_C_COMPILE_FLAGS ${NWCHEM_C_COMPILE_FLAGS} -fopenmp)
    set(NWCHEM_Fortran_COMPILE_FLAGS ${NWCHEM_Fortran_COMPILE_FLAGS} -fopenmp)
    set(NWCHEM_COMPILEDEFS "${NWCHEM_COMPILEDEFS} -DUSE_OPENMP ")
    set(NWCHEM_LINKER_FLAGS "-fopenmp")
  endif()


elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel" AND "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
# Using Intel Compilers

set(NWCHEM_Fortran_COMPILE_FLAGS ${NWCHEM_Fortran_COMPILE_FLAGS} -i8 -align -fpp -qopt-report-file=stder -fimf-arch-consistency=true -finline-limit=250)
set(NWCHEM_Fortran_OPTIMIZE_FLAGS ${NWCHEM_Fortran_OPTIMIZE_FLAGS} -O3 -unroll -ip -xHost)

set(NWCHEM_Fortran_DEBUG_FLAGS -O2 -g -fp-model source)

set(NWCHEM_C_COMPILE_FLAGS ${NWCHEM_C_COMPILE_FLAGS} -xHost -ftz -qopt-report-file=stderr -qopt-report-phase=vec)
set(NWCHEM_C_OPTIMIZE_FLAGS ${NWCHEM_C_OPTIMIZE_FLAGS} -O3 -ip -no-prec-div)

set(NWCHEM_COMPILEDEFS "${NWCHEM_COMPILEDEFS} -DIFCV8 -DIFCLINUX ")

  if(USE_OPENMP)
    set(NWCHEM_C_COMPILE_FLAGS ${NWCHEM_C_COMPILE_FLAGS} -qopenmp -qopt-report-phase:openmp)
    set(NWCHEM_Fortran_COMPILE_FLAGS ${NWCHEM_Fortran_COMPILE_FLAGS} -qno-openmp-offload -qopenmp -qopt-report-phase=openmp)
    set(NWCHEM_COMPILEDEFS "${NWCHEM_COMPILEDEFS} -DUSE_OPENMP ")
    set(NWCHEM_LINKER_FLAGS "-qopenmp")
  else()
    set(NWCHEM_Fortran_COMPILE_FLAGS ${NWCHEM_Fortran_COMPILE_FLAGS} -qno-openmp)
  endif()


elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "XL" AND "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "XL")
# Using IBM XL Compilers

set(NWCHEM_Fortran_COMPILE_FLAGS ${NWCHEM_Fortran_COMPILE_FLAGS} -q64 -qextname -qfixed -NQ40000 -NT80000 -qmaxmem=8192 -qxlf77=leadzero -qintsize=8)
set(NWCHEM_Fortran_OPTIMIZE_FLAGS ${NWCHEM_Fortran_OPTIMIZE_FLAGS} -O3 -qstrict -qcache=auto)

set(NWCHEM_C_COMPILE_FLAGS ${NWCHEM_C_COMPILE_FLAGS}  -q64 -qlanglvl=extended)
#set(NWCHEM_C_OPTIMIZE_FLAGS ${NWCHEM_C_OPTIMIZE_FLAGS} )

set(NWCHEM_Fortran_DEBUG_FLAGS -O2 -g)

set(NWCHEM_COMPILEDEFS "${NWCHEM_COMPILEDEFS} -DXLFLINUX -DCHKUNDFLW ")

else()
  message(FATAL_ERROR "A mix of GNU/Intel/IBM XL C and Fortran Compilers cannot be specified!")
endif()

#message("NWCHEM_Fortran_COMPILE_FLAGS = ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_OPTIMIZE_FLAGS}")

#set(CMAKE_Fortran_FLAGS "${NWCHEM_Fortran_COMPILE_FLAGS}")
#set(CMAKE_C_FLAGS "${NWCHEM_C_COMPILE_FLAGS} ${NWCHEM_C_OPTIMIZE_FLAGS}")

set(NWCHEM_COMPILEDEFS "${NWCHEM_COMPILEDEFS} -DPARALLEL_DIAG ")
add_definitions("${NWCHEM_COMPILEDEFS}")

if(USE_NOIO)
   add_definitions(-DNOIO -DEAFHACK)
endif()

if(USE_DEBUG)
  set(NWCHEM_Fortran_COMPILE_FLAGS ${NWCHEM_Fortran_COMPILE_FLAGS} -g)
  set(NWCHEM_C_COMPILE_FLAGS ${NWCHEM_C_COMPILE_FLAGS} -g)
endif()

if (SCALAPACK_LIBRARIES)
   add_definitions(-DSCALAPACK)
else()
   #Set to empty since this variable is used in nwchem linker
   set(SCALAPACK_LIBRARIES "")
endif()

