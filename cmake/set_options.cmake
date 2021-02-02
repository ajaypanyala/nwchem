
if(NOT "${CMAKE_HOST_SYSTEM_NAME}" STREQUAL "Linux" AND NOT "${CMAKE_HOST_SYSTEM_NAME}" STREQUAL "Darwin")
  message(FATAL_ERROR "NWChem cmake build only supports Linux and MacOSX as of now.")
endif()
if(NOT "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" AND NOT "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" AND NOT "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "XL")
  message(FATAL_ERROR "NWChem cmake build only supports building with GNU, Intel and IBM XL compilers only as of now.")
endif()


if(DEFINED NWCHEM_INSTALL_DIR)
    set(NWCHEM_INSTALL_DIR ${NWCHEM_INSTALL_DIR})
else()
    set(NWCHEM_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/nwchem-install)
endif()

if(NOT USE_MPI)
    message(FATAL_ERROR "Must specify USE_MPI and have a working MPI installation to compile NWCHEM")
else()
    set(USE_MPI ON)
endif()

# NWPW
if (USE_MPIF)
   set(USE_MPIF ON)
endif()

if (USE_MPIF4)
   set(USE_MPIF4 ON)
endif()

if (USE_MLIB)
   set(USE_MLIB ON)
endif()


# TCE Options

if(MA_TRANS_BLOCKED)
   set(MA_TRANS_BLOCKED ON)
endif()

if(IPCCSD)
   set(IPCCSD ON)
endif()

if(EACCSD)
   set(EACCSD ON)
endif()

if(CCSDTLR)
   set(CCSDTLR ON)
endif()

if(CCSDTQ)
   set(CCSDTQ ON)
endif()

if(MRCC_METHODS)
   set(MRCC_METHODS ON)
endif()

#if(DEFINED ENV{USE_PSTAT} OR USE_PSTAT)
if(USE_PSTAT)
   set(USE_PSTAT ON)
endif()

#if(DEFINED ENV{TCE_NEW_OPENMP} OR TCE_NEW_OPENMP)
if(TCE_NEW_OPENMP)
   set(TCE_NEW_OPENMP ON)
endif()

#if(DEFINED ENV{OPENMP_OFFLOAD} OR OPENMP_OFFLOAD)
if(OPENMP_OFFLOAD)
   set(OPENMP_OFFLOAD ON)
endif()

if(TCE_CUDA)
   set(TCE_CUDA ON)
   INCLUDE(FindCUDA)
   if(CUDA_TOOLKIT_INCLUDE)
        set(CUDA_TOOLKIT_INCLUDE ${CUDA_TOOLKIT_INCLUDE})
    else()
        message(WARNING "TCE_CUDA option is enabled, but CUDA_TOOLKIT_INCLUDE is not set.")
    endif()
endif()

# General NWCHEM Options
if(USE_OFFLOAD)
   set(USE_OFFLOAD ON)
endif()

if(USE_OPENMP)
   set(USE_OPENMP ON)
endif()

if(NWCHEM_MODULES)
    set(NWCHEM_MODULES ${NWCHEM_MODULES})
# elseif(DEFINED ENV{NWCHEM_MODULES})
#     set(NWCHEM_MODULES $ENV{NWCHEM_MODULES})
endif()

if(ARMCI_NETWORK)
    set(ARMCI_NETWORK ${ARMCI_NETWORK})
# elseif(DEFINED ENV{ARMCI_NETWORK})
#     set(ARMCI_NETWORK $ENV{ARMCI_NETWORK})
endif()

if(MPI_INCLUDE_PATH)
    set(MPI_INCLUDE_PATH ${MPI_INCLUDE_PATH})
# elseif(DEFINED ENV{MPI_INCLUDE_PATH})
#     set(MPI_INCLUDE_PATH $ENV{MPI_INCLUDE_PATH})
endif()

if(MPI_LIBRARY_PATH)
    set(MPI_LIBRARY_PATH ${MPI_LIBRARY_PATH})
# elseif(DEFINED ENV{MPI_LIBRARY_PATH})
#     set(MPI_LIBRARY_PATH $ENV{MPI_LIBRARY_PATH})
endif()

if(MPI_LIBRARIES)
    set(MPI_LIBRARIES ${MPI_LIBRARIES})
# elseif(DEFINED ENV{MPI_LIBRARIES})
#     set(MPI_LIBRARIES $ENV{MPI_LIBRARIES})
endif()

if(BLAS_LIBRARIES)
    set(BLAS_LIBRARIES ${BLAS_LIBRARIES})
# elseif(DEFINED ENV{BLAS_LIBRARIES})
#     set(BLAS_LIBRARIES $ENV{BLAS_LIBRARIES})
endif()

if(LAPACK_LIBRARIES)
    set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES})
# elseif(DEFINED ENV{LAPACK_LIBRARIES})
#     set(LAPACK_LIBRARIES $ENV{LAPACK_LIBRARIES})
endif()

if(SCALAPACK_LIBRARIES)
    set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARIES})
# else
# if(DEFINED ENV{SCALAPACK_LIBRARIES})
#     set(SCALAPACK_LIBRARIES $ENV{SCALAPACK_LIBRARIES})
endif()

# if(DEFINED ENV{BLAS_SIZE})
#    set(BLAS_SIZE $ENV{BLAS_SIZE})
# else
if(BLAS_SIZE)
   set(BLAS_SIZE ${BLAS_SIZE})
endif()


