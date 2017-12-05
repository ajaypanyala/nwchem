
set(CMAKE_BUILD_TYPE Release)

# Compilers (assuming the compilers are in the PATH)
set(CMAKE_C_COMPILER icc)
set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_Fortran_COMPILER ifort)
#set(MPI_C_COMPILER mpicc)
#set(MPI_CXX_COMPILER mpicxx)

# MPI 
set(MPI_INCLUDE_PATH /usr/lib/openmpi/include/)
set(MPI_LIBRARY_PATH /usr/lib/openmpi/lib/)
set(MPI_LIBRARIES "-lmpi_f77 -lmpi -ldl -lhwloc")

# BLAS, LAPACK & SCALAPACK. Assume 8-byte for now.
set(BLAS_LIBRARIES "-mkl -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl" CACHE STRING "BLAS linker flags")
set(LAPACK_LIBRARIES "${BLAS_LIBRARIES}" CACHE STRING "LAPACK linker flags")
set(SCALAPACK_LIBRARIES "-mkl -lmkl_scalapack_ilp64 -lmkl_blacs_openmpi_ilp64 -lmkl_intel_thread -lpthread -lm -ldl" CACHE STRING "SCALAPACK linker flags")


# General NWCHEM Options
set(NWCHEM_MODULES ALL)  
#set(NWCHEM_INSTALL_DIR "")      
#set(ARMCI_NETWORK OPENIB)
option(USE_OPENMP "OpenMP" ON)
option(USE_MPI "Enable MPI" ON)
set(USE_OFFLOAD "OFFLOAD" OFF) #Intel only ?

#set(ARMCI_DEFAULT_SHMMAX_UBOUND 65536) #Always set to 131072 in tools/GNUMakefile
#set(INTEL_64ALIGN 1)  #Set when USE_OFFLOAD is set
#set(BLAS_SIZE 8)
#USE_NOIO
#export USE_CPPRESERVE=y
#export USE_NOFSCHECK=y

# NWPW Options
option(USE_MPIF "MPIF" ON)
option(USE_MPIF4 "MPIF4" ON)
option(USE_MLIB "USE_MLIB" OFF)

# TCE Options
option(MA_TRANS_BLOCKED "MA Blocked" OFF)
option(IPCCSD "Build IPCCSD" ON)
option(EACCSD "Build EACCSD" ON)
option(MRCC_METHODS "Build MRCC" ON)	
option(CCSDTLR "CCSDTLR" OFF)
option(CCSDTQ "CCSDTQ" OFF)
option(TCE_CUDA "TCE CUDA" OFF)
set(GPU_ARCH sm_35)
option(CCSD_T_GPU_OLD "CCSD_T w/o ttlg" OFF)
# set(CUDA_INCLUDE "")
option(USE_PSTAT "PSTAT" OFF)
option(TCE_NEW_OPENMP "TCE_NEW_OPENMP" OFF)
option(OPENMP_OFFLOAD "OPENMP_OFFLOAD" OFF)
