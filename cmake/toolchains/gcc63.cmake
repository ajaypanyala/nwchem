
set(CMAKE_BUILD_TYPE Release)

# Compilers (assuming the compilers are in the PATH)
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_Fortran_COMPILER gfortran)
#set(MPI_C_COMPILER mpicc)
#set(MPI_CXX_COMPILER mpicxx)

# MPI 
set(MPI_INCLUDE_PATH /usr/lib/x86_64-linux-gnu/openmpi/include)
set(MPI_LIBRARY_PATH /usr/lib/x86_64-linux-gnu/openmpi/lib)
set(MPI_LIBRARIES "-lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi")

# General NWCHEM Options
set(NWCHEM_MODULES ALL)  
#set(NWCHEM_INSTALL_DIR "")      
#set(ARMCI_NETWORK OPENIB)
option(USE_OPENMP "OpenMP" ON)
option(USE_MPI "Enable MPI" ON)
set(USE_OFFLOAD "OFFLOAD" OFF) 

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
