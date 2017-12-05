
cmake >= 3.8

Only LINUX64 target - gcc, intel, xl compilers - c,c++,fortran compilers should all be either GNU or Intel or XL - cannot mix

GCC > 4.8.2 and gcc LT 7 and intel fortran/c >=16

if using cuda - tensorgen lib needs c++11 - only gcc 4.8 - 5.4 - use cuda8 - do not use gcc 6 since cuda 8 does not support it

USE_OFFLOAD = OFF (USE_OFFLOAD - knl not supported)

USE_DEBUG support added - add to toolchain file

-qsave added from somewhere ?

ddscf_unopt, nwdft/rt_tddft, hessian, cphf -> -ffixed-line-length-none

CUDA_INCLUDE="-I/share/apps/cuda/8.0.44/include"
CUDA_LIBS="/share/apps/cuda/8.0.44/lib64/libcublas.so /share/apps/cuda/8.0.44/lib64/libcudart.so"
EXTRA_LIBS="src/tce/ttlg/libttlg.a -lstdc++"


Handle
-------

nwdft - scf_dft, so_dft, spec



Fortran

OPT -> band_energy.F task_band.F
---------------------------------

ifort  -c -i8 -align -fpp -qopt-report-file=stderr -qopenmp -qopt-report-phase=openmp -qno-openmp-offload -fimf-arch-consistency=true -finline-limit=250 -O3  -unroll  -ip -xHost -I.  -I/home/panyala/software/nwchem-devel/src/include -I/home/panyala/software/nwchem-devel/src/tools/install/include -DEXT_INT -DLINUX -DLINUX64 -DUSE_OPENMP  -DIFCV8 -DIFCLINUX -DSCALAPACK -DNOIO -DEAFHACK -DPARALLEL_DIAG   band_energy.F

gfortran  -c -m64 -ffast-math  -Warray-bounds -fopenmp -fdefault-integer-8 -finline-functions -O2  -Wuninitialized -fno-aggressive-loop-optimizations -O3  -mfpmath=sse  -fno-tr\
ee-dominator-opts  -ffast-math  -fprefetch-loop-arrays  -ftree-vectorize   -fopt-info-vec -mtune=native -I.  -I/home/panyala/software/nwchem-devel/src/include -I/home/panyala/s\
oftware/nwchem-devel/src/tools/install/include -DGFORTRAN -DCHKUNDFLW -DGCC4 -DGCC46 -DUSE_OPENMP -DEXT_INT -DLINUX -DLINUX64 -DPARALLEL_DIAG   task_band.F

IBM XL
-------
gfortran  -c -m64 -ffast-math  -Warray-bounds -fopenmp -fdefault-integer-8 -O2  -Wuninitialized -fno-aggressive-loop-optimizations -I.  -I/ccs/home/panyala/code/nwchem-devel/src/include -I/ccs/home/panyala/code/nwchem-devel/src/tools/install/include -DGFORTRAN -DCHKUNDFLW -DGCC4 -DGCC46 -DUSE_OPENMP -DEXT_INT -DLINUX -DLINUX64 -DPARALLEL_DIAG   task_band.F

xlf  -c -q64 -qextname -qfixed  -NQ40000 -NT80000 -qmaxmem=8192 -qxlf77=leadzero -qintsize=8 -qsave -O3 -qstrict -qcache=auto   -I.  -I/ccs/home/ngawande/Summit/nwchem/nwchem-regb/src/include -I/ccs/home/ngawande/Summit/nwchem/nwchem-regb/src/tools/install/include -WF,"-DEXT_INT -DLINUX -DLINUX64 -DXLFLINUX -DCHKUNDFLW -DPARALLEL_DIAG"   geom_hnd.F


UNOPT -> band_input.F psp_generator_input.F
---------------------------------

ifort  -c -i8 -align -fpp -qopt-report-file=stderr -qopenmp -qopt-report-phase=openmp -qno-openmp-offload -fimf-arch-consistency=true -finline-limit=250 -O2 -g -fp-model source  -I.  -I/home/panyala/software/nwchem-devel/src/include -I/home/panyala/software/nwchem-devel/src/tools/install/include -DEXT_INT -DLINUX -DLINUX64 -DUSE_OPENMP  -DIFCV8 -DIFCLINUX -DSCALAPACK -DNOIO -DEAFHACK -DPARALLEL_DIAG   band_input.F

gfortran  -c -m64 -ffast-math  -Warray-bounds -fopenmp -fdefault-integer-8 -finline-functions -O2 -g -fno-aggressive-loop-optimizations -g -O   -I.  -I/home/panyala/software/nw\
chem-devel/src/include -I/home/panyala/software/nwchem-devel/src/tools/install/include -DGFORTRAN -DCHKUNDFLW -DGCC4 -DGCC46 -DUSE_OPENMP -DEXT_INT -DLINUX -DLINUX64 -DPARALLEL\
_DIAG   band_input.F

IBM XL
-------

gfortran  -c -m64 -ffast-math  -Warray-bounds -fopenmp -fdefault-integer-8 -O2 -g -fno-aggressive-loop-optimizations  -I.  -I/ccs/home/panyala/code/nwchem-devel/src/include -I/ccs/home/panyala/code/nwchem-devel/src/tools/install/include -DGFORTRAN -DCHKUNDFLW -DGCC4 -DGCC46 -DUSE_OPENMP -DEXT_INT -DLINUX -DLINUX64 -DPARALLEL_DIAG   band_input.F

xlf  -c -q64 -qextname -qfixed  -NQ40000 -NT80000 -qmaxmem=8192 -qxlf77=leadzero -qintsize=8 -O2 -g  -I.  -I/ccs/home/ngawande/Summit/nwchem/nwchem-regb/src/include -I/ccs/home/ngawande/Summit/nwchem/nwchem-regb/src/tools/install/include -WF,"-DEXT_INT -DLINUX -DLINUX64 -DXLFLINUX -DCHKUNDFLW -DPARALLEL_DIAG"   psp_generator_input.F



----
C
----

OPT-> paw_atom.c, paw_ion.c
---------------------------------

icc -c -I.  -I/home/panyala/software/nwchem-devel/src/include -I/home/panyala/software/nwchem-devel/src/tools/install/include -DEXT_INT -DLINUX -DLINUX64 -DUSE_OPENMP  -DIFCV8 -DIFCLINUX -DSCALAPACK -DNOIO -DEAFHACK -DPARALLEL_DIAG  -xHost -ftz -qopt-report-phase=vec  -qopt-report-file=stderr -qopenmp -qopt-report-phase:openmp -O3 -ip -no-prec-div -o paw_ion.o paw_ion.c

gcc -c -I.  -I/home/panyala/software/nwchem-devel/src/include -I/home/panyala/software/nwchem-devel/src/tools/install/include -DGFORTRAN -DCHKUNDFLW -DGCC4 -DGCC46 -DUSE_OPENMP\
 -DEXT_INT -DLINUX -DLINUX64 -DPARALLEL_DIAG  -m64 -Wall -O3 -funroll-loops -ffast-math -fopenmp -O1 -o paw_ion.o paw_ion.c

IBM XL
-------
xlc -c -I.  -I/ccs/home/ngawande/Summit/nwchem/nwchem-regb/src/include -I/ccs/home/ngawande/Summit/nwchem/nwchem-regb/src/tools/install/include -DEXT_INT -DLINUX -DLINUX64 -DXLFLINUX -DCHKUNDFLW -DPARALLEL_DIAG  -q64 -qlanglvl=extended  -o paw_atom.o paw_atom.c

gcc -c -I.  -I/ccs/home/panyala/code/nwchem-devel/src/include -I/ccs/home/panyala/code/nwchem-devel/src/tools/install/include -DGFORTRAN -DCHKUNDFLW -DGCC4 -DGCC46 -DUSE_OPENMP -DEXT_INT -DLINUX -DLINUX64 -DPARALLEL_DIAG  -m64 -Wall -m64  -O  -o paw_ion.o paw_ion.c




UNOPT -> dirac.c, pspsolve.c
---------------------------------

icc -c -I.  -I/home/panyala/software/nwchem-devel/src/include -I/home/panyala/software/nwchem-devel/src/tools/install/include -DEXT_INT -DLINUX -DLINUX64 -DUSE_OPENMP  -DIFCV8 -DIFCLINUX -DSCALAPACK -DNOIO -DEAFHACK -DPARALLEL_DIAG  -xHost -ftz -qopt-report-phase=vec  -qopt-report-file=stderr -qopenmp -qopt-report-phase:openmp -g  -o dirac.o dirac.c


 gcc -c -I.  -I/home/panyala/software/nwchem-devel/src/include -I/home/panyala/software/nwchem-devel/src/tools/install/include -DGFORTRAN -DCHKUNDFLW -DGCC4 -DGCC46 -DUSE_OPENMP\
 -DEXT_INT -DLINUX -DLINUX64 -DPARALLEL_DIAG  -m64 -Wall -O3 -funroll-loops -ffast-math -fopenmp -g  -o dirac.o dirac.c

IBM XL
-------

xlc -c -I.  -I/ccs/home/ngawande/Summit/nwchem/nwchem-regb/src/include -I/ccs/home/ngawande/Summit/nwchem/nwchem-regb/src/tools/install/include -DEXT_INT -DLINUX -DLINUX64 -DXLFLINUX -DCHKUNDFLW -DPARALLEL_DIAG  -q64 -qlanglvl=extended -g  -o dirac.o dirac.c

gcc -c -I.  -I/ccs/home/panyala/code/nwchem-devel/src/include -I/ccs/home/panyala/code/nwchem-devel/src/tools/install/include -DGFORTRAN -DCHKUNDFLW -DGCC4 -DGCC46 -DUSE_OPENMP -DEXT_INT -DLINUX -DLINUX64 -DPARALLEL_DIAG  -m64 -Wall -m64  -O -g  -o dirac.o dirac.c



---

Install XCode Dev Tools to get clang

xcode-select --install
 
Install brew from https://brew.sh/
 
brew install lzlib wget flex bison doxygen autoconf automake libtool
 
Get cmake from  https://cmake.org/files/v3.9/cmake-3.9.1.tar.gz
 
CC=clang CXX=clang++ ./configure --prefix=/path-to-install-cmake
 
Then install GCC + openmpi using the following
 
wget http://mirrors-usa.go-parts.com/gcc/releases/gcc-6.3.0/gcc-6.3.0.tar.gz
tar xf gcc-6.3.0.tar.gz
cd gcc-6.3.0
./contrib/download_prerequisites
./configure --prefix=/home/user/gcc-6.3 --disable-multilib --enable-languages=c,c++,fortran
make -j16
make install
cd ..
 
 
wget https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.0.tar.gz
tar xf openmpi-2.1.0.tar.gz
echo "Building depend mpi"
cd openmpi-2.1.0
./configure --prefix=$1/openmpi-2.1 --enable-mpi-cxx --enable-mpi-fortran
make -j 16
make install
cd ../
 



---
on summitdev replace $NWCHEM_TOP/src/nwpw/nwpwlib/md_xs/FEFF6L/feff6Lib.f with feff6Lib.f-summitdev

summitdev ibm xl

config/makefile.h

line 2100 - 2158 handles the compile options for LINUX64 target using ibm xl compilers.
on summitdev - xl compilers have to be loaded and used for building NWCHEM. ga will not build correctly when using gcc on summitdev


?       stubs_cm.F  (contains for ALL - generated stubs.F, util/util_nwchem_version.F, uitl_ga_version.F)
M       tce/mbpt2/mbpt2_t4.F
?       tools/autotools
?       tools/ga-5-6
?       util/util_module_avail.F -- copied code to stubs_cm.F

fft/ - macros not used
peigs - opt or not ?

split cmake files
-----------------
nwdft - did not handle macros defined using nm on object files.
NWints - done
nwpw
tce



NOTES
-=====

Did not handle all cases - only tce
- nwpw cmakelists for example
- obj= obj_optimize= in gnumake file - compile with and w/o -Ox flags - i assume -Ox on everywhere (disable opt on ex ibm for some files ibm xl compiler gives internal compile error)
- USES_BLAS = not doing this either

- Setup autotools build - copy setup from GA
- Assume Intel compilers > v15


find . -name CMakeLists.txt -exec sed -i "s/${NWC_LIBRARIES})/)/g" {} \;

__LP64__
_LP64
These macros are defined, with value 1, if (and only if) the compilation is for a target where long int and pointer both use 64-bits and int uses 32-bit.

-m32
-m64
Generate code for a 32-bit or 64-bit environment. The 32-bit environment sets int, long and pointer to 32 bits and generates code that runs on any i386 system. The 64-bit environment sets int to 32 bits and long and pointer to 64 bits and generates code for AMD's x86-64 architecture. For darwin only the -m64 option turns off the -fno-pic and -mdynamic-no-pic options.


make nwchem_config generates:

nwchem-6.6/src/config: NWCHEM_CONFIG
nwchem-6.6/src/config: nwchem_config.h
nwchem-6.6/src: stubs.F
nwchem-6.6/src/util: util_module_avail.F


THe main makefile in src/

imports config/makefile.h - which in turn include nwchem_config.h and defines NWSUBDIRS = NW_CORE_SUBDIRS and NW_MODULE_SUBDIRS (nwchem_config.h)

NWSUBDIRS = $(NW_CORE_SUBDIRS) $(NW_MODULE_SUBDIRS)

NW_CORE_SUBDIRS = tools include basis geom inp input  \
	pstat rtdb task symmetry util peigs perfm bq cons $(CORE_SUBDIRS_EXTRA = blas lapack )

NW_MODULE_SUBDIRS = NWints atomscf ddscf gradients moints nwdft rimp2 stepper driver optim cphf ccsd vib mcscf prepar esp hessian selci dplot mp2_grad qhop property solvation nwpw fft python analyz nwmd cafe space drdy vscf qmmm qmd etrans tce bq mm cons perfm dntmc dangchang ccca

CORE_LIBS = -lnwcutil -lga -larmci -lpeigs -lperfm -lcons -lbq -lnwcutil

NW_MODULE_LIBS = -lccsd -lmcscf -lselci -lmp2 -lmoints -lstepper -ldriver -loptim -lnwdft -lgradients -lcphf -lesp -lddscf -ldangchang -lguess -lhessian -lvib -lnwcutil -lrimp2 -lproperty -lsolvation -lnwints -lprepar -lnwmd -lnwpw -lofpw -lpaw -lpspw -lband -lnwpwlib -lcafe -lspace -lanalyze -lqhop -lpfft -ldplot -lnwpython -ldrdy -lvscf -lqmmm -lqmd -letrans -lpspw -ltce -lbq -lmm -lcons -lperfm -ldntmc -lccca

-DNWX_LIBRARIES=/opt/libraries/nwchem-6.6/src/tools/install/lib/libga.a

70 modules containing fortran sources, except these 14, the remaining 55 modules were built in 15min
'nwxc','uccsdt','lucia','ccsd','property','qhop','symmetry','nwpw','ddscf','peigs','NWints','selci','rimp2','smd'




----------------------------------------


PROBLEMATIC SOURCES:
--------------------

What about RDMFT - hangs ?

COMPILE WITH FIXES
------------------
tce/mrcc/tce_mrcc_energy.F - compiles if you set macro -DMRCC_METHODS in mrcc CMakeLists.txt (export MRCC_METHODS=y in bashrc)
ddscf-fast- fastj.F,xlm_poles.F,newfmm.F,solver.F -> add defs in cmakefile -DFMM_LMAX=30 -compiles

peigs/CMAKELISTS - DTCGMSG had to be defined else llog() method will be redfined (libga) during link time


nwpwlib/Parallel - disable parallel-tcgmsg.F for now
ifeq ($(USE_MPIF),y) Compile Parallel-mpi.o
else  OBJ_OPTIMIZE += Parallel-tcgmsg.o


REMOVE/ALTER
-------------
util/util_nwchem_version.F
util/util_ga_version.F
util/util_module_avail.F -- Manually copied from built nwchem version to util folder

are generated and called in nwchem.F - hard code version numbers for now

BUGS
====

modify and compile ../src/nwxc/nwxc_eval.F to meet calls to nwxc_eval_df() [that are invoked if nwxc is used] in if-else conditions from everywhere in nwchem
In input/input_parse.F and task -> task.F, task_energy.F, task_dynamics.F - comment calls to unused modules like argos, lucia, smd, diana, etc

./nwpw/nwpwlib/io/second.F - #define TCGMSG forcibly at top causing redeclaration if USE_MPI macro is defined.

NOT COMPILING
-------------
tce/crop/tce_crop.F - shift error - The shapes of the array expressions do not conform. 135,138,391,394

not linked = but used - input rtdb pstat geom


----------------------------------

ifort -i8 -align -fpp -qopt-report-file=stderr -qopenmp -qopt-report-phase=openmp -qno-openmp-offload -fimf-arch-consistency=true -finline-limit=250 -O2 -g -fp-model source  -I.  -I/opt/libraries/nwchem-6.6/src/include -I/opt/libraries/nwchem-6.6/src/tools/install/include -DEXT_INT -DLINUX -DLINUX64 -DUSE_OPENMP  -DIFCV8 -DIFCLINUX -DSCALAPACK -DPARALLEL_DIAG -DCOMPILATION_DATE="'`date +%a_%b_%d_%H:%M:%S_%Y`'" -DCOMPILATION_DIR="'/opt/libraries/nwchem-6.6'" -DNWCHEM_BRANCH="'6.6'"  -c -o nwchem.o nwchem.F
ifort -i8 -align -fpp -qopt-report-file=stderr -qopenmp -qopt-report-phase=openmp -qno-openmp-offload -fimf-arch-consistency=true -finline-limit=250 -O2 -g -fp-model source  -I.  -I/opt/libraries/nwchem-6.6/src/include -I/opt/libraries/nwchem-6.6/src/tools/install/include -DEXT_INT -DLINUX -DLINUX64 -DUSE_OPENMP  -DIFCV8 -DIFCLINUX -DSCALAPACK -DPARALLEL_DIAG -DCOMPILATION_DATE="'`date +%a_%b_%d_%H:%M:%S_%Y`'" -DCOMPILATION_DIR="'/opt/libraries/nwchem-6.6'" -DNWCHEM_BRANCH="'6.6'"  -c -o stubs.o stubs.F
make[1]: Leaving directory '/opt/libraries/nwchem-6.6/src'
ifort -i8 -align -fpp -qopt-report-file=stderr -qopenmp -qopt-report-phase=openmp -qno-openmp-offload -fimf-arch-consistency=true -finline-limit=250 -O2 -g -fp-model source   -L/opt/libraries/nwchem-6.6/lib/LINUX64 -L/opt/libraries/nwchem-6.6/src/tools/install/lib  -o /opt/libraries/nwchem-6.6/bin/LINUX64/nwchem nwchem.o stubs.o -lnwctask -lccsd -lmcscf -lselci -lmp2 -lmoints -lstepper -ldriver -loptim -lnwdft -lgradients -lcphf -lesp -lddscf -ldangchang -lguess -lhessian -lvib -lnwcutil -lrimp2 -lproperty -lsolvation -lnwints -lprepar -lnwmd -lnwpw -lofpw -lpaw -lpspw -lband -lnwpwlib -lcafe -lspace -lanalyze -lqhop -lpfft -ldplot -ldrdy -lvscf -lqmmm -lqmd -letrans -lpspw -ltce -lbq -lmm -lcons -lperfm -ldntmc -lccca -lnwcutil -lga -larmci -lpeigs -lperfm -lcons -lbq -lnwcutil   -mkl -lmkl_scalapack_ilp64 -lmkl_blacs_intelmpi_ilp64 -lmkl_intel_thread -lpthread -lm   -mkl -lmkl_lapack95_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm   -L/usr//lib -L/usr/lib/openmpi/lib -lmpi_f90 -lmpi_f77 -lmpi -ldl -lhwloc   -libverbs   -lrt  -lpthread
/bin/rm -f nwchem.o stubs.o


---

ifort  -c -i8 -align -fpp -qopt-report-file=stderr -qopenmp -qopt-report-phase=openmp -qno-openmp-offload -fimf-arch-consistency=true -finline-limit=250 -O2 -g -fp-model source  -DEXT_INT -DLINUX -DLINUX64 -DUSE_OPENMP  -DIFCV8 -DIFCLINUX -DSCALAPACK -DPARALLEL_DIAG   sym_mol.F


TODO
====

src/ccsd module
flags and other optimization options

NWCHEM MODULES - try smallQM first
tce-cuda/KNL
python stuff
set use_offload to INTEL_64ALIGN=1 only when using intel compiler ?

#set(ARMCI_DEFAULT_SHMMAX_UBOUND 65536) #Always set to 131072 in tools/GNUMakefile
#set(INTEL_64ALIGN 1)  #Set when USE_OFFLOAD is set
#set(BLAS_SIZE 8)
#USE_NOIO
#export USE_CPPRESERVE=y
#export USE_NOFSCHECK=y


NWCHEM DEPENDENCIES
===================
We build GA, blas and scalapack are also provided as nwc modules if the user does not have them.

cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchains/intel-openmpi-mkl.cmake -DNWCHEM_SUPERBUILD=OFF -DNWCHEM_INSTALL_DIR=.

#All vars can be provided in .bashrc or toolchain file

#when using env vars c,cpp,fort compilers need to be provided. optionally install-dir can be provided or set in env vars.
cmake .. -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DNWCHEM_SUPERBUILD=OFF -DNWCHEM_INSTALL_DIR=.



export TAMM_PROC_COUNT=1
export NWCHEM_PROC_COUNT=2
export NWCHEM_TOP=/opt/libraries/nwchem-devel
export NWCHEM_TARGET=LINUX64
#export ARMCI_NETWORK=OPENIB
#export ARMCI_NETWORK=MPI-PR                                                                  
export ARMCI_DEFAULT_SHMMAX_UBOUND=65536
export INTEL_64ALIGN=1
export USE_OPENMP=y
export NWCHEM_MODULES=all
export SCALAPACK_LIB="-mkl=parallel -lmkl_scalapack_ilp64 -lmkl_blacs_openmpi_ilp64 -lmkl_intel_thread -lpthread -lm"
export SCALAPACK="$SCALAPACK_LIB"
export LAPACK_LIB="-mkl -lmkl_blacs_openmpi_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm"
export BLAS_LIB="$LAPACK_LIB"
export BLASOPT="$LAPACK_LIB"
export USE_SCALAPACK=y
export SCALAPACK_SIZE=8
export BLAS_SIZE=8
export LAPACK_SIZE=8
export USE_MPI=y
export USE_MPIF=y
export USE_MPIF4=y
export PYTHONHOME=/home/panyala/anaconda2/
export PYTHONVERSION=2.7
export USE_PYTHONCONFIG=Y
export IPCCSD=y
export EACCSD=y
export MRCC_METHODS=y

unset MPI_INCLUDE
unset MPI_LIB
unset LIBMPI
unset GA_DEV


export NWCHEM_INSTALL_DIR=/opt/libraries/nwchem-cmake-install
export SCALAPACK_LIBRARIES="-mkl=parallel -lmkl_scalapack_ilp64 -lmkl_blacs_openmpi_ilp64 -lmkl_intel_thread -lpthread -lm"
export BLAS_LIBRARIES="-mkl -lmkl_blacs_openmpi_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm"
export LAPACK_LIBRARIES="-mkl -lmkl_blacs_openmpi_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm"

export MPI_INCLUDE_PATH=/usr/lib/openmpi/include/
export MPI_LIBRARIES="-lmpi_f77 -lmpi -ldl -lhwloc"

#Not needed, but since its passed to GA, we keep it
export MPI_LIBRARY_PATH=/usr/lib/openmpi/lib/



==========
SVN Repo error

sqlite3 .svn/wc.db "select * from work_queue"
sqlite3 .svn/wc.db "delete from work_queue"
svn cleanup
