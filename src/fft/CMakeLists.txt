
set(FFT_SOURCE_FILES
./pfftwrap/iessl.F
./pfftwrap/nwfft3d.F
./pfft1.0/fftpack/dradb2.F
./pfft1.0/fftpack/dpssf2.F
./pfft1.0/fftpack/drffti.F
./pfft1.0/fftpack/dcffti.F
./pfft1.0/fftpack/dpssf5.F
./pfft1.0/fftpack/dpssf4.F
./pfft1.0/fftpack/dpssb2.F
./pfft1.0/fftpack/drftf1.F
./pfft1.0/fftpack/dpssb5.F
./pfft1.0/fftpack/dradb4.F
./pfft1.0/fftpack/dradb3.F
./pfft1.0/fftpack/dradbg.F
./pfft1.0/fftpack/drfti1.F
./pfft1.0/fftpack/dpssf3.F
./pfft1.0/fftpack/dcfti1.F
./pfft1.0/fftpack/dradf4.F
./pfft1.0/fftpack/dcftb1.F
./pfft1.0/fftpack/drftb1.F
./pfft1.0/fftpack/dcftf1.F
./pfft1.0/fftpack/dpssb4.F
./pfft1.0/fftpack/dpssb3.F
./pfft1.0/fftpack/dradb5.F
./pfft1.0/fftpack/drfftb.F
./pfft1.0/fftpack/dradf3.F
./pfft1.0/fftpack/dradfg.F
./pfft1.0/fftpack/dcfftf.F
./pfft1.0/fftpack/dpssb.F
./pfft1.0/fftpack/dradf5.F
./pfft1.0/fftpack/dradf2.F
./pfft1.0/fftpack/dcfftb.F
./pfft1.0/fftpack/drfftf.F
./pfft1.0/fftpack/dpssf.F
./pfft1.0/pfft/wrkspc2.F
./pfft1.0/pfft/igdiff.F
./pfft1.0/pfft/icopy.F
./pfft1.0/pfft/blktrn.F
./pfft1.0/pfft/lstind.F
./pfft1.0/pfft/mapchk.F
./pfft1.0/pfft/pfft1.F
./pfft1.0/pfft/pfft1i.F
./pfft1.0/pfft/pairup.F
./pfft1.0/pfft/mrfftr.F
./pfft1.0/pfft/prfft.F
./pfft1.0/pfft/fillup.F
./pfft1.0/pfft/prffti.F
./pfft1.0/pfft/pfft3d.F
./pfft1.0/pfft/fillup2.F
./pfft1.0/pfft/mfftr1.F
./pfft1.0/pfft/mfftc1.F
./pfft1.0/pfft/prfft3.F
./pfft1.0/pfft/chk3d.F
./pfft1.0/pfft/drain2.F
./pfft1.0/pfft/dcopy.F
./pfft1.0/pfft/drain.F
./pfft1.0/pfft/mapfil.F
./pfft1.0/pfft/wrkspc.F
./pfft1.0/pfft/prlst2.F
./pfft1.0/pfft/gexit.F
./pfft1.0/pfft/lstcnt.F
./pfft1.0/pfft/mcfftr.F
./pfft1.0/pfft/xstop.F
./pfft1.0/pfft/pfft.F
)

#-D${IO_STYLE} -DIBM ${GENDEFS} -DTCGMSG -DFFTPACK -DEMSLFFT
add_library(pfft STATIC ${FFT_SOURCE_FILES})
target_include_directories(pfft PUBLIC ${NWC_INCLUDES})
target_compile_options(pfft PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(pfft PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

