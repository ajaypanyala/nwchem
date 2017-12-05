
set(SELCI_SOURCE_FILES
./slect.F
./mkdele.F
./iodiff.F
./mkinta.F
./bserch.F
./eirerj.F
#./ointer.parthack.F
./prtri.F
./makehs.F
./fcidim.F
./conf.F
./getij.F
./anal.F
./mkindx.F
./rdhead.F
./restin.F
./loadg.F
./restk.F
./readh.F
./ptidy.F
./loadab.F
#./ointer.mod.F
./orbanal.F
./putinb.F
./popcnt32.F
./inicij.F
./makef.F
./makeh.F
./tester.F
./stool.F
./eijeji.F
./selci_select.F
./initci.F
./load.F
./mpdele.F
./pkcon.F
./makehd.F
./pconls.F
./gtijkl.F
./e0mp2.F
./iprtri.F
./loadw.F
./eijkl.F
./minit.F
./mknsneed.F
./rdhint.F
./rd1.F
./eijil.F
./selci_hpp.F
./ointer.F
./input.F
./selci.F
./makhdb.F
./sigma.F
./upkcon.F
./eijkj.F
./opencc.F
./mkvc.F
./loada.F
./ci.F
./davids.F
./eij.F
axb.f
axbt.f
gettim.f
idmin.f
isum.f
mkindb.f
mxma.f
mxmb.f
oeq.f
prtcon.f
rdconf.f
rdhcon.f
rdhwmt.f
selci_check.f
selci_icopy.f
selci_output.f
selci_renorm.f
selci_sread.f
selci_swrite.f
selci_yacobi.f
tqlrat.f
tred1.f
vadd.f
wrtcon.f
wtconf.f
wthcon.f

)

add_library(selciC OBJECT ./couple.c ./cputm.c) #spindim.c
target_include_directories(selciC PUBLIC ${NWC_INCLUDES})
target_compile_options(selciC PUBLIC ${NWCHEM_C_COMPILE_FLAGS})
target_compile_options(selciC PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})


add_library(selci STATIC ${SELCI_SOURCE_FILES} $<TARGET_OBJECTS:selciC>)
#add_dependencies(selci ${NWC_DEPENDENCIES})
target_include_directories(selci PUBLIC ${NWC_INCLUDES})
target_compile_options(selci PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(selci PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

