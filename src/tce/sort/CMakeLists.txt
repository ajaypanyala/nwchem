
set(SORT_SOURCE_FILES
#tce_sort.F
#tce_sort6.F
# tce_sortacc6.F
# tce_sort4_jeff_goldstone.F
# tce_sort2.F
# tce_sortacc2.F
# tce_sort4_jeff_torun.F
#tce_sortacc.F
# tce_sort4.F
# tce_sortacc4.F
# new_sort8.F

tce_sort0.F
tce_sort1.f  
tce_sortacc0.f
new_sort2.F
new_sort4.F
new_sort6.F
tce_sort8.F
tce_sortacc8.F
)


add_library(tcesortC OBJECT tce_sort_4kg.c)
target_include_directories(tcesortC PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_options(tcesortC PUBLIC ${NWCHEM_C_COMPILE_FLAGS})
target_compile_options(tcesortC PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})

add_library(tce_sort STATIC ${SORT_SOURCE_FILES} $<TARGET_OBJECTS:tcesortC>)
#add_dependencies(tce_sort ${NWC_DEPENDENCIES})
target_include_directories(tce_sort PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})

target_compile_options(tce_sort PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_sort PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
