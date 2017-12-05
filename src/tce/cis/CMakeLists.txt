
set(CIS_SOURCE_FILES
cis.F
cis_xdiagon.F
cis_xguess.F
cis_x1.F
cis_xtidy.F
)

add_library(tce_cis STATIC ${CIS_SOURCE_FILES})
#add_dependencies(tce_cis ${NWC_DEPENDENCIES})
target_include_directories(tce_cis PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_cis PUBLIC DEBUG_PRINT)

# From GNUMAKEFILE - NOT HANDLED 
# For unknown reasons the IBM xlf compiler generated an ICE on Power7 systems
# while trying to optimize this code. This problem might go away in future.
# Date 01/11/2012.
# ifeq ($(NWCHEM_TARGET),$(findstring $(NWCHEM_TARGET),LAPI64 IBM64))
#   OBJ += cis_xguess.o
# else
#   OBJ_OPTIMIZE += cis_xguess.o
target_compile_options(tce_cis PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_cis PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})
