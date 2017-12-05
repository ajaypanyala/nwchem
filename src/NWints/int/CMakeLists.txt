set(NWINTS_INT_SOURCE_FILES
showxyz.F
hf1_er.F
hferi_gen.F
hf1.F
hf1mke3.f
hf1mkr.f
hf2oi.f
hf3mkr.f
hf3OIs.f
hfabc.f
hfctr3.f
hfkei.f
hfnai.f
hf_print.f
hfreord_pq.f
hf3PEabc.F
hfmkr.F
hfset.F
case_impl.F
hf2mkr.F
hfnai_er.F
hf1_tran.F
igamma.F
hfmke.F
matchNxyz.F
int_order.F
hf2.F
hferi.F
getNxyz.F
defNxyz.F
hf3pot.F
hf1set3.F
)

add_library(nwints_int OBJECT ${NWINTS_INT_SOURCE_FILES})
target_include_directories(nwints_int PUBLIC ${NWC_INCLUDES} ${NWINTS_INCLUDES})
target_compile_options(nwints_int PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(nwints_int PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})