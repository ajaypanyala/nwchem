
set(CCSD_T_SOURCE_FILES
ccsd2_t.F
offl_ccsd_t_doubles_l.F
ccsd_t.F
ccsdt2_q_right.F
ccsd_t_singles_l_6dts.F
ccsd_t_doubles_l_6dts.F
cr_ccsd_t_D.F
ccsd_t_doubles_l.F
ccsd_t_dot.F
ccsd_t_dts_kernels_omp.F
cr_qcisd_t.F
ccsd_t_restart.F
cr_ccsd_t_E.F
ccsd_t_singles_l.F
ccsd_t_6dts.F
offl_ccsd_t_singles_l_1.F
tce_hashnsort.F
lambda_ccsd_t_left.F
m2_dot_m1t1.F
offl_ccsd_t_singles_l.F
ccsd_t_kernels_omp.F
ccsdt2_q_left.F
lr_ccsd_tq.F
ccsd_t_doubles.F
ccsd2_q_left.F
ccsd_t_pstat.F
ccsdt2_q.F
ccsd_t_transpose.F
ccsd2_t_left.F
ccsd2_q.F
ccsd2_q_right.F
cr_ccsd_t_N.F
ccsd_t_singles.F
ccsd_t_neword.F
lambda_ccsd_t.F
)

if(TCE_NEW_OPENMP)
    list(APPEND CCSD_T_SOURCE_FILES ./ccsd_t_omp.F)
endif()



#TODO: Double check the GPU part
if(TCE_CUDA)  
    INCLUDE(FindCUDA)
    #set(CCSD_T_CUDA_SOURCE_FILES ccsd_t_singles_gpu.F ccsd_t_doubles_gpu.F 
    #            hybrid.c memory.cu sd_t_total.cu ccsd_t_gpu.F)
    if(CUDA_TOOLKIT_INCLUDE) 
	    list(APPEND TCE_INCLUDES ${CUDA_TOOLKIT_INCLUDE})
    else()
	    message(WARNING "CUDA_TOOLKIT_INCLUDE not set.")
    endif()

    add_library(tce_ccsd_t_gpu_f STATIC ccsd_t_singles_gpu.F ccsd_t_doubles_gpu.F ccsd_t_gpu.F)
    target_include_directories(tce_ccsd_t_gpu_f PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
    target_compile_definitions(tce_ccsd_t_gpu_f PUBLIC DEBUG_PRINT)
    target_compile_options(tce_ccsd_t_gpu_f PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
    target_compile_options(tce_ccsd_t_gpu_f PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

    add_library(tce_ccsd_t_gpu_c STATIC hybrid.c) 
    target_include_directories(tce_ccsd_t_gpu_c PUBLIC ${NWC_INCLUDES})
    target_compile_options(tce_ccsd_t_gpu_c PUBLIC ${NWCHEM_C_COMPILE_FLAGS})
    target_compile_options(tce_ccsd_t_gpu_c PUBLIC ${NWCHEM_C_OPTIMIZE_FLAGS})

if(CCSD_T_GPU_OLD)
    cuda_add_library(tce_ccsd_t_cuda STATIC memory.cu sd_t_total.cu direct_sd1.cu direct_sd2.cu)
else()
    cuda_add_library(tce_ccsd_t_cuda STATIC memory.cu sd_t_total_ttlg.cu)
endif()
    target_include_directories(tce_ccsd_t_cuda PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES} ${PROJECT_SOURCE_DIR}/src/tce/ttlg ${PROJECT_SOURCE_DIR}/src/tce/ttlg/includes)
    target_compile_definitions(tce_ccsd_t_cuda PUBLIC DEBUG_PRINT)
    target_compile_options(tce_ccsd_t_cuda PRIVATE -Xcompiler -fPIC -std=c++11 -Xptxas --warn-on-spills -arch=${GPU_ARCH})

    set(TCE_GPU_LIBS tce_ccsd_t_cuda tce_ccsd_t_gpu_f tce_ccsd_t_gpu_c)
    set(TCE_GPU_LIBS ${TCE_GPU_LIBS} PARENT_SCOPE)

    if(CCSD_T_GPU_OLD)
        set(TCE_CUDA_LIBS ${CUDA_cublas_LIBRARY} ${CUDA_CUDART_LIBRARY} -lstdc++)    
    else()
        set(TCE_CUDA_LIBS ${CUDA_cublas_LIBRARY} ${CUDA_CUDART_LIBRARY} ttlg -lstdc++)
    endif()
    set(TCE_CUDA_LIBS ${TCE_CUDA_LIBS} PARENT_SCOPE)
endif()

add_library(tce_ccsd_t_un OBJECT cr_ccsd_t.F lr_ccsd_t.F)
target_include_directories(tce_ccsd_t_un PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_ccsd_t_un PUBLIC DEBUG_PRINT)
target_compile_options(tce_ccsd_t_un PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS} ${NWCHEM_Fortran_DEBUG_FLAGS})

add_library(tce_ccsd_t STATIC ${CCSD_T_SOURCE_FILES} $<TARGET_OBJECTS:tce_ccsd_t_un>)
#add_dependencies(tce_ccsd_t ${NWC_DEPENDENCIES})
target_include_directories(tce_ccsd_t PUBLIC ${NWC_INCLUDES} ${TCE_INCLUDES})
target_compile_definitions(tce_ccsd_t PUBLIC DEBUG_PRINT)

if(OPENMP_OFFLOAD)
    target_compile_definitions(tce_ccsd_t PUBLIC OPENMP_OFFLOAD)
    target_compile_definitions(tce_ccsd_t_un PUBLIC OPENMP_OFFLOAD)
endif()

target_compile_options(tce_ccsd_t PUBLIC ${NWCHEM_Fortran_COMPILE_FLAGS})
target_compile_options(tce_ccsd_t PUBLIC ${NWCHEM_Fortran_OPTIMIZE_FLAGS})

