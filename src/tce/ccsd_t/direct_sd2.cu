#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <locale.h>
#include <algorithm>
#include "header.h"
using namespace std;

typedef long Integer;

#define NWCHEM_TCE_CCSD_T
#define CEIL(a, b) 		(((a) + (b) - 1) / (b))

// created by tc_gen_definition_new()
#define D2_1_SIZE_SLICE_1_G 16
#define D2_1_SIZE_SLICE_1_A 16
#define D2_1_SIZE_SLICE_1_D 4
#define D2_1_SIZE_SLICE_1_E 1
#define D2_1_SIZE_SLICE_1_F 8
#define D2_1_SIZE_SLICE_1_C 8
#define D2_1_SIZE_SLICE_1_B 1

#define D2_1_SIZE_INT_UNIT_1 D2_1_SIZE_SLICE_1_G

#define D2_1_SIZE_TB_1_X 	D2_1_SIZE_SLICE_1_A * D2_1_SIZE_SLICE_1_E
#define D2_1_SIZE_TB_1_Y 	D2_1_SIZE_SLICE_1_F * D2_1_SIZE_SLICE_1_B
#define D2_1_SIZE_REG_1_X 	D2_1_SIZE_SLICE_1_D
#define D2_1_SIZE_REG_1_Y 	D2_1_SIZE_SLICE_1_C

#ifndef NWCHEM_TCE_CCSD_T
#define CEIL(a, b) 		(((a) + (b) - 1) / (b))
#endif

// created by tc_gen_code_Kernel()
__global__ void d2_1_kernel__1_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_1_SIZE_SLICE_1_A;
	int idx_e = threadIdx.x / D2_1_SIZE_SLICE_1_A;
	int idx_f = threadIdx.y % D2_1_SIZE_SLICE_1_F;
	int idx_b = threadIdx.y / D2_1_SIZE_SLICE_1_F;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_1_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_1_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_1_SIZE_SLICE_1_C + (blk_idx_d * D2_1_SIZE_SLICE_1_D + (blk_idx_e * D2_1_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_1_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'y', 't2', ['g', 'f', 'c', 'b']], [16, 'STR_SD2_V2_H7', 'x', 'v2', ['g', 'a', 'd', 'e']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_1_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'f', 'c', 'b']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_f < rng_f
			sm_a[threadIdx.x][threadIdx.y + ll * 8] = dev_t2[(blk_idx_f * D2_1_SIZE_SLICE_1_F + idx_f + (blk_idx_c * D2_1_SIZE_SLICE_1_C + ll + (blk_idx_b * D2_1_SIZE_SLICE_1_B + 0) * size_c) * size_f) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'a', 'd', 'e']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_f + 0 < rng_a
			sm_b[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_v2[(blk_idx_a * D2_1_SIZE_SLICE_1_A + idx_f + 0 + (blk_idx_d * D2_1_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_1_SIZE_SLICE_1_E + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_f + 8 < rng_a
			// Exception: Full-Full
			sm_b[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_v2[(blk_idx_a * D2_1_SIZE_SLICE_1_A + idx_f + 8 + (blk_idx_d * D2_1_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_1_SIZE_SLICE_1_E + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_1_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 0];
			temp_bv[1] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 8];
			temp_bv[2] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 16];
			temp_bv[3] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 24];
			temp_bv[4] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 32];
			temp_bv[5] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 40];
			temp_bv[6] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 48];
			temp_bv[7] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_b[ll][idx_a + (idx_e) * D2_1_SIZE_SLICE_1_A + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_1_kernel__2_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_1_SIZE_SLICE_1_A;
	int idx_e = threadIdx.x / D2_1_SIZE_SLICE_1_A;
	int idx_f = threadIdx.y % D2_1_SIZE_SLICE_1_F;
	int idx_b = threadIdx.y / D2_1_SIZE_SLICE_1_F;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_1_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_1_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_1_SIZE_SLICE_1_C + (blk_idx_d * D2_1_SIZE_SLICE_1_D + (blk_idx_e * D2_1_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_1_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'y', 't2', ['g', 'f', 'c', 'b']], [16, 'STR_SD2_V2_H7', 'x', 'v2', ['g', 'a', 'd', 'e']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_1_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_1_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_1_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'f', 'c', 'b']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_f < rng_f
			sm_a[threadIdx.x][threadIdx.y + ll * 8] = dev_t2[(blk_idx_f * D2_1_SIZE_SLICE_1_F + idx_f + (blk_idx_c * D2_1_SIZE_SLICE_1_C + ll + (blk_idx_b * D2_1_SIZE_SLICE_1_B + 0) * size_c) * size_f) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_1_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'a', 'd', 'e']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_f + 0 < rng_a
			sm_b[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_v2[(blk_idx_a * D2_1_SIZE_SLICE_1_A + idx_f + 0 + (blk_idx_d * D2_1_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_1_SIZE_SLICE_1_E + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_f + 8 < rng_a
			if (threadIdx.x + l < size_internal) 
			sm_b[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_v2[(blk_idx_a * D2_1_SIZE_SLICE_1_A + idx_f + 8 + (blk_idx_d * D2_1_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_1_SIZE_SLICE_1_E + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_1_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 0];
			temp_bv[1] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 8];
			temp_bv[2] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 16];
			temp_bv[3] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 24];
			temp_bv[4] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 32];
			temp_bv[5] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 40];
			temp_bv[6] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 48];
			temp_bv[7] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_b[ll][idx_a + (idx_e) * D2_1_SIZE_SLICE_1_A + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_1_kernel__3_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_1_SIZE_SLICE_1_A;
	int idx_e = threadIdx.x / D2_1_SIZE_SLICE_1_A;
	int idx_f = threadIdx.y % D2_1_SIZE_SLICE_1_F;
	int idx_b = threadIdx.y / D2_1_SIZE_SLICE_1_F;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_1_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_1_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_1_SIZE_SLICE_1_C + (blk_idx_d * D2_1_SIZE_SLICE_1_D + (blk_idx_e * D2_1_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_1_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_1_SIZE_SLICE_1_A)) >= D2_1_SIZE_SLICE_1_A)
	{
		rng_a = D2_1_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_1_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_1_SIZE_SLICE_1_B)) >= D2_1_SIZE_SLICE_1_B)
	{
		rng_b = D2_1_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_1_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_1_SIZE_SLICE_1_C)) >= D2_1_SIZE_SLICE_1_C)
	{
		rng_c = D2_1_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_1_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_1_SIZE_SLICE_1_D)) >= D2_1_SIZE_SLICE_1_D)
	{
		rng_d = D2_1_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_1_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_1_SIZE_SLICE_1_E)) >= D2_1_SIZE_SLICE_1_E)
	{
		rng_e = D2_1_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_1_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_1_SIZE_SLICE_1_F)) >= D2_1_SIZE_SLICE_1_F)
	{
		rng_f = D2_1_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_1_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'y', 't2', ['g', 'f', 'c', 'b']], [16, 'STR_SD2_V2_H7', 'x', 'v2', ['g', 'a', 'd', 'e']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_1_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_f < rng_f && 0 < rng_b)
		for (int ll = 0; ll < rng_c; ll++)
		{
			// ['g', 'f', 'c', 'b']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_f < rng_f
			sm_a[threadIdx.x][threadIdx.y + ll * 8] = dev_t2[(blk_idx_f * D2_1_SIZE_SLICE_1_F + idx_f + (blk_idx_c * D2_1_SIZE_SLICE_1_C + ll + (blk_idx_b * D2_1_SIZE_SLICE_1_B + 0) * size_c) * size_f) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_f < rng_a && 0 < rng_e)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'a', 'd', 'e']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_f + 0 < rng_a
			sm_b[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_v2[(blk_idx_a * D2_1_SIZE_SLICE_1_A + idx_f + 0 + (blk_idx_d * D2_1_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_1_SIZE_SLICE_1_E + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_f + 8 < rng_a
			if (idx_f + 8 < rng_a) 
			sm_b[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_v2[(blk_idx_a * D2_1_SIZE_SLICE_1_A + idx_f + 8 + (blk_idx_d * D2_1_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_1_SIZE_SLICE_1_E + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_1_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 0];
			temp_bv[1] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 8];
			temp_bv[2] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 16];
			temp_bv[3] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 24];
			temp_bv[4] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 32];
			temp_bv[5] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 40];
			temp_bv[6] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 48];
			temp_bv[7] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_b[ll][idx_a + (idx_e) * D2_1_SIZE_SLICE_1_A + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_e < rng_e && idx_f < rng_f && idx_b < rng_b)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_c && j < rng_d)
			{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_1_kernel__4_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_1_SIZE_SLICE_1_A;
	int idx_e = threadIdx.x / D2_1_SIZE_SLICE_1_A;
	int idx_f = threadIdx.y % D2_1_SIZE_SLICE_1_F;
	int idx_b = threadIdx.y / D2_1_SIZE_SLICE_1_F;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_1_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_1_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_1_SIZE_SLICE_1_C + (blk_idx_d * D2_1_SIZE_SLICE_1_D + (blk_idx_e * D2_1_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_1_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_1_SIZE_SLICE_1_A)) >= D2_1_SIZE_SLICE_1_A)
	{
		rng_a = D2_1_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_1_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_1_SIZE_SLICE_1_B)) >= D2_1_SIZE_SLICE_1_B)
	{
		rng_b = D2_1_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_1_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_1_SIZE_SLICE_1_C)) >= D2_1_SIZE_SLICE_1_C)
	{
		rng_c = D2_1_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_1_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_1_SIZE_SLICE_1_D)) >= D2_1_SIZE_SLICE_1_D)
	{
		rng_d = D2_1_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_1_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_1_SIZE_SLICE_1_E)) >= D2_1_SIZE_SLICE_1_E)
	{
		rng_e = D2_1_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_1_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_1_SIZE_SLICE_1_F)) >= D2_1_SIZE_SLICE_1_F)
	{
		rng_f = D2_1_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_1_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'y', 't2', ['g', 'f', 'c', 'b']], [16, 'STR_SD2_V2_H7', 'x', 'v2', ['g', 'a', 'd', 'e']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_1_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_1_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_f < rng_f && 0 < rng_b && threadIdx.x < D2_1_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_c; ll++)
		{
			// ['g', 'f', 'c', 'b']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_f < rng_f
			sm_a[threadIdx.x][threadIdx.y + ll * 8] = dev_t2[(blk_idx_f * D2_1_SIZE_SLICE_1_F + idx_f + (blk_idx_c * D2_1_SIZE_SLICE_1_C + ll + (blk_idx_b * D2_1_SIZE_SLICE_1_B + 0) * size_c) * size_f) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_f < rng_a && 0 < rng_e && threadIdx.x < D2_1_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'a', 'd', 'e']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_f + 0 < rng_a
			sm_b[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_v2[(blk_idx_a * D2_1_SIZE_SLICE_1_A + idx_f + 0 + (blk_idx_d * D2_1_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_1_SIZE_SLICE_1_E + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_f + 8 < rng_a
			if (threadIdx.x + l < size_internal && idx_f + 8 < rng_a) 
			sm_b[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_v2[(blk_idx_a * D2_1_SIZE_SLICE_1_A + idx_f + 8 + (blk_idx_d * D2_1_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_1_SIZE_SLICE_1_E + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_1_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 0];
			temp_bv[1] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 8];
			temp_bv[2] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 16];
			temp_bv[3] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 24];
			temp_bv[4] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 32];
			temp_bv[5] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 40];
			temp_bv[6] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 48];
			temp_bv[7] = sm_a[ll][idx_f + (idx_b) * D2_1_SIZE_SLICE_1_F + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_b[ll][idx_a + (idx_e) * D2_1_SIZE_SLICE_1_A + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_e < rng_e && idx_f < rng_f && idx_b < rng_b)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_c && j < rng_d)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// written by tc_interface.tc_gen_code_interface_Header()
void sd_t_d2_1_cogent(int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, double* t3, double* host_t2, double* host_v2, int cond_kernel_1, int opt_register_transpose)
{
	int num_thread_blocks_kernel_1;

	double* dev_t3;
	double* dev_t2;
	double* dev_v2;


	num_thread_blocks_kernel_1 = CEIL(size_a, D2_1_SIZE_SLICE_1_A) * CEIL(size_b, D2_1_SIZE_SLICE_1_B) * CEIL(size_c, D2_1_SIZE_SLICE_1_C) * CEIL(size_d, D2_1_SIZE_SLICE_1_D) * CEIL(size_e, D2_1_SIZE_SLICE_1_E) * CEIL(size_f, D2_1_SIZE_SLICE_1_F);

#ifdef NWCHEM_TCE_CCSD_T
    dev_t3 = t3_d;
#else
    // cudaMalloc()
	cudaMalloc((void**) &dev_t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f);
#endif

    cudaMalloc((void**) &dev_t2, sizeof(double) * size_b * size_c * size_f * size_g);
	cudaMalloc((void**) &dev_v2, sizeof(double) * size_e * size_d * size_a * size_g);

#ifndef NWCHEM_TCE_CCSD_T
	// cudaMemcpy()
	cudaMemcpy(dev_t3, t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f, cudaMemcpyHostToDevice);
#endif
    cudaMemcpy(dev_t2, host_t2, sizeof(double) * size_b * size_c * size_f * size_g, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2, host_v2, sizeof(double) * size_e * size_d * size_a * size_g, cudaMemcpyHostToDevice);

#ifndef NWCHEM_TCE_CCSD_T
	// Related to Kernels
	// There are 1 Basic Kernels
	long long int tmp_operations = 2 * (long long int)(size_a * size_b * size_c * size_d * size_e * size_f) * size_g;
    
    printf ("========================================= fusedKernels =============================================\n");
	printf ("		Grid Size  : %6d (1D)\n", num_thread_blocks_kernel_1);
	printf ("		Block-size : %2d, %2d (2D)\n", D2_1_SIZE_TB_1_X, D2_1_SIZE_TB_1_Y);
	printf ("		Reg.-size  : %2d, %2d (2D)\n", D2_1_SIZE_REG_1_X, D2_1_SIZE_REG_1_Y);
	printf ("		A thread deals with (%d x %d) elements (basically)\n", D2_1_SIZE_TB_1_X * D2_1_SIZE_REG_1_X, D2_1_SIZE_TB_1_Y * D2_1_SIZE_REG_1_Y);
	printf ("		# of Operations: %lld\n", tmp_operations);
	printf ("====================================================================================================\n");
#endif
    
    dim3 gridsize_1(num_thread_blocks_kernel_1);
	dim3 blocksize_1(D2_1_SIZE_TB_1_X, D2_1_SIZE_TB_1_Y);

	int stride_output_a = 1;
	int stride_output_b = stride_output_a * size_a;
	int stride_output_c = stride_output_b * size_b;
	int stride_output_d = stride_output_c * size_c;
	int stride_output_e = stride_output_d * size_d;
	int stride_output_f = stride_output_e * size_e;

	int stride_reg_x_1 = stride_output_d;
	int stride_reg_y_1 = stride_output_c;

	int size_internal = size_g;

	int stride_int_t2 = 1;
	int stride_int_v2 = 1;

	// Decision Tree for Kernel Types
	// No Chance to Utilize the Register Transpose
	if (size_a % D2_1_SIZE_SLICE_1_A == 0 && size_b % D2_1_SIZE_SLICE_1_B == 0 && size_c % D2_1_SIZE_SLICE_1_C == 0 && size_d % D2_1_SIZE_SLICE_1_D == 0 && size_e % D2_1_SIZE_SLICE_1_E == 0 && size_f % D2_1_SIZE_SLICE_1_F == 0)
	{
		// [2] Extenral Index: Full
		if (size_g % D2_1_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Full && Internal: Full
			//printf ("External: Full, Internal: Full\n");
			d2_1_kernel__1_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_1_SIZE_SLICE_1_A), CEIL(size_b, D2_1_SIZE_SLICE_1_B), CEIL(size_c, D2_1_SIZE_SLICE_1_C), CEIL(size_d, D2_1_SIZE_SLICE_1_D), CEIL(size_e, D2_1_SIZE_SLICE_1_E), CEIL(size_f, D2_1_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Full && Internal: Partial
			//printf ("External: Full, Internal: Partial\n");
			d2_1_kernel__2_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_1_SIZE_SLICE_1_A), CEIL(size_b, D2_1_SIZE_SLICE_1_B), CEIL(size_c, D2_1_SIZE_SLICE_1_C), CEIL(size_d, D2_1_SIZE_SLICE_1_D), CEIL(size_e, D2_1_SIZE_SLICE_1_E), CEIL(size_f, D2_1_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}
	else
	{
		// [2] Extenral Index: Partial
		if (size_g % D2_1_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Partial && Internal: Full
			//printf ("External: Partial, Internal: Full\n");
			d2_1_kernel__3_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_1_SIZE_SLICE_1_A), CEIL(size_b, D2_1_SIZE_SLICE_1_B), CEIL(size_c, D2_1_SIZE_SLICE_1_C), CEIL(size_d, D2_1_SIZE_SLICE_1_D), CEIL(size_e, D2_1_SIZE_SLICE_1_E), CEIL(size_f, D2_1_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Partial && Internal: Partial
			//printf ("External: Partial, Internal: Partial\n");
			d2_1_kernel__4_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_1_SIZE_SLICE_1_A), CEIL(size_b, D2_1_SIZE_SLICE_1_B), CEIL(size_c, D2_1_SIZE_SLICE_1_C), CEIL(size_d, D2_1_SIZE_SLICE_1_D), CEIL(size_e, D2_1_SIZE_SLICE_1_E), CEIL(size_f, D2_1_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}

#ifndef NWCHEM_TCE_CCSD_T
	// Copy the Result from Device to Host
	cudaMemcpy(t3, dev_t3, sizeof(double) * (size_a * size_b * size_c * size_d * size_e * size_f), cudaMemcpyDeviceToHost);

	// cudaFree()
    cudaFree(dev_t3);
#endif
    
    cudaFree(dev_t2);	cudaFree(dev_v2);
}

// created by tc_gen_definition_new()
#define D2_2_SIZE_SLICE_1_G 16
#define D2_2_SIZE_SLICE_1_A 16
#define D2_2_SIZE_SLICE_1_F 4
#define D2_2_SIZE_SLICE_1_B 1
#define D2_2_SIZE_SLICE_1_C 8
#define D2_2_SIZE_SLICE_1_D 8
#define D2_2_SIZE_SLICE_1_E 1

#define D2_2_SIZE_INT_UNIT_1 D2_2_SIZE_SLICE_1_G

#define D2_2_SIZE_TB_1_X 	D2_2_SIZE_SLICE_1_A * D2_2_SIZE_SLICE_1_B
#define D2_2_SIZE_TB_1_Y 	D2_2_SIZE_SLICE_1_C * D2_2_SIZE_SLICE_1_E
#define D2_2_SIZE_REG_1_X 	D2_2_SIZE_SLICE_1_F
#define D2_2_SIZE_REG_1_Y 	D2_2_SIZE_SLICE_1_D

// created by tc_gen_code_Kernel()
__global__ void d2_2_kernel__1_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_2_SIZE_SLICE_1_A;
	int idx_b = threadIdx.x / D2_2_SIZE_SLICE_1_A;
	int idx_c = threadIdx.y % D2_2_SIZE_SLICE_1_C;
	int idx_e = threadIdx.y / D2_2_SIZE_SLICE_1_C;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_2_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_2_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_2_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_2_SIZE_SLICE_1_D + (blk_idx_e * D2_2_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_2_SIZE_SLICE_1_F) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'f', 'b', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'c', 'd', 'e']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_2_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'f', 'b', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_f * D2_2_SIZE_SLICE_1_F + ll + (blk_idx_b * D2_2_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_2_SIZE_SLICE_1_A + idx_c + 0) * size_b) * size_f) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			// Exception: Full-Full
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_f * D2_2_SIZE_SLICE_1_F + ll + (blk_idx_b * D2_2_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_2_SIZE_SLICE_1_A + idx_c + 8) * size_b) * size_f) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'c', 'd', 'e']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_c < rng_c
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_c * D2_2_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_2_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_2_SIZE_SLICE_1_E + 0) * size_d) * size_c) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_2_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 0];
			temp_bv[1] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 8];
			temp_bv[2] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 16];
			temp_bv[3] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 24];
			temp_bv[4] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 32];
			temp_bv[5] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 40];
			temp_bv[6] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 48];
			temp_bv[7] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_b + (idx_a) * D2_2_SIZE_SLICE_1_B + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_2_kernel__2_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_2_SIZE_SLICE_1_A;
	int idx_b = threadIdx.x / D2_2_SIZE_SLICE_1_A;
	int idx_c = threadIdx.y % D2_2_SIZE_SLICE_1_C;
	int idx_e = threadIdx.y / D2_2_SIZE_SLICE_1_C;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_2_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_2_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_2_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_2_SIZE_SLICE_1_D + (blk_idx_e * D2_2_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_2_SIZE_SLICE_1_F) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'f', 'b', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'c', 'd', 'e']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_2_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_2_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_2_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'f', 'b', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_f * D2_2_SIZE_SLICE_1_F + ll + (blk_idx_b * D2_2_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_2_SIZE_SLICE_1_A + idx_c + 0) * size_b) * size_f) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			if (threadIdx.x + l < size_internal) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_f * D2_2_SIZE_SLICE_1_F + ll + (blk_idx_b * D2_2_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_2_SIZE_SLICE_1_A + idx_c + 8) * size_b) * size_f) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_2_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'c', 'd', 'e']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_c < rng_c
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_c * D2_2_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_2_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_2_SIZE_SLICE_1_E + 0) * size_d) * size_c) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_2_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 0];
			temp_bv[1] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 8];
			temp_bv[2] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 16];
			temp_bv[3] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 24];
			temp_bv[4] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 32];
			temp_bv[5] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 40];
			temp_bv[6] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 48];
			temp_bv[7] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_b + (idx_a) * D2_2_SIZE_SLICE_1_B + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_2_kernel__3_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_2_SIZE_SLICE_1_A;
	int idx_b = threadIdx.x / D2_2_SIZE_SLICE_1_A;
	int idx_c = threadIdx.y % D2_2_SIZE_SLICE_1_C;
	int idx_e = threadIdx.y / D2_2_SIZE_SLICE_1_C;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_2_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_2_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_2_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_2_SIZE_SLICE_1_D + (blk_idx_e * D2_2_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_2_SIZE_SLICE_1_F) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_2_SIZE_SLICE_1_A)) >= D2_2_SIZE_SLICE_1_A)
	{
		rng_a = D2_2_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_2_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_2_SIZE_SLICE_1_B)) >= D2_2_SIZE_SLICE_1_B)
	{
		rng_b = D2_2_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_2_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_2_SIZE_SLICE_1_C)) >= D2_2_SIZE_SLICE_1_C)
	{
		rng_c = D2_2_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_2_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_2_SIZE_SLICE_1_D)) >= D2_2_SIZE_SLICE_1_D)
	{
		rng_d = D2_2_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_2_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_2_SIZE_SLICE_1_E)) >= D2_2_SIZE_SLICE_1_E)
	{
		rng_e = D2_2_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_2_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_2_SIZE_SLICE_1_F)) >= D2_2_SIZE_SLICE_1_F)
	{
		rng_f = D2_2_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_2_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'f', 'b', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'c', 'd', 'e']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_2_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (0 < rng_b && idx_c < rng_a)
		for (int ll = 0; ll < rng_f; ll++)
		{
			// ['g', 'f', 'b', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_f * D2_2_SIZE_SLICE_1_F + ll + (blk_idx_b * D2_2_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_2_SIZE_SLICE_1_A + idx_c + 0) * size_b) * size_f) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			if (0 < rng_b && idx_c + 8 < rng_a) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_f * D2_2_SIZE_SLICE_1_F + ll + (blk_idx_b * D2_2_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_2_SIZE_SLICE_1_A + idx_c + 8) * size_b) * size_f) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_c < rng_c && 0 < rng_e)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'c', 'd', 'e']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_c < rng_c
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_c * D2_2_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_2_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_2_SIZE_SLICE_1_E + 0) * size_d) * size_c) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_2_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 0];
			temp_bv[1] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 8];
			temp_bv[2] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 16];
			temp_bv[3] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 24];
			temp_bv[4] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 32];
			temp_bv[5] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 40];
			temp_bv[6] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 48];
			temp_bv[7] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_b + (idx_a) * D2_2_SIZE_SLICE_1_B + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_b < rng_b && idx_c < rng_c && idx_e < rng_e)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_d && j < rng_f)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_2_kernel__4_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_2_SIZE_SLICE_1_A;
	int idx_b = threadIdx.x / D2_2_SIZE_SLICE_1_A;
	int idx_c = threadIdx.y % D2_2_SIZE_SLICE_1_C;
	int idx_e = threadIdx.y / D2_2_SIZE_SLICE_1_C;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_2_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_2_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_2_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_2_SIZE_SLICE_1_D + (blk_idx_e * D2_2_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_2_SIZE_SLICE_1_F) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_2_SIZE_SLICE_1_A)) >= D2_2_SIZE_SLICE_1_A)
	{
		rng_a = D2_2_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_2_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_2_SIZE_SLICE_1_B)) >= D2_2_SIZE_SLICE_1_B)
	{
		rng_b = D2_2_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_2_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_2_SIZE_SLICE_1_C)) >= D2_2_SIZE_SLICE_1_C)
	{
		rng_c = D2_2_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_2_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_2_SIZE_SLICE_1_D)) >= D2_2_SIZE_SLICE_1_D)
	{
		rng_d = D2_2_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_2_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_2_SIZE_SLICE_1_E)) >= D2_2_SIZE_SLICE_1_E)
	{
		rng_e = D2_2_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_2_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_2_SIZE_SLICE_1_F)) >= D2_2_SIZE_SLICE_1_F)
	{
		rng_f = D2_2_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_2_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'f', 'b', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'c', 'd', 'e']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_2_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_2_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (0 < rng_b && idx_c < rng_a && threadIdx.x < D2_2_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_f; ll++)
		{
			// ['g', 'f', 'b', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_f * D2_2_SIZE_SLICE_1_F + ll + (blk_idx_b * D2_2_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_2_SIZE_SLICE_1_A + idx_c + 0) * size_b) * size_f) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			if (threadIdx.x + l < size_internal && 0 < rng_b && idx_c + 8 < rng_a) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_f * D2_2_SIZE_SLICE_1_F + ll + (blk_idx_b * D2_2_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_2_SIZE_SLICE_1_A + idx_c + 8) * size_b) * size_f) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_c < rng_c && 0 < rng_e && threadIdx.x < D2_2_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'c', 'd', 'e']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_c < rng_c
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_c * D2_2_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_2_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_2_SIZE_SLICE_1_E + 0) * size_d) * size_c) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_2_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 0];
			temp_bv[1] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 8];
			temp_bv[2] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 16];
			temp_bv[3] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 24];
			temp_bv[4] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 32];
			temp_bv[5] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 40];
			temp_bv[6] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 48];
			temp_bv[7] = sm_b[ll][idx_c + (idx_e) * D2_2_SIZE_SLICE_1_C + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_b + (idx_a) * D2_2_SIZE_SLICE_1_B + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_b < rng_b && idx_c < rng_c && idx_e < rng_e)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_d && j < rng_f)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// written by tc_interface.tc_gen_code_interface_Header()
void sd_t_d2_2_cogent(int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, double* t3, double* host_t2, double* host_v2, int cond_kernel_1, int opt_register_transpose)
{
	int num_thread_blocks_kernel_1;

	double* dev_t3;
	double* dev_t2;
	double* dev_v2;


	num_thread_blocks_kernel_1 = CEIL(size_a, D2_2_SIZE_SLICE_1_A) * CEIL(size_b, D2_2_SIZE_SLICE_1_B) * CEIL(size_c, D2_2_SIZE_SLICE_1_C) * CEIL(size_d, D2_2_SIZE_SLICE_1_D) * CEIL(size_e, D2_2_SIZE_SLICE_1_E) * CEIL(size_f, D2_2_SIZE_SLICE_1_F);

#ifdef NWCHEM_TCE_CCSD_T
    dev_t3 = t3_d;
#else
    // cudaMalloc()
	cudaMalloc((void**) &dev_t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f);
#endif

    cudaMalloc((void**) &dev_t2, sizeof(double) * size_a * size_b * size_f * size_g);
	cudaMalloc((void**) &dev_v2, sizeof(double) * size_e * size_d * size_c * size_g);

#ifndef NWCHEM_TCE_CCSD_T
	// cudaMemcpy()
	cudaMemcpy(dev_t3, t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f, cudaMemcpyHostToDevice);
#endif
    
    cudaMemcpy(dev_t2, host_t2, sizeof(double) * size_a * size_b * size_f * size_g, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2, host_v2, sizeof(double) * size_e * size_d * size_c * size_g, cudaMemcpyHostToDevice);

#ifndef NWCHEM_TCE_CCSD_T
	// Related to Kernels
	// There are 1 Basic Kernels
	long long int tmp_operations = 2 * (long long int)(size_a * size_b * size_c * size_d * size_e * size_f) * size_g;
    
    printf ("========================================= fusedKernels =============================================\n");
	printf ("		Grid Size  : %6d (1D)\n", num_thread_blocks_kernel_1);
	printf ("		Block-size : %2d, %2d (2D)\n", D2_2_SIZE_TB_1_X, D2_2_SIZE_TB_1_Y);
	printf ("		Reg.-size  : %2d, %2d (2D)\n", D2_2_SIZE_REG_1_X, D2_2_SIZE_REG_1_Y);
	printf ("		A thread deals with (%d x %d) elements (basically)\n", D2_2_SIZE_TB_1_X * D2_2_SIZE_REG_1_X, D2_2_SIZE_TB_1_Y * D2_2_SIZE_REG_1_Y);
	printf ("		# of Operations: %lld\n", tmp_operations);
	printf ("====================================================================================================\n");
#endif
    
    dim3 gridsize_1(num_thread_blocks_kernel_1);
	dim3 blocksize_1(D2_2_SIZE_TB_1_X, D2_2_SIZE_TB_1_Y);

	int stride_output_a = 1;
	int stride_output_b = stride_output_a * size_a;
	int stride_output_c = stride_output_b * size_b;
	int stride_output_d = stride_output_c * size_c;
	int stride_output_e = stride_output_d * size_d;
	int stride_output_f = stride_output_e * size_e;

	int stride_reg_x_1 = stride_output_f;
	int stride_reg_y_1 = stride_output_d;

	int size_internal = size_g;

	int stride_int_t2 = 1;
	int stride_int_v2 = 1;

	// Decision Tree for Kernel Types
	// No Chance to Utilize the Register Transpose
	if (size_a % D2_2_SIZE_SLICE_1_A == 0 && size_b % D2_2_SIZE_SLICE_1_B == 0 && size_c % D2_2_SIZE_SLICE_1_C == 0 && size_d % D2_2_SIZE_SLICE_1_D == 0 && size_e % D2_2_SIZE_SLICE_1_E == 0 && size_f % D2_2_SIZE_SLICE_1_F == 0)
	{
		// [2] Extenral Index: Full
		if (size_g % D2_2_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Full && Internal: Full
			//printf ("External: Full, Internal: Full\n");
			d2_2_kernel__1_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_2_SIZE_SLICE_1_A), CEIL(size_b, D2_2_SIZE_SLICE_1_B), CEIL(size_c, D2_2_SIZE_SLICE_1_C), CEIL(size_d, D2_2_SIZE_SLICE_1_D), CEIL(size_e, D2_2_SIZE_SLICE_1_E), CEIL(size_f, D2_2_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Full && Internal: Partial
			//printf ("External: Full, Internal: Partial\n");
			d2_2_kernel__2_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_2_SIZE_SLICE_1_A), CEIL(size_b, D2_2_SIZE_SLICE_1_B), CEIL(size_c, D2_2_SIZE_SLICE_1_C), CEIL(size_d, D2_2_SIZE_SLICE_1_D), CEIL(size_e, D2_2_SIZE_SLICE_1_E), CEIL(size_f, D2_2_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}
	else
	{
		// [2] Extenral Index: Partial
		if (size_g % D2_2_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Partial && Internal: Full
			//printf ("External: Partial, Internal: Full\n");
			d2_2_kernel__3_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_2_SIZE_SLICE_1_A), CEIL(size_b, D2_2_SIZE_SLICE_1_B), CEIL(size_c, D2_2_SIZE_SLICE_1_C), CEIL(size_d, D2_2_SIZE_SLICE_1_D), CEIL(size_e, D2_2_SIZE_SLICE_1_E), CEIL(size_f, D2_2_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Partial && Internal: Partial
			//printf ("External: Partial, Internal: Partial\n");
			d2_2_kernel__4_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_2_SIZE_SLICE_1_A), CEIL(size_b, D2_2_SIZE_SLICE_1_B), CEIL(size_c, D2_2_SIZE_SLICE_1_C), CEIL(size_d, D2_2_SIZE_SLICE_1_D), CEIL(size_e, D2_2_SIZE_SLICE_1_E), CEIL(size_f, D2_2_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
    }
    
#ifndef NWCHEM_TCE_CCSD_T
	// Copy the Result from Device to Host
	cudaMemcpy(t3, dev_t3, sizeof(double) * (size_a * size_b * size_c * size_d * size_e * size_f), cudaMemcpyDeviceToHost);

	// cudaFree()
    cudaFree(dev_t3);
#endif

    cudaFree(dev_t2);	cudaFree(dev_v2);
}

// created by tc_gen_definition_new()
#define D2_3_SIZE_SLICE_1_G 16
#define D2_3_SIZE_SLICE_1_A 16
#define D2_3_SIZE_SLICE_1_F 4
#define D2_3_SIZE_SLICE_1_C 1
#define D2_3_SIZE_SLICE_1_B 8
#define D2_3_SIZE_SLICE_1_D 8
#define D2_3_SIZE_SLICE_1_E 1

#define D2_3_SIZE_INT_UNIT_1 D2_3_SIZE_SLICE_1_G

#define D2_3_SIZE_TB_1_X 	D2_3_SIZE_SLICE_1_A * D2_3_SIZE_SLICE_1_C
#define D2_3_SIZE_TB_1_Y 	D2_3_SIZE_SLICE_1_B * D2_3_SIZE_SLICE_1_E
#define D2_3_SIZE_REG_1_X 	D2_3_SIZE_SLICE_1_F
#define D2_3_SIZE_REG_1_Y 	D2_3_SIZE_SLICE_1_D

// created by tc_gen_code_Kernel()
__global__ void d2_3_kernel__1_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_3_SIZE_SLICE_1_A;
	int idx_c = threadIdx.x / D2_3_SIZE_SLICE_1_A;
	int idx_b = threadIdx.y % D2_3_SIZE_SLICE_1_B;
	int idx_e = threadIdx.y / D2_3_SIZE_SLICE_1_B;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_3_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_3_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_3_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_3_SIZE_SLICE_1_D + (blk_idx_e * D2_3_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_3_SIZE_SLICE_1_F) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'f', 'c', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'b', 'd', 'e']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_3_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'f', 'c', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_f * D2_3_SIZE_SLICE_1_F + ll + (blk_idx_c * D2_3_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_3_SIZE_SLICE_1_A + idx_b + 0) * size_c) * size_f) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			// Exception: Full-Full
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_f * D2_3_SIZE_SLICE_1_F + ll + (blk_idx_c * D2_3_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_3_SIZE_SLICE_1_A + idx_b + 8) * size_c) * size_f) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'b', 'd', 'e']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_b < rng_b
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_b * D2_3_SIZE_SLICE_1_B + idx_b + (blk_idx_d * D2_3_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_3_SIZE_SLICE_1_E + 0) * size_d) * size_b) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_3_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 0];
			temp_bv[1] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 8];
			temp_bv[2] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 16];
			temp_bv[3] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 24];
			temp_bv[4] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 32];
			temp_bv[5] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 40];
			temp_bv[6] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 48];
			temp_bv[7] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_c + (idx_a) * D2_3_SIZE_SLICE_1_C + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_3_kernel__2_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_3_SIZE_SLICE_1_A;
	int idx_c = threadIdx.x / D2_3_SIZE_SLICE_1_A;
	int idx_b = threadIdx.y % D2_3_SIZE_SLICE_1_B;
	int idx_e = threadIdx.y / D2_3_SIZE_SLICE_1_B;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_3_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_3_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_3_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_3_SIZE_SLICE_1_D + (blk_idx_e * D2_3_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_3_SIZE_SLICE_1_F) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'f', 'c', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'b', 'd', 'e']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_3_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_3_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_3_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'f', 'c', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_f * D2_3_SIZE_SLICE_1_F + ll + (blk_idx_c * D2_3_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_3_SIZE_SLICE_1_A + idx_b + 0) * size_c) * size_f) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			if (threadIdx.x + l < size_internal) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_f * D2_3_SIZE_SLICE_1_F + ll + (blk_idx_c * D2_3_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_3_SIZE_SLICE_1_A + idx_b + 8) * size_c) * size_f) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_3_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'b', 'd', 'e']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_b < rng_b
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_b * D2_3_SIZE_SLICE_1_B + idx_b + (blk_idx_d * D2_3_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_3_SIZE_SLICE_1_E + 0) * size_d) * size_b) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_3_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 0];
			temp_bv[1] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 8];
			temp_bv[2] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 16];
			temp_bv[3] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 24];
			temp_bv[4] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 32];
			temp_bv[5] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 40];
			temp_bv[6] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 48];
			temp_bv[7] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_c + (idx_a) * D2_3_SIZE_SLICE_1_C + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_3_kernel__3_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_3_SIZE_SLICE_1_A;
	int idx_c = threadIdx.x / D2_3_SIZE_SLICE_1_A;
	int idx_b = threadIdx.y % D2_3_SIZE_SLICE_1_B;
	int idx_e = threadIdx.y / D2_3_SIZE_SLICE_1_B;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_3_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_3_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_3_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_3_SIZE_SLICE_1_D + (blk_idx_e * D2_3_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_3_SIZE_SLICE_1_F) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_3_SIZE_SLICE_1_A)) >= D2_3_SIZE_SLICE_1_A)
	{
		rng_a = D2_3_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_3_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_3_SIZE_SLICE_1_B)) >= D2_3_SIZE_SLICE_1_B)
	{
		rng_b = D2_3_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_3_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_3_SIZE_SLICE_1_C)) >= D2_3_SIZE_SLICE_1_C)
	{
		rng_c = D2_3_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_3_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_3_SIZE_SLICE_1_D)) >= D2_3_SIZE_SLICE_1_D)
	{
		rng_d = D2_3_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_3_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_3_SIZE_SLICE_1_E)) >= D2_3_SIZE_SLICE_1_E)
	{
		rng_e = D2_3_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_3_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_3_SIZE_SLICE_1_F)) >= D2_3_SIZE_SLICE_1_F)
	{
		rng_f = D2_3_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_3_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'f', 'c', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'b', 'd', 'e']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_3_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (0 < rng_c && idx_b < rng_a)
		for (int ll = 0; ll < rng_f; ll++)
		{
			// ['g', 'f', 'c', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_f * D2_3_SIZE_SLICE_1_F + ll + (blk_idx_c * D2_3_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_3_SIZE_SLICE_1_A + idx_b + 0) * size_c) * size_f) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			if (0 < rng_c && idx_b + 8 < rng_a) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_f * D2_3_SIZE_SLICE_1_F + ll + (blk_idx_c * D2_3_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_3_SIZE_SLICE_1_A + idx_b + 8) * size_c) * size_f) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_b < rng_b && 0 < rng_e)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'b', 'd', 'e']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_b < rng_b
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_b * D2_3_SIZE_SLICE_1_B + idx_b + (blk_idx_d * D2_3_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_3_SIZE_SLICE_1_E + 0) * size_d) * size_b) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_3_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 0];
			temp_bv[1] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 8];
			temp_bv[2] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 16];
			temp_bv[3] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 24];
			temp_bv[4] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 32];
			temp_bv[5] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 40];
			temp_bv[6] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 48];
			temp_bv[7] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_c + (idx_a) * D2_3_SIZE_SLICE_1_C + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_c < rng_c && idx_b < rng_b && idx_e < rng_e)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_d && j < rng_f)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_3_kernel__4_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_3_SIZE_SLICE_1_A;
	int idx_c = threadIdx.x / D2_3_SIZE_SLICE_1_A;
	int idx_b = threadIdx.y % D2_3_SIZE_SLICE_1_B;
	int idx_e = threadIdx.y / D2_3_SIZE_SLICE_1_B;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_3_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_3_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_3_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_3_SIZE_SLICE_1_D + (blk_idx_e * D2_3_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_3_SIZE_SLICE_1_F) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_3_SIZE_SLICE_1_A)) >= D2_3_SIZE_SLICE_1_A)
	{
		rng_a = D2_3_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_3_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_3_SIZE_SLICE_1_B)) >= D2_3_SIZE_SLICE_1_B)
	{
		rng_b = D2_3_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_3_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_3_SIZE_SLICE_1_C)) >= D2_3_SIZE_SLICE_1_C)
	{
		rng_c = D2_3_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_3_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_3_SIZE_SLICE_1_D)) >= D2_3_SIZE_SLICE_1_D)
	{
		rng_d = D2_3_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_3_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_3_SIZE_SLICE_1_E)) >= D2_3_SIZE_SLICE_1_E)
	{
		rng_e = D2_3_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_3_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_3_SIZE_SLICE_1_F)) >= D2_3_SIZE_SLICE_1_F)
	{
		rng_f = D2_3_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_3_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'f', 'c', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'b', 'd', 'e']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_3_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_3_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (0 < rng_c && idx_b < rng_a && threadIdx.x < D2_3_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_f; ll++)
		{
			// ['g', 'f', 'c', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_f * D2_3_SIZE_SLICE_1_F + ll + (blk_idx_c * D2_3_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_3_SIZE_SLICE_1_A + idx_b + 0) * size_c) * size_f) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			if (threadIdx.x + l < size_internal && 0 < rng_c && idx_b + 8 < rng_a) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_f * D2_3_SIZE_SLICE_1_F + ll + (blk_idx_c * D2_3_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_3_SIZE_SLICE_1_A + idx_b + 8) * size_c) * size_f) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_b < rng_b && 0 < rng_e && threadIdx.x < D2_3_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'b', 'd', 'e']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_b < rng_b
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_b * D2_3_SIZE_SLICE_1_B + idx_b + (blk_idx_d * D2_3_SIZE_SLICE_1_D + ll + (blk_idx_e * D2_3_SIZE_SLICE_1_E + 0) * size_d) * size_b) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_3_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 0];
			temp_bv[1] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 8];
			temp_bv[2] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 16];
			temp_bv[3] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 24];
			temp_bv[4] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 32];
			temp_bv[5] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 40];
			temp_bv[6] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 48];
			temp_bv[7] = sm_b[ll][idx_b + (idx_e) * D2_3_SIZE_SLICE_1_B + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_c + (idx_a) * D2_3_SIZE_SLICE_1_C + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_c < rng_c && idx_b < rng_b && idx_e < rng_e)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_d && j < rng_f)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// written by tc_interface.tc_gen_code_interface_Header()
void sd_t_d2_3_cogent(int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, double* t3, double* host_t2, double* host_v2, int cond_kernel_1, int opt_register_transpose)
{
	int num_thread_blocks_kernel_1;

	double* dev_t3;
	double* dev_t2;
	double* dev_v2;

    num_thread_blocks_kernel_1 = CEIL(size_a, D2_3_SIZE_SLICE_1_A) * CEIL(size_b, D2_3_SIZE_SLICE_1_B) * CEIL(size_c, D2_3_SIZE_SLICE_1_C) * CEIL(size_d, D2_3_SIZE_SLICE_1_D) * CEIL(size_e, D2_3_SIZE_SLICE_1_E) * CEIL(size_f, D2_3_SIZE_SLICE_1_F);

#ifdef NWCHEM_TCE_CCSD_T
    dev_t3 = t3_d;
#else
	// cudaMalloc()
	cudaMalloc((void**) &dev_t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f);
#endif 

    cudaMalloc((void**) &dev_t2, sizeof(double) * size_a * size_c * size_f * size_g);
	cudaMalloc((void**) &dev_v2, sizeof(double) * size_e * size_d * size_b * size_g);

#ifndef NWCHEM_TCE_CCSD_T
	// cudaMemcpy()
    cudaMemcpy(dev_t3, t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f, cudaMemcpyHostToDevice);
#endif    

	cudaMemcpy(dev_t2, host_t2, sizeof(double) * size_a * size_c * size_f * size_g, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2, host_v2, sizeof(double) * size_e * size_d * size_b * size_g, cudaMemcpyHostToDevice);

#ifndef NWCHEM_TCE_CCSD_T
	// Related to Kernels
	// There are 1 Basic Kernels
	long long int tmp_operations = 2 * (long long int)(size_a * size_b * size_c * size_d * size_e * size_f) * size_g;
    
    printf ("========================================= fusedKernels =============================================\n");
	printf ("		Grid Size  : %6d (1D)\n", num_thread_blocks_kernel_1);
	printf ("		Block-size : %2d, %2d (2D)\n", D2_3_SIZE_TB_1_X, D2_3_SIZE_TB_1_Y);
	printf ("		Reg.-size  : %2d, %2d (2D)\n", D2_3_SIZE_REG_1_X, D2_3_SIZE_REG_1_Y);
	printf ("		A thread deals with (%d x %d) elements (basically)\n", D2_3_SIZE_TB_1_X * D2_3_SIZE_REG_1_X, D2_3_SIZE_TB_1_Y * D2_3_SIZE_REG_1_Y);
	printf ("		# of Operations: %lld\n", tmp_operations);
	printf ("====================================================================================================\n");
#endif
    
    dim3 gridsize_1(num_thread_blocks_kernel_1);
	dim3 blocksize_1(D2_3_SIZE_TB_1_X, D2_3_SIZE_TB_1_Y);

	int stride_output_a = 1;
	int stride_output_b = stride_output_a * size_a;
	int stride_output_c = stride_output_b * size_b;
	int stride_output_d = stride_output_c * size_c;
	int stride_output_e = stride_output_d * size_d;
	int stride_output_f = stride_output_e * size_e;

	int stride_reg_x_1 = stride_output_f;
	int stride_reg_y_1 = stride_output_d;

	int size_internal = size_g;

	int stride_int_t2 = 1;
	int stride_int_v2 = 1;

	// Decision Tree for Kernel Types
	// No Chance to Utilize the Register Transpose
	if (size_a % D2_3_SIZE_SLICE_1_A == 0 && size_b % D2_3_SIZE_SLICE_1_B == 0 && size_c % D2_3_SIZE_SLICE_1_C == 0 && size_d % D2_3_SIZE_SLICE_1_D == 0 && size_e % D2_3_SIZE_SLICE_1_E == 0 && size_f % D2_3_SIZE_SLICE_1_F == 0)
	{
		// [2] Extenral Index: Full
		if (size_g % D2_3_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Full && Internal: Full
			//printf ("External: Full, Internal: Full\n");
			d2_3_kernel__1_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_3_SIZE_SLICE_1_A), CEIL(size_b, D2_3_SIZE_SLICE_1_B), CEIL(size_c, D2_3_SIZE_SLICE_1_C), CEIL(size_d, D2_3_SIZE_SLICE_1_D), CEIL(size_e, D2_3_SIZE_SLICE_1_E), CEIL(size_f, D2_3_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Full && Internal: Partial
			//printf ("External: Full, Internal: Partial\n");
			d2_3_kernel__2_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_3_SIZE_SLICE_1_A), CEIL(size_b, D2_3_SIZE_SLICE_1_B), CEIL(size_c, D2_3_SIZE_SLICE_1_C), CEIL(size_d, D2_3_SIZE_SLICE_1_D), CEIL(size_e, D2_3_SIZE_SLICE_1_E), CEIL(size_f, D2_3_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}
	else
	{
		// [2] Extenral Index: Partial
		if (size_g % D2_3_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Partial && Internal: Full
			//printf ("External: Partial, Internal: Full\n");
			d2_3_kernel__3_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_3_SIZE_SLICE_1_A), CEIL(size_b, D2_3_SIZE_SLICE_1_B), CEIL(size_c, D2_3_SIZE_SLICE_1_C), CEIL(size_d, D2_3_SIZE_SLICE_1_D), CEIL(size_e, D2_3_SIZE_SLICE_1_E), CEIL(size_f, D2_3_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Partial && Internal: Partial
			//printf ("External: Partial, Internal: Partial\n");
			d2_3_kernel__4_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_3_SIZE_SLICE_1_A), CEIL(size_b, D2_3_SIZE_SLICE_1_B), CEIL(size_c, D2_3_SIZE_SLICE_1_C), CEIL(size_d, D2_3_SIZE_SLICE_1_D), CEIL(size_e, D2_3_SIZE_SLICE_1_E), CEIL(size_f, D2_3_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}

#ifndef NWCHEM_TCE_CCSD_T
	// Copy the Result from Device to Host
	cudaMemcpy(t3, dev_t3, sizeof(double) * (size_a * size_b * size_c * size_d * size_e * size_f), cudaMemcpyDeviceToHost);

	// cudaFree()
    cudaFree(dev_t3);
#endif
    
    cudaFree(dev_t2);	cudaFree(dev_v2);
}

// created by tc_gen_definition_new()
#define D2_4_SIZE_SLICE_1_G 16
#define D2_4_SIZE_SLICE_1_A 16
#define D2_4_SIZE_SLICE_1_D 4
#define D2_4_SIZE_SLICE_1_F 1
#define D2_4_SIZE_SLICE_1_E 8
#define D2_4_SIZE_SLICE_1_C 8
#define D2_4_SIZE_SLICE_1_B 1

#define D2_4_SIZE_INT_UNIT_1 D2_4_SIZE_SLICE_1_G

#define D2_4_SIZE_TB_1_X 	D2_4_SIZE_SLICE_1_A * D2_4_SIZE_SLICE_1_F
#define D2_4_SIZE_TB_1_Y 	D2_4_SIZE_SLICE_1_E * D2_4_SIZE_SLICE_1_B
#define D2_4_SIZE_REG_1_X 	D2_4_SIZE_SLICE_1_D
#define D2_4_SIZE_REG_1_Y 	D2_4_SIZE_SLICE_1_C

// created by tc_gen_code_Kernel()
__global__ void d2_4_kernel__1_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_4_SIZE_SLICE_1_A;
	int idx_f = threadIdx.x / D2_4_SIZE_SLICE_1_A;
	int idx_e = threadIdx.y % D2_4_SIZE_SLICE_1_E;
	int idx_b = threadIdx.y / D2_4_SIZE_SLICE_1_E;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_4_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_4_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_4_SIZE_SLICE_1_C + (blk_idx_d * D2_4_SIZE_SLICE_1_D + (blk_idx_e * D2_4_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_4_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'y', 't2', ['g', 'e', 'c', 'b']], [16, 'STR_SD2_V2_H7', 'x', 'v2', ['g', 'a', 'd', 'f']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_4_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'e', 'c', 'b']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_e < rng_e
			sm_a[threadIdx.x][threadIdx.y + ll * 8] = dev_t2[(blk_idx_e * D2_4_SIZE_SLICE_1_E + idx_e + (blk_idx_c * D2_4_SIZE_SLICE_1_C + ll + (blk_idx_b * D2_4_SIZE_SLICE_1_B + 0) * size_c) * size_e) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'a', 'd', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_e + 0 < rng_a
			sm_b[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_v2[(blk_idx_a * D2_4_SIZE_SLICE_1_A + idx_e + 0 + (blk_idx_d * D2_4_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_4_SIZE_SLICE_1_F + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_e + 8 < rng_a
			// Exception: Full-Full
			sm_b[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_v2[(blk_idx_a * D2_4_SIZE_SLICE_1_A + idx_e + 8 + (blk_idx_d * D2_4_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_4_SIZE_SLICE_1_F + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_4_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 0];
			temp_bv[1] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 8];
			temp_bv[2] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 16];
			temp_bv[3] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 24];
			temp_bv[4] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 32];
			temp_bv[5] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 40];
			temp_bv[6] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 48];
			temp_bv[7] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_b[ll][idx_a + (idx_f) * D2_4_SIZE_SLICE_1_A + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_4_kernel__2_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_4_SIZE_SLICE_1_A;
	int idx_f = threadIdx.x / D2_4_SIZE_SLICE_1_A;
	int idx_e = threadIdx.y % D2_4_SIZE_SLICE_1_E;
	int idx_b = threadIdx.y / D2_4_SIZE_SLICE_1_E;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_4_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_4_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_4_SIZE_SLICE_1_C + (blk_idx_d * D2_4_SIZE_SLICE_1_D + (blk_idx_e * D2_4_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_4_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'y', 't2', ['g', 'e', 'c', 'b']], [16, 'STR_SD2_V2_H7', 'x', 'v2', ['g', 'a', 'd', 'f']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_4_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_4_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_4_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'e', 'c', 'b']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_e < rng_e
			sm_a[threadIdx.x][threadIdx.y + ll * 8] = dev_t2[(blk_idx_e * D2_4_SIZE_SLICE_1_E + idx_e + (blk_idx_c * D2_4_SIZE_SLICE_1_C + ll + (blk_idx_b * D2_4_SIZE_SLICE_1_B + 0) * size_c) * size_e) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_4_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'a', 'd', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_e + 0 < rng_a
			sm_b[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_v2[(blk_idx_a * D2_4_SIZE_SLICE_1_A + idx_e + 0 + (blk_idx_d * D2_4_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_4_SIZE_SLICE_1_F + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_e + 8 < rng_a
			if (threadIdx.x + l < size_internal) 
			sm_b[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_v2[(blk_idx_a * D2_4_SIZE_SLICE_1_A + idx_e + 8 + (blk_idx_d * D2_4_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_4_SIZE_SLICE_1_F + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_4_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 0];
			temp_bv[1] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 8];
			temp_bv[2] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 16];
			temp_bv[3] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 24];
			temp_bv[4] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 32];
			temp_bv[5] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 40];
			temp_bv[6] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 48];
			temp_bv[7] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_b[ll][idx_a + (idx_f) * D2_4_SIZE_SLICE_1_A + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_4_kernel__3_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_4_SIZE_SLICE_1_A;
	int idx_f = threadIdx.x / D2_4_SIZE_SLICE_1_A;
	int idx_e = threadIdx.y % D2_4_SIZE_SLICE_1_E;
	int idx_b = threadIdx.y / D2_4_SIZE_SLICE_1_E;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_4_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_4_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_4_SIZE_SLICE_1_C + (blk_idx_d * D2_4_SIZE_SLICE_1_D + (blk_idx_e * D2_4_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_4_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_4_SIZE_SLICE_1_A)) >= D2_4_SIZE_SLICE_1_A)
	{
		rng_a = D2_4_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_4_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_4_SIZE_SLICE_1_B)) >= D2_4_SIZE_SLICE_1_B)
	{
		rng_b = D2_4_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_4_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_4_SIZE_SLICE_1_C)) >= D2_4_SIZE_SLICE_1_C)
	{
		rng_c = D2_4_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_4_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_4_SIZE_SLICE_1_D)) >= D2_4_SIZE_SLICE_1_D)
	{
		rng_d = D2_4_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_4_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_4_SIZE_SLICE_1_E)) >= D2_4_SIZE_SLICE_1_E)
	{
		rng_e = D2_4_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_4_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_4_SIZE_SLICE_1_F)) >= D2_4_SIZE_SLICE_1_F)
	{
		rng_f = D2_4_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_4_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'y', 't2', ['g', 'e', 'c', 'b']], [16, 'STR_SD2_V2_H7', 'x', 'v2', ['g', 'a', 'd', 'f']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_4_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_e < rng_e && 0 < rng_b)
		for (int ll = 0; ll < rng_c; ll++)
		{
			// ['g', 'e', 'c', 'b']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_e < rng_e
			sm_a[threadIdx.x][threadIdx.y + ll * 8] = dev_t2[(blk_idx_e * D2_4_SIZE_SLICE_1_E + idx_e + (blk_idx_c * D2_4_SIZE_SLICE_1_C + ll + (blk_idx_b * D2_4_SIZE_SLICE_1_B + 0) * size_c) * size_e) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_e < rng_a && 0 < rng_f)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'a', 'd', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_e + 0 < rng_a
			sm_b[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_v2[(blk_idx_a * D2_4_SIZE_SLICE_1_A + idx_e + 0 + (blk_idx_d * D2_4_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_4_SIZE_SLICE_1_F + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_e + 8 < rng_a
			if (idx_e + 8 < rng_a) 
			sm_b[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_v2[(blk_idx_a * D2_4_SIZE_SLICE_1_A + idx_e + 8 + (blk_idx_d * D2_4_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_4_SIZE_SLICE_1_F + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_4_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 0];
			temp_bv[1] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 8];
			temp_bv[2] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 16];
			temp_bv[3] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 24];
			temp_bv[4] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 32];
			temp_bv[5] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 40];
			temp_bv[6] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 48];
			temp_bv[7] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_b[ll][idx_a + (idx_f) * D2_4_SIZE_SLICE_1_A + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_f < rng_f && idx_e < rng_e && idx_b < rng_b)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_c && j < rng_d)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_4_kernel__4_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_4_SIZE_SLICE_1_A;
	int idx_f = threadIdx.x / D2_4_SIZE_SLICE_1_A;
	int idx_e = threadIdx.y % D2_4_SIZE_SLICE_1_E;
	int idx_b = threadIdx.y / D2_4_SIZE_SLICE_1_E;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_4_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_4_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_4_SIZE_SLICE_1_C + (blk_idx_d * D2_4_SIZE_SLICE_1_D + (blk_idx_e * D2_4_SIZE_SLICE_1_E + idx_e + (blk_idx_f * D2_4_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_4_SIZE_SLICE_1_A)) >= D2_4_SIZE_SLICE_1_A)
	{
		rng_a = D2_4_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_4_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_4_SIZE_SLICE_1_B)) >= D2_4_SIZE_SLICE_1_B)
	{
		rng_b = D2_4_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_4_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_4_SIZE_SLICE_1_C)) >= D2_4_SIZE_SLICE_1_C)
	{
		rng_c = D2_4_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_4_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_4_SIZE_SLICE_1_D)) >= D2_4_SIZE_SLICE_1_D)
	{
		rng_d = D2_4_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_4_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_4_SIZE_SLICE_1_E)) >= D2_4_SIZE_SLICE_1_E)
	{
		rng_e = D2_4_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_4_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_4_SIZE_SLICE_1_F)) >= D2_4_SIZE_SLICE_1_F)
	{
		rng_f = D2_4_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_4_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'y', 't2', ['g', 'e', 'c', 'b']], [16, 'STR_SD2_V2_H7', 'x', 'v2', ['g', 'a', 'd', 'f']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_4_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_4_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_e < rng_e && 0 < rng_b && threadIdx.x < D2_4_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_c; ll++)
		{
			// ['g', 'e', 'c', 'b']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_e < rng_e
			sm_a[threadIdx.x][threadIdx.y + ll * 8] = dev_t2[(blk_idx_e * D2_4_SIZE_SLICE_1_E + idx_e + (blk_idx_c * D2_4_SIZE_SLICE_1_C + ll + (blk_idx_b * D2_4_SIZE_SLICE_1_B + 0) * size_c) * size_e) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_e < rng_a && 0 < rng_f && threadIdx.x < D2_4_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'a', 'd', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_e + 0 < rng_a
			sm_b[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_v2[(blk_idx_a * D2_4_SIZE_SLICE_1_A + idx_e + 0 + (blk_idx_d * D2_4_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_4_SIZE_SLICE_1_F + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_e + 8 < rng_a
			if (threadIdx.x + l < size_internal && idx_e + 8 < rng_a) 
			sm_b[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_v2[(blk_idx_a * D2_4_SIZE_SLICE_1_A + idx_e + 8 + (blk_idx_d * D2_4_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_4_SIZE_SLICE_1_F + 0) * size_d) * size_a) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_4_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 0];
			temp_bv[1] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 8];
			temp_bv[2] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 16];
			temp_bv[3] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 24];
			temp_bv[4] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 32];
			temp_bv[5] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 40];
			temp_bv[6] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 48];
			temp_bv[7] = sm_a[ll][idx_e + (idx_b) * D2_4_SIZE_SLICE_1_E + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_b[ll][idx_a + (idx_f) * D2_4_SIZE_SLICE_1_A + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_f < rng_f && idx_e < rng_e && idx_b < rng_b)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_c && j < rng_d)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// written by tc_interface.tc_gen_code_interface_Header()
void sd_t_d2_4_cogent(int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, double* t3, double* host_t2, double* host_v2, int cond_kernel_1, int opt_register_transpose)
{
	int num_thread_blocks_kernel_1;

	double* dev_t3;
	double* dev_t2;
	double* dev_v2;

	num_thread_blocks_kernel_1 = CEIL(size_a, D2_4_SIZE_SLICE_1_A) * CEIL(size_b, D2_4_SIZE_SLICE_1_B) * CEIL(size_c, D2_4_SIZE_SLICE_1_C) * CEIL(size_d, D2_4_SIZE_SLICE_1_D) * CEIL(size_e, D2_4_SIZE_SLICE_1_E) * CEIL(size_f, D2_4_SIZE_SLICE_1_F);
    
#ifdef NWCHEM_TCE_CCSD_T
    dev_t3 = t3_d;
#else
    // cudaMalloc()
    cudaMalloc((void**) &dev_t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f);
#endif
	cudaMalloc((void**) &dev_t2, sizeof(double) * size_b * size_c * size_e * size_g);
	cudaMalloc((void**) &dev_v2, sizeof(double) * size_f * size_d * size_a * size_g);

#ifndef NWCHEM_TCE_CCSD_T
	// cudaMemcpy()
    cudaMemcpy(dev_t3, t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f, cudaMemcpyHostToDevice);
#endif

	cudaMemcpy(dev_t2, host_t2, sizeof(double) * size_b * size_c * size_e * size_g, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2, host_v2, sizeof(double) * size_f * size_d * size_a * size_g, cudaMemcpyHostToDevice);

#ifndef NWCHEM_TCE_CCSD_T
	// Related to Kernels
	// There are 1 Basic Kernels
	long long int tmp_operations = 2 * (long long int)(size_a * size_b * size_c * size_d * size_e * size_f) * size_g;
    
    printf ("========================================= fusedKernels =============================================\n");
	printf ("		Grid Size  : %6d (1D)\n", num_thread_blocks_kernel_1);
	printf ("		Block-size : %2d, %2d (2D)\n", D2_4_SIZE_TB_1_X, D2_4_SIZE_TB_1_Y);
	printf ("		Reg.-size  : %2d, %2d (2D)\n", D2_4_SIZE_REG_1_X, D2_4_SIZE_REG_1_Y);
	printf ("		A thread deals with (%d x %d) elements (basically)\n", D2_4_SIZE_TB_1_X * D2_4_SIZE_REG_1_X, D2_4_SIZE_TB_1_Y * D2_4_SIZE_REG_1_Y);
	printf ("		# of Operations: %lld\n", tmp_operations);
	printf ("====================================================================================================\n");
#endif

    dim3 gridsize_1(num_thread_blocks_kernel_1);
	dim3 blocksize_1(D2_4_SIZE_TB_1_X, D2_4_SIZE_TB_1_Y);

	int stride_output_a = 1;
	int stride_output_b = stride_output_a * size_a;
	int stride_output_c = stride_output_b * size_b;
	int stride_output_d = stride_output_c * size_c;
	int stride_output_e = stride_output_d * size_d;
	int stride_output_f = stride_output_e * size_e;

	int stride_reg_x_1 = stride_output_d;
	int stride_reg_y_1 = stride_output_c;

	int size_internal = size_g;

	int stride_int_t2 = 1;
	int stride_int_v2 = 1;

	// Decision Tree for Kernel Types
	// No Chance to Utilize the Register Transpose
	if (size_a % D2_4_SIZE_SLICE_1_A == 0 && size_b % D2_4_SIZE_SLICE_1_B == 0 && size_c % D2_4_SIZE_SLICE_1_C == 0 && size_d % D2_4_SIZE_SLICE_1_D == 0 && size_e % D2_4_SIZE_SLICE_1_E == 0 && size_f % D2_4_SIZE_SLICE_1_F == 0)
	{
		// [2] Extenral Index: Full
		if (size_g % D2_4_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Full && Internal: Full
			//printf ("External: Full, Internal: Full\n");
			d2_4_kernel__1_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_4_SIZE_SLICE_1_A), CEIL(size_b, D2_4_SIZE_SLICE_1_B), CEIL(size_c, D2_4_SIZE_SLICE_1_C), CEIL(size_d, D2_4_SIZE_SLICE_1_D), CEIL(size_e, D2_4_SIZE_SLICE_1_E), CEIL(size_f, D2_4_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Full && Internal: Partial
			//printf ("External: Full, Internal: Partial\n");
			d2_4_kernel__2_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_4_SIZE_SLICE_1_A), CEIL(size_b, D2_4_SIZE_SLICE_1_B), CEIL(size_c, D2_4_SIZE_SLICE_1_C), CEIL(size_d, D2_4_SIZE_SLICE_1_D), CEIL(size_e, D2_4_SIZE_SLICE_1_E), CEIL(size_f, D2_4_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}
	else
	{
		// [2] Extenral Index: Partial
		if (size_g % D2_4_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Partial && Internal: Full
			//printf ("External: Partial, Internal: Full\n");
			d2_4_kernel__3_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_4_SIZE_SLICE_1_A), CEIL(size_b, D2_4_SIZE_SLICE_1_B), CEIL(size_c, D2_4_SIZE_SLICE_1_C), CEIL(size_d, D2_4_SIZE_SLICE_1_D), CEIL(size_e, D2_4_SIZE_SLICE_1_E), CEIL(size_f, D2_4_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Partial && Internal: Partial
			//printf ("External: Partial, Internal: Partial\n");
			d2_4_kernel__4_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_4_SIZE_SLICE_1_A), CEIL(size_b, D2_4_SIZE_SLICE_1_B), CEIL(size_c, D2_4_SIZE_SLICE_1_C), CEIL(size_d, D2_4_SIZE_SLICE_1_D), CEIL(size_e, D2_4_SIZE_SLICE_1_E), CEIL(size_f, D2_4_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}

#ifndef NWCHEM_TCE_CCSD_T
	// Copy the Result from Device to Host
	cudaMemcpy(t3, dev_t3, sizeof(double) * (size_a * size_b * size_c * size_d * size_e * size_f), cudaMemcpyDeviceToHost);

	// cudaFree()
    cudaFree(dev_t3);
#endif

    cudaFree(dev_t2);	cudaFree(dev_v2);
}

// created by tc_gen_definition_new()
#define D2_5_SIZE_SLICE_1_G 16
#define D2_5_SIZE_SLICE_1_A 16
#define D2_5_SIZE_SLICE_1_E 4
#define D2_5_SIZE_SLICE_1_B 1
#define D2_5_SIZE_SLICE_1_C 8
#define D2_5_SIZE_SLICE_1_D 8
#define D2_5_SIZE_SLICE_1_F 1

#define D2_5_SIZE_INT_UNIT_1 D2_5_SIZE_SLICE_1_G

#define D2_5_SIZE_TB_1_X 	D2_5_SIZE_SLICE_1_A * D2_5_SIZE_SLICE_1_B
#define D2_5_SIZE_TB_1_Y 	D2_5_SIZE_SLICE_1_C * D2_5_SIZE_SLICE_1_F
#define D2_5_SIZE_REG_1_X 	D2_5_SIZE_SLICE_1_E
#define D2_5_SIZE_REG_1_Y 	D2_5_SIZE_SLICE_1_D

// created by tc_gen_code_Kernel()
__global__ void d2_5_kernel__1_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_5_SIZE_SLICE_1_A;
	int idx_b = threadIdx.x / D2_5_SIZE_SLICE_1_A;
	int idx_c = threadIdx.y % D2_5_SIZE_SLICE_1_C;
	int idx_f = threadIdx.y / D2_5_SIZE_SLICE_1_C;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_5_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_5_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_5_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_5_SIZE_SLICE_1_D + (blk_idx_e * D2_5_SIZE_SLICE_1_E + (blk_idx_f * D2_5_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'e', 'b', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'c', 'd', 'f']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_5_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'e', 'b', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_e * D2_5_SIZE_SLICE_1_E + ll + (blk_idx_b * D2_5_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_5_SIZE_SLICE_1_A + idx_c + 0) * size_b) * size_e) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			// Exception: Full-Full
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_e * D2_5_SIZE_SLICE_1_E + ll + (blk_idx_b * D2_5_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_5_SIZE_SLICE_1_A + idx_c + 8) * size_b) * size_e) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'c', 'd', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_c < rng_c
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_c * D2_5_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_5_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_5_SIZE_SLICE_1_F + 0) * size_d) * size_c) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_5_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 0];
			temp_bv[1] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 8];
			temp_bv[2] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 16];
			temp_bv[3] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 24];
			temp_bv[4] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 32];
			temp_bv[5] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 40];
			temp_bv[6] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 48];
			temp_bv[7] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_b + (idx_a) * D2_5_SIZE_SLICE_1_B + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_5_kernel__2_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_5_SIZE_SLICE_1_A;
	int idx_b = threadIdx.x / D2_5_SIZE_SLICE_1_A;
	int idx_c = threadIdx.y % D2_5_SIZE_SLICE_1_C;
	int idx_f = threadIdx.y / D2_5_SIZE_SLICE_1_C;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_5_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_5_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_5_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_5_SIZE_SLICE_1_D + (blk_idx_e * D2_5_SIZE_SLICE_1_E + (blk_idx_f * D2_5_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'e', 'b', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'c', 'd', 'f']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_5_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_5_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_5_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'e', 'b', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_e * D2_5_SIZE_SLICE_1_E + ll + (blk_idx_b * D2_5_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_5_SIZE_SLICE_1_A + idx_c + 0) * size_b) * size_e) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			if (threadIdx.x + l < size_internal) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_e * D2_5_SIZE_SLICE_1_E + ll + (blk_idx_b * D2_5_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_5_SIZE_SLICE_1_A + idx_c + 8) * size_b) * size_e) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_5_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'c', 'd', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_c < rng_c
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_c * D2_5_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_5_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_5_SIZE_SLICE_1_F + 0) * size_d) * size_c) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_5_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 0];
			temp_bv[1] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 8];
			temp_bv[2] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 16];
			temp_bv[3] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 24];
			temp_bv[4] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 32];
			temp_bv[5] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 40];
			temp_bv[6] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 48];
			temp_bv[7] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_b + (idx_a) * D2_5_SIZE_SLICE_1_B + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_5_kernel__3_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_5_SIZE_SLICE_1_A;
	int idx_b = threadIdx.x / D2_5_SIZE_SLICE_1_A;
	int idx_c = threadIdx.y % D2_5_SIZE_SLICE_1_C;
	int idx_f = threadIdx.y / D2_5_SIZE_SLICE_1_C;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_5_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_5_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_5_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_5_SIZE_SLICE_1_D + (blk_idx_e * D2_5_SIZE_SLICE_1_E + (blk_idx_f * D2_5_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_5_SIZE_SLICE_1_A)) >= D2_5_SIZE_SLICE_1_A)
	{
		rng_a = D2_5_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_5_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_5_SIZE_SLICE_1_B)) >= D2_5_SIZE_SLICE_1_B)
	{
		rng_b = D2_5_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_5_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_5_SIZE_SLICE_1_C)) >= D2_5_SIZE_SLICE_1_C)
	{
		rng_c = D2_5_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_5_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_5_SIZE_SLICE_1_D)) >= D2_5_SIZE_SLICE_1_D)
	{
		rng_d = D2_5_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_5_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_5_SIZE_SLICE_1_E)) >= D2_5_SIZE_SLICE_1_E)
	{
		rng_e = D2_5_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_5_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_5_SIZE_SLICE_1_F)) >= D2_5_SIZE_SLICE_1_F)
	{
		rng_f = D2_5_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_5_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'e', 'b', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'c', 'd', 'f']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_5_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (0 < rng_b && idx_c < rng_a)
		for (int ll = 0; ll < rng_e; ll++)
		{
			// ['g', 'e', 'b', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_e * D2_5_SIZE_SLICE_1_E + ll + (blk_idx_b * D2_5_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_5_SIZE_SLICE_1_A + idx_c + 0) * size_b) * size_e) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			if (0 < rng_b && idx_c + 8 < rng_a) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_e * D2_5_SIZE_SLICE_1_E + ll + (blk_idx_b * D2_5_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_5_SIZE_SLICE_1_A + idx_c + 8) * size_b) * size_e) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_c < rng_c && 0 < rng_f)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'c', 'd', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_c < rng_c
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_c * D2_5_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_5_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_5_SIZE_SLICE_1_F + 0) * size_d) * size_c) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_5_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 0];
			temp_bv[1] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 8];
			temp_bv[2] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 16];
			temp_bv[3] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 24];
			temp_bv[4] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 32];
			temp_bv[5] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 40];
			temp_bv[6] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 48];
			temp_bv[7] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_b + (idx_a) * D2_5_SIZE_SLICE_1_B + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_b < rng_b && idx_c < rng_c && idx_f < rng_f)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_d && j < rng_e)
			{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_5_kernel__4_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_5_SIZE_SLICE_1_A;
	int idx_b = threadIdx.x / D2_5_SIZE_SLICE_1_A;
	int idx_c = threadIdx.y % D2_5_SIZE_SLICE_1_C;
	int idx_f = threadIdx.y / D2_5_SIZE_SLICE_1_C;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_5_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_5_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_5_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_5_SIZE_SLICE_1_D + (blk_idx_e * D2_5_SIZE_SLICE_1_E + (blk_idx_f * D2_5_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_5_SIZE_SLICE_1_A)) >= D2_5_SIZE_SLICE_1_A)
	{
		rng_a = D2_5_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_5_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_5_SIZE_SLICE_1_B)) >= D2_5_SIZE_SLICE_1_B)
	{
		rng_b = D2_5_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_5_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_5_SIZE_SLICE_1_C)) >= D2_5_SIZE_SLICE_1_C)
	{
		rng_c = D2_5_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_5_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_5_SIZE_SLICE_1_D)) >= D2_5_SIZE_SLICE_1_D)
	{
		rng_d = D2_5_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_5_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_5_SIZE_SLICE_1_E)) >= D2_5_SIZE_SLICE_1_E)
	{
		rng_e = D2_5_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_5_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_5_SIZE_SLICE_1_F)) >= D2_5_SIZE_SLICE_1_F)
	{
		rng_f = D2_5_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_5_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'e', 'b', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'c', 'd', 'f']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_5_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_5_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (0 < rng_b && idx_c < rng_a && threadIdx.x < D2_5_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_e; ll++)
		{
			// ['g', 'e', 'b', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_e * D2_5_SIZE_SLICE_1_E + ll + (blk_idx_b * D2_5_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_5_SIZE_SLICE_1_A + idx_c + 0) * size_b) * size_e) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			if (threadIdx.x + l < size_internal && 0 < rng_b && idx_c + 8 < rng_a) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_e * D2_5_SIZE_SLICE_1_E + ll + (blk_idx_b * D2_5_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_5_SIZE_SLICE_1_A + idx_c + 8) * size_b) * size_e) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_c < rng_c && 0 < rng_f && threadIdx.x < D2_5_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'c', 'd', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_c < rng_c
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_c * D2_5_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_5_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_5_SIZE_SLICE_1_F + 0) * size_d) * size_c) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_5_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 0];
			temp_bv[1] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 8];
			temp_bv[2] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 16];
			temp_bv[3] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 24];
			temp_bv[4] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 32];
			temp_bv[5] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 40];
			temp_bv[6] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 48];
			temp_bv[7] = sm_b[ll][idx_c + (idx_f) * D2_5_SIZE_SLICE_1_C + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_b + (idx_a) * D2_5_SIZE_SLICE_1_B + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_b < rng_b && idx_c < rng_c && idx_f < rng_f)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_d && j < rng_e)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// written by tc_interface.tc_gen_code_interface_Header()
void sd_t_d2_5_cogent(int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, double* t3, double* host_t2, double* host_v2, int cond_kernel_1, int opt_register_transpose)
{
	int num_thread_blocks_kernel_1;

	double* dev_t3;
	double* dev_t2;
	double* dev_v2;

	num_thread_blocks_kernel_1 = CEIL(size_a, D2_5_SIZE_SLICE_1_A) * CEIL(size_b, D2_5_SIZE_SLICE_1_B) * CEIL(size_c, D2_5_SIZE_SLICE_1_C) * CEIL(size_d, D2_5_SIZE_SLICE_1_D) * CEIL(size_e, D2_5_SIZE_SLICE_1_E) * CEIL(size_f, D2_5_SIZE_SLICE_1_F);

#ifdef NWCHEM_TCE_CCSD_T
    dev_t3 = t3_d;
#else
    // cudaMalloc()
    cudaMalloc((void**) &dev_t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f);
#endif    

	cudaMalloc((void**) &dev_t2, sizeof(double) * size_a * size_b * size_e * size_g);
	cudaMalloc((void**) &dev_v2, sizeof(double) * size_f * size_d * size_c * size_g);

#ifndef NWCHEM_TCE_CCSD_T
	// cudaMemcpy()
	cudaMemcpy(dev_t3, t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f, cudaMemcpyHostToDevice);
#endif

    cudaMemcpy(dev_t2, host_t2, sizeof(double) * size_a * size_b * size_e * size_g, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2, host_v2, sizeof(double) * size_f * size_d * size_c * size_g, cudaMemcpyHostToDevice);

#ifndef NWCHEM_TCE_CCSD_T
	// Related to Kernels
	// There are 1 Basic Kernels
	long long int tmp_operations = 2 * (long long int)(size_a * size_b * size_c * size_d * size_e * size_f) * size_g;
    
    printf ("========================================= fusedKernels =============================================\n");
	printf ("		Grid Size  : %6d (1D)\n", num_thread_blocks_kernel_1);
	printf ("		Block-size : %2d, %2d (2D)\n", D2_5_SIZE_TB_1_X, D2_5_SIZE_TB_1_Y);
	printf ("		Reg.-size  : %2d, %2d (2D)\n", D2_5_SIZE_REG_1_X, D2_5_SIZE_REG_1_Y);
	printf ("		A thread deals with (%d x %d) elements (basically)\n", D2_5_SIZE_TB_1_X * D2_5_SIZE_REG_1_X, D2_5_SIZE_TB_1_Y * D2_5_SIZE_REG_1_Y);
	printf ("		# of Operations: %lld\n", tmp_operations);
	printf ("====================================================================================================\n");
#endif
    
    dim3 gridsize_1(num_thread_blocks_kernel_1);
	dim3 blocksize_1(D2_5_SIZE_TB_1_X, D2_5_SIZE_TB_1_Y);

	int stride_output_a = 1;
	int stride_output_b = stride_output_a * size_a;
	int stride_output_c = stride_output_b * size_b;
	int stride_output_d = stride_output_c * size_c;
	int stride_output_e = stride_output_d * size_d;
	int stride_output_f = stride_output_e * size_e;

	int stride_reg_x_1 = stride_output_e;
	int stride_reg_y_1 = stride_output_d;

	int size_internal = size_g;

	int stride_int_t2 = 1;
	int stride_int_v2 = 1;

	// Decision Tree for Kernel Types
	// No Chance to Utilize the Register Transpose
	if (size_a % D2_5_SIZE_SLICE_1_A == 0 && size_b % D2_5_SIZE_SLICE_1_B == 0 && size_c % D2_5_SIZE_SLICE_1_C == 0 && size_d % D2_5_SIZE_SLICE_1_D == 0 && size_e % D2_5_SIZE_SLICE_1_E == 0 && size_f % D2_5_SIZE_SLICE_1_F == 0)
	{
		// [2] Extenral Index: Full
		if (size_g % D2_5_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Full && Internal: Full
			//printf ("External: Full, Internal: Full\n");
			d2_5_kernel__1_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_5_SIZE_SLICE_1_A), CEIL(size_b, D2_5_SIZE_SLICE_1_B), CEIL(size_c, D2_5_SIZE_SLICE_1_C), CEIL(size_d, D2_5_SIZE_SLICE_1_D), CEIL(size_e, D2_5_SIZE_SLICE_1_E), CEIL(size_f, D2_5_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Full && Internal: Partial
			//printf ("External: Full, Internal: Partial\n");
			d2_5_kernel__2_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_5_SIZE_SLICE_1_A), CEIL(size_b, D2_5_SIZE_SLICE_1_B), CEIL(size_c, D2_5_SIZE_SLICE_1_C), CEIL(size_d, D2_5_SIZE_SLICE_1_D), CEIL(size_e, D2_5_SIZE_SLICE_1_E), CEIL(size_f, D2_5_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}
	else
	{
		// [2] Extenral Index: Partial
		if (size_g % D2_5_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Partial && Internal: Full
			//printf ("External: Partial, Internal: Full\n");
			d2_5_kernel__3_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_5_SIZE_SLICE_1_A), CEIL(size_b, D2_5_SIZE_SLICE_1_B), CEIL(size_c, D2_5_SIZE_SLICE_1_C), CEIL(size_d, D2_5_SIZE_SLICE_1_D), CEIL(size_e, D2_5_SIZE_SLICE_1_E), CEIL(size_f, D2_5_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Partial && Internal: Partial
			//printf ("External: Partial, Internal: Partial\n");
			d2_5_kernel__4_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_5_SIZE_SLICE_1_A), CEIL(size_b, D2_5_SIZE_SLICE_1_B), CEIL(size_c, D2_5_SIZE_SLICE_1_C), CEIL(size_d, D2_5_SIZE_SLICE_1_D), CEIL(size_e, D2_5_SIZE_SLICE_1_E), CEIL(size_f, D2_5_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}

#ifndef NWCHEM_TCE_CCSD_T
	// Copy the Result from Device to Host
	cudaMemcpy(t3, dev_t3, sizeof(double) * (size_a * size_b * size_c * size_d * size_e * size_f), cudaMemcpyDeviceToHost);

	// cudaFree()
    cudaFree(dev_t3);
#endif

    cudaFree(dev_t2);	cudaFree(dev_v2);

	// Shoule be Fixed
	// HostFree

}

// created by tc_gen_definition_new()
#define D2_6_SIZE_SLICE_1_G 16
#define D2_6_SIZE_SLICE_1_A 16
#define D2_6_SIZE_SLICE_1_E 4
#define D2_6_SIZE_SLICE_1_C 1
#define D2_6_SIZE_SLICE_1_B 8
#define D2_6_SIZE_SLICE_1_D 8
#define D2_6_SIZE_SLICE_1_F 1

#define D2_6_SIZE_INT_UNIT_1 D2_6_SIZE_SLICE_1_G

#define D2_6_SIZE_TB_1_X 	D2_6_SIZE_SLICE_1_A * D2_6_SIZE_SLICE_1_C
#define D2_6_SIZE_TB_1_Y 	D2_6_SIZE_SLICE_1_B * D2_6_SIZE_SLICE_1_F
#define D2_6_SIZE_REG_1_X 	D2_6_SIZE_SLICE_1_E
#define D2_6_SIZE_REG_1_Y 	D2_6_SIZE_SLICE_1_D

// created by tc_gen_code_Kernel()
__global__ void d2_6_kernel__1_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_6_SIZE_SLICE_1_A;
	int idx_c = threadIdx.x / D2_6_SIZE_SLICE_1_A;
	int idx_b = threadIdx.y % D2_6_SIZE_SLICE_1_B;
	int idx_f = threadIdx.y / D2_6_SIZE_SLICE_1_B;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_6_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_6_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_6_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_6_SIZE_SLICE_1_D + (blk_idx_e * D2_6_SIZE_SLICE_1_E + (blk_idx_f * D2_6_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'e', 'c', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'b', 'd', 'f']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_6_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'e', 'c', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_e * D2_6_SIZE_SLICE_1_E + ll + (blk_idx_c * D2_6_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_6_SIZE_SLICE_1_A + idx_b + 0) * size_c) * size_e) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			// Exception: Full-Full
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_e * D2_6_SIZE_SLICE_1_E + ll + (blk_idx_c * D2_6_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_6_SIZE_SLICE_1_A + idx_b + 8) * size_c) * size_e) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'b', 'd', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_b < rng_b
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_b * D2_6_SIZE_SLICE_1_B + idx_b + (blk_idx_d * D2_6_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_6_SIZE_SLICE_1_F + 0) * size_d) * size_b) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_6_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 0];
			temp_bv[1] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 8];
			temp_bv[2] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 16];
			temp_bv[3] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 24];
			temp_bv[4] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 32];
			temp_bv[5] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 40];
			temp_bv[6] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 48];
			temp_bv[7] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_c + (idx_a) * D2_6_SIZE_SLICE_1_C + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_6_kernel__2_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_6_SIZE_SLICE_1_A;
	int idx_c = threadIdx.x / D2_6_SIZE_SLICE_1_A;
	int idx_b = threadIdx.y % D2_6_SIZE_SLICE_1_B;
	int idx_f = threadIdx.y / D2_6_SIZE_SLICE_1_B;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_6_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_6_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_6_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_6_SIZE_SLICE_1_D + (blk_idx_e * D2_6_SIZE_SLICE_1_E + (blk_idx_f * D2_6_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'e', 'c', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'b', 'd', 'f']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_6_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_6_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_6_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'e', 'c', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_e * D2_6_SIZE_SLICE_1_E + ll + (blk_idx_c * D2_6_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_6_SIZE_SLICE_1_A + idx_b + 0) * size_c) * size_e) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			if (threadIdx.x + l < size_internal) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_e * D2_6_SIZE_SLICE_1_E + ll + (blk_idx_c * D2_6_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_6_SIZE_SLICE_1_A + idx_b + 8) * size_c) * size_e) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_6_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'b', 'd', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_b < rng_b
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_b * D2_6_SIZE_SLICE_1_B + idx_b + (blk_idx_d * D2_6_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_6_SIZE_SLICE_1_F + 0) * size_d) * size_b) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_6_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 0];
			temp_bv[1] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 8];
			temp_bv[2] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 16];
			temp_bv[3] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 24];
			temp_bv[4] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 32];
			temp_bv[5] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 40];
			temp_bv[6] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 48];
			temp_bv[7] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_c + (idx_a) * D2_6_SIZE_SLICE_1_C + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_6_kernel__3_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_6_SIZE_SLICE_1_A;
	int idx_c = threadIdx.x / D2_6_SIZE_SLICE_1_A;
	int idx_b = threadIdx.y % D2_6_SIZE_SLICE_1_B;
	int idx_f = threadIdx.y / D2_6_SIZE_SLICE_1_B;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_6_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_6_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_6_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_6_SIZE_SLICE_1_D + (blk_idx_e * D2_6_SIZE_SLICE_1_E + (blk_idx_f * D2_6_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_6_SIZE_SLICE_1_A)) >= D2_6_SIZE_SLICE_1_A)
	{
		rng_a = D2_6_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_6_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_6_SIZE_SLICE_1_B)) >= D2_6_SIZE_SLICE_1_B)
	{
		rng_b = D2_6_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_6_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_6_SIZE_SLICE_1_C)) >= D2_6_SIZE_SLICE_1_C)
	{
		rng_c = D2_6_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_6_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_6_SIZE_SLICE_1_D)) >= D2_6_SIZE_SLICE_1_D)
	{
		rng_d = D2_6_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_6_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_6_SIZE_SLICE_1_E)) >= D2_6_SIZE_SLICE_1_E)
	{
		rng_e = D2_6_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_6_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_6_SIZE_SLICE_1_F)) >= D2_6_SIZE_SLICE_1_F)
	{
		rng_f = D2_6_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_6_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'e', 'c', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'b', 'd', 'f']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_6_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (0 < rng_c && idx_b < rng_a)
		for (int ll = 0; ll < rng_e; ll++)
		{
			// ['g', 'e', 'c', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_e * D2_6_SIZE_SLICE_1_E + ll + (blk_idx_c * D2_6_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_6_SIZE_SLICE_1_A + idx_b + 0) * size_c) * size_e) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			if (0 < rng_c && idx_b + 8 < rng_a) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_e * D2_6_SIZE_SLICE_1_E + ll + (blk_idx_c * D2_6_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_6_SIZE_SLICE_1_A + idx_b + 8) * size_c) * size_e) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_b < rng_b && 0 < rng_f)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'b', 'd', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_b < rng_b
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_b * D2_6_SIZE_SLICE_1_B + idx_b + (blk_idx_d * D2_6_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_6_SIZE_SLICE_1_F + 0) * size_d) * size_b) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_6_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 0];
			temp_bv[1] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 8];
			temp_bv[2] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 16];
			temp_bv[3] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 24];
			temp_bv[4] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 32];
			temp_bv[5] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 40];
			temp_bv[6] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 48];
			temp_bv[7] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_c + (idx_a) * D2_6_SIZE_SLICE_1_C + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_c < rng_c && idx_b < rng_b && idx_f < rng_f)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_d && j < rng_e)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_6_kernel__4_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_6_SIZE_SLICE_1_A;
	int idx_c = threadIdx.x / D2_6_SIZE_SLICE_1_A;
	int idx_b = threadIdx.y % D2_6_SIZE_SLICE_1_B;
	int idx_f = threadIdx.y / D2_6_SIZE_SLICE_1_B;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_6_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_6_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_6_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_6_SIZE_SLICE_1_D + (blk_idx_e * D2_6_SIZE_SLICE_1_E + (blk_idx_f * D2_6_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_6_SIZE_SLICE_1_A)) >= D2_6_SIZE_SLICE_1_A)
	{
		rng_a = D2_6_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_6_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_6_SIZE_SLICE_1_B)) >= D2_6_SIZE_SLICE_1_B)
	{
		rng_b = D2_6_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_6_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_6_SIZE_SLICE_1_C)) >= D2_6_SIZE_SLICE_1_C)
	{
		rng_c = D2_6_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_6_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_6_SIZE_SLICE_1_D)) >= D2_6_SIZE_SLICE_1_D)
	{
		rng_d = D2_6_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_6_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_6_SIZE_SLICE_1_E)) >= D2_6_SIZE_SLICE_1_E)
	{
		rng_e = D2_6_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_6_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_6_SIZE_SLICE_1_F)) >= D2_6_SIZE_SLICE_1_F)
	{
		rng_f = D2_6_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_6_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'e', 'c', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'b', 'd', 'f']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_6_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_6_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (0 < rng_c && idx_b < rng_a && threadIdx.x < D2_6_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_e; ll++)
		{
			// ['g', 'e', 'c', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_e * D2_6_SIZE_SLICE_1_E + ll + (blk_idx_c * D2_6_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_6_SIZE_SLICE_1_A + idx_b + 0) * size_c) * size_e) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			if (threadIdx.x + l < size_internal && 0 < rng_c && idx_b + 8 < rng_a) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_e * D2_6_SIZE_SLICE_1_E + ll + (blk_idx_c * D2_6_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_6_SIZE_SLICE_1_A + idx_b + 8) * size_c) * size_e) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_b < rng_b && 0 < rng_f && threadIdx.x < D2_6_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'b', 'd', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_b < rng_b
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_b * D2_6_SIZE_SLICE_1_B + idx_b + (blk_idx_d * D2_6_SIZE_SLICE_1_D + ll + (blk_idx_f * D2_6_SIZE_SLICE_1_F + 0) * size_d) * size_b) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_6_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 0];
			temp_bv[1] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 8];
			temp_bv[2] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 16];
			temp_bv[3] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 24];
			temp_bv[4] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 32];
			temp_bv[5] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 40];
			temp_bv[6] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 48];
			temp_bv[7] = sm_b[ll][idx_b + (idx_f) * D2_6_SIZE_SLICE_1_B + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_c + (idx_a) * D2_6_SIZE_SLICE_1_C + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_c < rng_c && idx_b < rng_b && idx_f < rng_f)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_d && j < rng_e)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// written by tc_interface.tc_gen_code_interface_Header()
void sd_t_d2_6_cogent(int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, double* t3, double* host_t2, double* host_v2, int cond_kernel_1, int opt_register_transpose)
{
	int num_thread_blocks_kernel_1;

	double* dev_t3;
	double* dev_t2;
	double* dev_v2;

	num_thread_blocks_kernel_1 = CEIL(size_a, D2_6_SIZE_SLICE_1_A) * CEIL(size_b, D2_6_SIZE_SLICE_1_B) * CEIL(size_c, D2_6_SIZE_SLICE_1_C) * CEIL(size_d, D2_6_SIZE_SLICE_1_D) * CEIL(size_e, D2_6_SIZE_SLICE_1_E) * CEIL(size_f, D2_6_SIZE_SLICE_1_F);

#ifdef NWCHEM_TCE_CCSD_T
    dev_t3 = t3_d;
#else
    // cudaMalloc()
    cudaMalloc((void**) &dev_t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f);
#endif
	cudaMalloc((void**) &dev_t2, sizeof(double) * size_a * size_c * size_e * size_g);
	cudaMalloc((void**) &dev_v2, sizeof(double) * size_f * size_d * size_b * size_g);

#ifndef NWCHEM_TCE_CCSD_T
	// cudaMemcpy()
	cudaMemcpy(dev_t3, t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f, cudaMemcpyHostToDevice);
#endif

    cudaMemcpy(dev_t2, host_t2, sizeof(double) * size_a * size_c * size_e * size_g, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2, host_v2, sizeof(double) * size_f * size_d * size_b * size_g, cudaMemcpyHostToDevice);

#ifndef NWCHEM_TCE_CCSD_T
	// Related to Kernels
	// There are 1 Basic Kernels
	long long int tmp_operations = 2 * (long long int)(size_a * size_b * size_c * size_d * size_e * size_f) * size_g;

    printf ("========================================= fusedKernels =============================================\n");
	printf ("		Grid Size  : %6d (1D)\n", num_thread_blocks_kernel_1);
	printf ("		Block-size : %2d, %2d (2D)\n", D2_6_SIZE_TB_1_X, D2_6_SIZE_TB_1_Y);
	printf ("		Reg.-size  : %2d, %2d (2D)\n", D2_6_SIZE_REG_1_X, D2_6_SIZE_REG_1_Y);
	printf ("		A thread deals with (%d x %d) elements (basically)\n", D2_6_SIZE_TB_1_X * D2_6_SIZE_REG_1_X, D2_6_SIZE_TB_1_Y * D2_6_SIZE_REG_1_Y);
	printf ("		# of Operations: %lld\n", tmp_operations);
	printf ("====================================================================================================\n");
#endif

    dim3 gridsize_1(num_thread_blocks_kernel_1);
	dim3 blocksize_1(D2_6_SIZE_TB_1_X, D2_6_SIZE_TB_1_Y);

	int stride_output_a = 1;
	int stride_output_b = stride_output_a * size_a;
	int stride_output_c = stride_output_b * size_b;
	int stride_output_d = stride_output_c * size_c;
	int stride_output_e = stride_output_d * size_d;
	int stride_output_f = stride_output_e * size_e;

	int stride_reg_x_1 = stride_output_e;
	int stride_reg_y_1 = stride_output_d;

	int size_internal = size_g;

	int stride_int_t2 = 1;
	int stride_int_v2 = 1;

	// Decision Tree for Kernel Types
	// No Chance to Utilize the Register Transpose
	if (size_a % D2_6_SIZE_SLICE_1_A == 0 && size_b % D2_6_SIZE_SLICE_1_B == 0 && size_c % D2_6_SIZE_SLICE_1_C == 0 && size_d % D2_6_SIZE_SLICE_1_D == 0 && size_e % D2_6_SIZE_SLICE_1_E == 0 && size_f % D2_6_SIZE_SLICE_1_F == 0)
	{
		// [2] Extenral Index: Full
		if (size_g % D2_6_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Full && Internal: Full
			//printf ("External: Full, Internal: Full\n");
			d2_6_kernel__1_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_6_SIZE_SLICE_1_A), CEIL(size_b, D2_6_SIZE_SLICE_1_B), CEIL(size_c, D2_6_SIZE_SLICE_1_C), CEIL(size_d, D2_6_SIZE_SLICE_1_D), CEIL(size_e, D2_6_SIZE_SLICE_1_E), CEIL(size_f, D2_6_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Full && Internal: Partial
			//printf ("External: Full, Internal: Partial\n");
			d2_6_kernel__2_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_6_SIZE_SLICE_1_A), CEIL(size_b, D2_6_SIZE_SLICE_1_B), CEIL(size_c, D2_6_SIZE_SLICE_1_C), CEIL(size_d, D2_6_SIZE_SLICE_1_D), CEIL(size_e, D2_6_SIZE_SLICE_1_E), CEIL(size_f, D2_6_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}
	else
	{
		// [2] Extenral Index: Partial
		if (size_g % D2_6_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Partial && Internal: Full
			//printf ("External: Partial, Internal: Full\n");
			d2_6_kernel__3_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_6_SIZE_SLICE_1_A), CEIL(size_b, D2_6_SIZE_SLICE_1_B), CEIL(size_c, D2_6_SIZE_SLICE_1_C), CEIL(size_d, D2_6_SIZE_SLICE_1_D), CEIL(size_e, D2_6_SIZE_SLICE_1_E), CEIL(size_f, D2_6_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Partial && Internal: Partial
			//printf ("External: Partial, Internal: Partial\n");
			d2_6_kernel__4_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_6_SIZE_SLICE_1_A), CEIL(size_b, D2_6_SIZE_SLICE_1_B), CEIL(size_c, D2_6_SIZE_SLICE_1_C), CEIL(size_d, D2_6_SIZE_SLICE_1_D), CEIL(size_e, D2_6_SIZE_SLICE_1_E), CEIL(size_f, D2_6_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}

#ifndef NWCHEM_TCE_CCSD_T
	// Copy the Result from Device to Host
	cudaMemcpy(t3, dev_t3, sizeof(double) * (size_a * size_b * size_c * size_d * size_e * size_f), cudaMemcpyDeviceToHost);

	// cudaFree()
    cudaFree(dev_t3);
#endif

    cudaFree(dev_t2);	cudaFree(dev_v2);
}

// created by tc_gen_definition_new()
#define D2_7_SIZE_SLICE_1_G 16
#define D2_7_SIZE_SLICE_1_A 16
#define D2_7_SIZE_SLICE_1_E 4
#define D2_7_SIZE_SLICE_1_F 1
#define D2_7_SIZE_SLICE_1_D 8
#define D2_7_SIZE_SLICE_1_C 8
#define D2_7_SIZE_SLICE_1_B 1

#define D2_7_SIZE_INT_UNIT_1 D2_7_SIZE_SLICE_1_G

#define D2_7_SIZE_TB_1_X 	D2_7_SIZE_SLICE_1_A * D2_7_SIZE_SLICE_1_F
#define D2_7_SIZE_TB_1_Y 	D2_7_SIZE_SLICE_1_D * D2_7_SIZE_SLICE_1_B
#define D2_7_SIZE_REG_1_X 	D2_7_SIZE_SLICE_1_E
#define D2_7_SIZE_REG_1_Y 	D2_7_SIZE_SLICE_1_C

#ifndef NWCHEM_TCE_CCSD_T
#define CEIL(a, b) 		(((a) + (b) - 1) / (b))
#endif

// created by tc_gen_code_Kernel()
__global__ void d2_7_kernel__1_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_7_SIZE_SLICE_1_A;
	int idx_f = threadIdx.x / D2_7_SIZE_SLICE_1_A;
	int idx_d = threadIdx.y % D2_7_SIZE_SLICE_1_D;
	int idx_b = threadIdx.y / D2_7_SIZE_SLICE_1_D;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_7_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_7_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_7_SIZE_SLICE_1_C + (blk_idx_d * D2_7_SIZE_SLICE_1_D + idx_d + (blk_idx_e * D2_7_SIZE_SLICE_1_E + (blk_idx_f * D2_7_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'y', 't2', ['g', 'd', 'c', 'b']], [16, 'STR_SD2_V2_H7', 'x', 'v2', ['g', 'a', 'e', 'f']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_7_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'd', 'c', 'b']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_d < rng_d
			sm_a[threadIdx.x][threadIdx.y + ll * 8] = dev_t2[(blk_idx_d * D2_7_SIZE_SLICE_1_D + idx_d + (blk_idx_c * D2_7_SIZE_SLICE_1_C + ll + (blk_idx_b * D2_7_SIZE_SLICE_1_B + 0) * size_c) * size_d) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'a', 'e', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_d + 0 < rng_a
			sm_b[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_v2[(blk_idx_a * D2_7_SIZE_SLICE_1_A + idx_d + 0 + (blk_idx_e * D2_7_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_7_SIZE_SLICE_1_F + 0) * size_e) * size_a) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_d + 8 < rng_a
			// Exception: Full-Full
			sm_b[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_v2[(blk_idx_a * D2_7_SIZE_SLICE_1_A + idx_d + 8 + (blk_idx_e * D2_7_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_7_SIZE_SLICE_1_F + 0) * size_e) * size_a) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_7_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 0];
			temp_bv[1] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 8];
			temp_bv[2] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 16];
			temp_bv[3] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 24];
			temp_bv[4] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 32];
			temp_bv[5] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 40];
			temp_bv[6] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 48];
			temp_bv[7] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_b[ll][idx_a + (idx_f) * D2_7_SIZE_SLICE_1_A + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_7_kernel__2_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_7_SIZE_SLICE_1_A;
	int idx_f = threadIdx.x / D2_7_SIZE_SLICE_1_A;
	int idx_d = threadIdx.y % D2_7_SIZE_SLICE_1_D;
	int idx_b = threadIdx.y / D2_7_SIZE_SLICE_1_D;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_7_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_7_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_7_SIZE_SLICE_1_C + (blk_idx_d * D2_7_SIZE_SLICE_1_D + idx_d + (blk_idx_e * D2_7_SIZE_SLICE_1_E + (blk_idx_f * D2_7_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'y', 't2', ['g', 'd', 'c', 'b']], [16, 'STR_SD2_V2_H7', 'x', 'v2', ['g', 'a', 'e', 'f']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_7_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_7_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_7_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'd', 'c', 'b']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_d < rng_d
			sm_a[threadIdx.x][threadIdx.y + ll * 8] = dev_t2[(blk_idx_d * D2_7_SIZE_SLICE_1_D + idx_d + (blk_idx_c * D2_7_SIZE_SLICE_1_C + ll + (blk_idx_b * D2_7_SIZE_SLICE_1_B + 0) * size_c) * size_d) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_7_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'a', 'e', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_d + 0 < rng_a
			sm_b[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_v2[(blk_idx_a * D2_7_SIZE_SLICE_1_A + idx_d + 0 + (blk_idx_e * D2_7_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_7_SIZE_SLICE_1_F + 0) * size_e) * size_a) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_d + 8 < rng_a
			if (threadIdx.x + l < size_internal) 
			sm_b[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_v2[(blk_idx_a * D2_7_SIZE_SLICE_1_A + idx_d + 8 + (blk_idx_e * D2_7_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_7_SIZE_SLICE_1_F + 0) * size_e) * size_a) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_7_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 0];
			temp_bv[1] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 8];
			temp_bv[2] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 16];
			temp_bv[3] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 24];
			temp_bv[4] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 32];
			temp_bv[5] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 40];
			temp_bv[6] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 48];
			temp_bv[7] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_b[ll][idx_a + (idx_f) * D2_7_SIZE_SLICE_1_A + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_7_kernel__3_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_7_SIZE_SLICE_1_A;
	int idx_f = threadIdx.x / D2_7_SIZE_SLICE_1_A;
	int idx_d = threadIdx.y % D2_7_SIZE_SLICE_1_D;
	int idx_b = threadIdx.y / D2_7_SIZE_SLICE_1_D;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_7_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_7_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_7_SIZE_SLICE_1_C + (blk_idx_d * D2_7_SIZE_SLICE_1_D + idx_d + (blk_idx_e * D2_7_SIZE_SLICE_1_E + (blk_idx_f * D2_7_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_7_SIZE_SLICE_1_A)) >= D2_7_SIZE_SLICE_1_A)
	{
		rng_a = D2_7_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_7_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_7_SIZE_SLICE_1_B)) >= D2_7_SIZE_SLICE_1_B)
	{
		rng_b = D2_7_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_7_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_7_SIZE_SLICE_1_C)) >= D2_7_SIZE_SLICE_1_C)
	{
		rng_c = D2_7_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_7_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_7_SIZE_SLICE_1_D)) >= D2_7_SIZE_SLICE_1_D)
	{
		rng_d = D2_7_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_7_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_7_SIZE_SLICE_1_E)) >= D2_7_SIZE_SLICE_1_E)
	{
		rng_e = D2_7_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_7_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_7_SIZE_SLICE_1_F)) >= D2_7_SIZE_SLICE_1_F)
	{
		rng_f = D2_7_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_7_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'y', 't2', ['g', 'd', 'c', 'b']], [16, 'STR_SD2_V2_H7', 'x', 'v2', ['g', 'a', 'e', 'f']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_7_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_d < rng_d && 0 < rng_b)
		for (int ll = 0; ll < rng_c; ll++)
		{
			// ['g', 'd', 'c', 'b']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_d < rng_d
			sm_a[threadIdx.x][threadIdx.y + ll * 8] = dev_t2[(blk_idx_d * D2_7_SIZE_SLICE_1_D + idx_d + (blk_idx_c * D2_7_SIZE_SLICE_1_C + ll + (blk_idx_b * D2_7_SIZE_SLICE_1_B + 0) * size_c) * size_d) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_d < rng_a && 0 < rng_f)
		for (int ll = 0; ll < rng_e; ll++)
		{
			// ['g', 'a', 'e', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_d + 0 < rng_a
			sm_b[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_v2[(blk_idx_a * D2_7_SIZE_SLICE_1_A + idx_d + 0 + (blk_idx_e * D2_7_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_7_SIZE_SLICE_1_F + 0) * size_e) * size_a) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_d + 8 < rng_a
			if (idx_d + 8 < rng_a) 
			sm_b[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_v2[(blk_idx_a * D2_7_SIZE_SLICE_1_A + idx_d + 8 + (blk_idx_e * D2_7_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_7_SIZE_SLICE_1_F + 0) * size_e) * size_a) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_7_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 0];
			temp_bv[1] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 8];
			temp_bv[2] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 16];
			temp_bv[3] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 24];
			temp_bv[4] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 32];
			temp_bv[5] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 40];
			temp_bv[6] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 48];
			temp_bv[7] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_b[ll][idx_a + (idx_f) * D2_7_SIZE_SLICE_1_A + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_f < rng_f && idx_d < rng_d && idx_b < rng_b)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_c && j < rng_e)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_7_kernel__4_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_7_SIZE_SLICE_1_A;
	int idx_f = threadIdx.x / D2_7_SIZE_SLICE_1_A;
	int idx_d = threadIdx.y % D2_7_SIZE_SLICE_1_D;
	int idx_b = threadIdx.y / D2_7_SIZE_SLICE_1_D;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_7_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_7_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_7_SIZE_SLICE_1_C + (blk_idx_d * D2_7_SIZE_SLICE_1_D + idx_d + (blk_idx_e * D2_7_SIZE_SLICE_1_E + (blk_idx_f * D2_7_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_7_SIZE_SLICE_1_A)) >= D2_7_SIZE_SLICE_1_A)
	{
		rng_a = D2_7_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_7_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_7_SIZE_SLICE_1_B)) >= D2_7_SIZE_SLICE_1_B)
	{
		rng_b = D2_7_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_7_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_7_SIZE_SLICE_1_C)) >= D2_7_SIZE_SLICE_1_C)
	{
		rng_c = D2_7_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_7_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_7_SIZE_SLICE_1_D)) >= D2_7_SIZE_SLICE_1_D)
	{
		rng_d = D2_7_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_7_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_7_SIZE_SLICE_1_E)) >= D2_7_SIZE_SLICE_1_E)
	{
		rng_e = D2_7_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_7_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_7_SIZE_SLICE_1_F)) >= D2_7_SIZE_SLICE_1_F)
	{
		rng_f = D2_7_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_7_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'y', 't2', ['g', 'd', 'c', 'b']], [16, 'STR_SD2_V2_H7', 'x', 'v2', ['g', 'a', 'e', 'f']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_7_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_7_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_d < rng_d && 0 < rng_b && threadIdx.x < D2_7_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_c; ll++)
		{
			// ['g', 'd', 'c', 'b']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_d < rng_d
			sm_a[threadIdx.x][threadIdx.y + ll * 8] = dev_t2[(blk_idx_d * D2_7_SIZE_SLICE_1_D + idx_d + (blk_idx_c * D2_7_SIZE_SLICE_1_C + ll + (blk_idx_b * D2_7_SIZE_SLICE_1_B + 0) * size_c) * size_d) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_d < rng_a && 0 < rng_f && threadIdx.x < D2_7_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_e; ll++)
		{
			// ['g', 'a', 'e', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_d + 0 < rng_a
			sm_b[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_v2[(blk_idx_a * D2_7_SIZE_SLICE_1_A + idx_d + 0 + (blk_idx_e * D2_7_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_7_SIZE_SLICE_1_F + 0) * size_e) * size_a) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_d + 8 < rng_a
			if (threadIdx.x + l < size_internal && idx_d + 8 < rng_a) 
			sm_b[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_v2[(blk_idx_a * D2_7_SIZE_SLICE_1_A + idx_d + 8 + (blk_idx_e * D2_7_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_7_SIZE_SLICE_1_F + 0) * size_e) * size_a) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_7_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 0];
			temp_bv[1] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 8];
			temp_bv[2] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 16];
			temp_bv[3] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 24];
			temp_bv[4] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 32];
			temp_bv[5] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 40];
			temp_bv[6] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 48];
			temp_bv[7] = sm_a[ll][idx_d + (idx_b) * D2_7_SIZE_SLICE_1_D + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_b[ll][idx_a + (idx_f) * D2_7_SIZE_SLICE_1_A + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_f < rng_f && idx_d < rng_d && idx_b < rng_b)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_c && j < rng_e)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// written by tc_interface.tc_gen_code_interface_Header()
void sd_t_d2_7_cogent(int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, double* t3, double* host_t2, double* host_v2, int cond_kernel_1, int opt_register_transpose)
{
	int num_thread_blocks_kernel_1;

	double* dev_t3;
	double* dev_t2;
	double* dev_v2;

	num_thread_blocks_kernel_1 = CEIL(size_a, D2_7_SIZE_SLICE_1_A) * CEIL(size_b, D2_7_SIZE_SLICE_1_B) * CEIL(size_c, D2_7_SIZE_SLICE_1_C) * CEIL(size_d, D2_7_SIZE_SLICE_1_D) * CEIL(size_e, D2_7_SIZE_SLICE_1_E) * CEIL(size_f, D2_7_SIZE_SLICE_1_F);

#ifdef NWCHEM_TCE_CCSD_T
    dev_t3 = t3_d;
#else
    // cudaMalloc()
    cudaMalloc((void**) &dev_t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f);
#endif
	cudaMalloc((void**) &dev_t2, sizeof(double) * size_b * size_c * size_d * size_g);
	cudaMalloc((void**) &dev_v2, sizeof(double) * size_f * size_e * size_a * size_g);

#ifndef NWCHEM_TCE_CCSD_T
	// cudaMemcpy()
	cudaMemcpy(dev_t3, t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f, cudaMemcpyHostToDevice);
#endif

    cudaMemcpy(dev_t2, host_t2, sizeof(double) * size_b * size_c * size_d * size_g, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2, host_v2, sizeof(double) * size_f * size_e * size_a * size_g, cudaMemcpyHostToDevice);

#ifndef NWCHEM_TCE_CCSD_T
	// Related to Kernels
	// There are 1 Basic Kernels
	long long int tmp_operations = 2 * (long long int)(size_a * size_b * size_c * size_d * size_e * size_f) * size_g;
    
    printf ("========================================= fusedKernels =============================================\n");
	printf ("		Grid Size  : %6d (1D)\n", num_thread_blocks_kernel_1);
	printf ("		Block-size : %2d, %2d (2D)\n", D2_7_SIZE_TB_1_X, D2_7_SIZE_TB_1_Y);
	printf ("		Reg.-size  : %2d, %2d (2D)\n", D2_7_SIZE_REG_1_X, D2_7_SIZE_REG_1_Y);
	printf ("		A thread deals with (%d x %d) elements (basically)\n", D2_7_SIZE_TB_1_X * D2_7_SIZE_REG_1_X, D2_7_SIZE_TB_1_Y * D2_7_SIZE_REG_1_Y);
	printf ("		# of Operations: %lld\n", tmp_operations);
	printf ("====================================================================================================\n");
#endif

    dim3 gridsize_1(num_thread_blocks_kernel_1);
	dim3 blocksize_1(D2_7_SIZE_TB_1_X, D2_7_SIZE_TB_1_Y);

	int stride_output_a = 1;
	int stride_output_b = stride_output_a * size_a;
	int stride_output_c = stride_output_b * size_b;
	int stride_output_d = stride_output_c * size_c;
	int stride_output_e = stride_output_d * size_d;
	int stride_output_f = stride_output_e * size_e;

	int stride_reg_x_1 = stride_output_e;
	int stride_reg_y_1 = stride_output_c;

	int size_internal = size_g;

	int stride_int_t2 = 1;
	int stride_int_v2 = 1;

	// Decision Tree for Kernel Types
	// No Chance to Utilize the Register Transpose
	if (size_a % D2_7_SIZE_SLICE_1_A == 0 && size_b % D2_7_SIZE_SLICE_1_B == 0 && size_c % D2_7_SIZE_SLICE_1_C == 0 && size_d % D2_7_SIZE_SLICE_1_D == 0 && size_e % D2_7_SIZE_SLICE_1_E == 0 && size_f % D2_7_SIZE_SLICE_1_F == 0)
	{
		// [2] Extenral Index: Full
		if (size_g % D2_7_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Full && Internal: Full
			//printf ("External: Full, Internal: Full\n");
			d2_7_kernel__1_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_7_SIZE_SLICE_1_A), CEIL(size_b, D2_7_SIZE_SLICE_1_B), CEIL(size_c, D2_7_SIZE_SLICE_1_C), CEIL(size_d, D2_7_SIZE_SLICE_1_D), CEIL(size_e, D2_7_SIZE_SLICE_1_E), CEIL(size_f, D2_7_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Full && Internal: Partial
			//printf ("External: Full, Internal: Partial\n");
			d2_7_kernel__2_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_7_SIZE_SLICE_1_A), CEIL(size_b, D2_7_SIZE_SLICE_1_B), CEIL(size_c, D2_7_SIZE_SLICE_1_C), CEIL(size_d, D2_7_SIZE_SLICE_1_D), CEIL(size_e, D2_7_SIZE_SLICE_1_E), CEIL(size_f, D2_7_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}
	else
	{
		// [2] Extenral Index: Partial
		if (size_g % D2_7_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Partial && Internal: Full
			//printf ("External: Partial, Internal: Full\n");
			d2_7_kernel__3_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_7_SIZE_SLICE_1_A), CEIL(size_b, D2_7_SIZE_SLICE_1_B), CEIL(size_c, D2_7_SIZE_SLICE_1_C), CEIL(size_d, D2_7_SIZE_SLICE_1_D), CEIL(size_e, D2_7_SIZE_SLICE_1_E), CEIL(size_f, D2_7_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Partial && Internal: Partial
			//printf ("External: Partial, Internal: Partial\n");
			d2_7_kernel__4_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_7_SIZE_SLICE_1_A), CEIL(size_b, D2_7_SIZE_SLICE_1_B), CEIL(size_c, D2_7_SIZE_SLICE_1_C), CEIL(size_d, D2_7_SIZE_SLICE_1_D), CEIL(size_e, D2_7_SIZE_SLICE_1_E), CEIL(size_f, D2_7_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}

#ifndef NWCHEM_TCE_CCSD_T
	// Copy the Result from Device to Host
	cudaMemcpy(t3, dev_t3, sizeof(double) * (size_a * size_b * size_c * size_d * size_e * size_f), cudaMemcpyDeviceToHost);

	// cudaFree()
    cudaFree(dev_t3);
#endif

    cudaFree(dev_t2);	cudaFree(dev_v2);
}

// created by tc_gen_definition_new()
#define D2_8_SIZE_SLICE_1_G 16
#define D2_8_SIZE_SLICE_1_A 16
#define D2_8_SIZE_SLICE_1_D 4
#define D2_8_SIZE_SLICE_1_B 1
#define D2_8_SIZE_SLICE_1_C 8
#define D2_8_SIZE_SLICE_1_E 8
#define D2_8_SIZE_SLICE_1_F 1

#define D2_8_SIZE_INT_UNIT_1 D2_8_SIZE_SLICE_1_G

#define D2_8_SIZE_TB_1_X 	D2_8_SIZE_SLICE_1_A * D2_8_SIZE_SLICE_1_B
#define D2_8_SIZE_TB_1_Y 	D2_8_SIZE_SLICE_1_C * D2_8_SIZE_SLICE_1_F
#define D2_8_SIZE_REG_1_X 	D2_8_SIZE_SLICE_1_D
#define D2_8_SIZE_REG_1_Y 	D2_8_SIZE_SLICE_1_E

// created by tc_gen_code_Kernel()
__global__ void d2_8_kernel__1_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_8_SIZE_SLICE_1_A;
	int idx_b = threadIdx.x / D2_8_SIZE_SLICE_1_A;
	int idx_c = threadIdx.y % D2_8_SIZE_SLICE_1_C;
	int idx_f = threadIdx.y / D2_8_SIZE_SLICE_1_C;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_8_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_8_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_8_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_8_SIZE_SLICE_1_D + (blk_idx_e * D2_8_SIZE_SLICE_1_E + (blk_idx_f * D2_8_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'd', 'b', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'c', 'e', 'f']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_8_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'd', 'b', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_d * D2_8_SIZE_SLICE_1_D + ll + (blk_idx_b * D2_8_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_8_SIZE_SLICE_1_A + idx_c + 0) * size_b) * size_d) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			// Exception: Full-Full
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_d * D2_8_SIZE_SLICE_1_D + ll + (blk_idx_b * D2_8_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_8_SIZE_SLICE_1_A + idx_c + 8) * size_b) * size_d) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'c', 'e', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_c < rng_c
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_c * D2_8_SIZE_SLICE_1_C + idx_c + (blk_idx_e * D2_8_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_8_SIZE_SLICE_1_F + 0) * size_e) * size_c) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_8_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 0];
			temp_bv[1] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 8];
			temp_bv[2] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 16];
			temp_bv[3] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 24];
			temp_bv[4] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 32];
			temp_bv[5] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 40];
			temp_bv[6] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 48];
			temp_bv[7] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_b + (idx_a) * D2_8_SIZE_SLICE_1_B + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_8_kernel__2_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_8_SIZE_SLICE_1_A;
	int idx_b = threadIdx.x / D2_8_SIZE_SLICE_1_A;
	int idx_c = threadIdx.y % D2_8_SIZE_SLICE_1_C;
	int idx_f = threadIdx.y / D2_8_SIZE_SLICE_1_C;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_8_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_8_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_8_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_8_SIZE_SLICE_1_D + (blk_idx_e * D2_8_SIZE_SLICE_1_E + (blk_idx_f * D2_8_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'd', 'b', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'c', 'e', 'f']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_8_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_8_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_8_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'd', 'b', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_d * D2_8_SIZE_SLICE_1_D + ll + (blk_idx_b * D2_8_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_8_SIZE_SLICE_1_A + idx_c + 0) * size_b) * size_d) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			if (threadIdx.x + l < size_internal) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_d * D2_8_SIZE_SLICE_1_D + ll + (blk_idx_b * D2_8_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_8_SIZE_SLICE_1_A + idx_c + 8) * size_b) * size_d) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_8_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'c', 'e', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_c < rng_c
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_c * D2_8_SIZE_SLICE_1_C + idx_c + (blk_idx_e * D2_8_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_8_SIZE_SLICE_1_F + 0) * size_e) * size_c) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_8_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 0];
			temp_bv[1] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 8];
			temp_bv[2] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 16];
			temp_bv[3] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 24];
			temp_bv[4] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 32];
			temp_bv[5] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 40];
			temp_bv[6] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 48];
			temp_bv[7] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_b + (idx_a) * D2_8_SIZE_SLICE_1_B + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_8_kernel__3_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_8_SIZE_SLICE_1_A;
	int idx_b = threadIdx.x / D2_8_SIZE_SLICE_1_A;
	int idx_c = threadIdx.y % D2_8_SIZE_SLICE_1_C;
	int idx_f = threadIdx.y / D2_8_SIZE_SLICE_1_C;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_8_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_8_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_8_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_8_SIZE_SLICE_1_D + (blk_idx_e * D2_8_SIZE_SLICE_1_E + (blk_idx_f * D2_8_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_8_SIZE_SLICE_1_A)) >= D2_8_SIZE_SLICE_1_A)
	{
		rng_a = D2_8_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_8_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_8_SIZE_SLICE_1_B)) >= D2_8_SIZE_SLICE_1_B)
	{
		rng_b = D2_8_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_8_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_8_SIZE_SLICE_1_C)) >= D2_8_SIZE_SLICE_1_C)
	{
		rng_c = D2_8_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_8_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_8_SIZE_SLICE_1_D)) >= D2_8_SIZE_SLICE_1_D)
	{
		rng_d = D2_8_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_8_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_8_SIZE_SLICE_1_E)) >= D2_8_SIZE_SLICE_1_E)
	{
		rng_e = D2_8_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_8_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_8_SIZE_SLICE_1_F)) >= D2_8_SIZE_SLICE_1_F)
	{
		rng_f = D2_8_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_8_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'd', 'b', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'c', 'e', 'f']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_8_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (0 < rng_b && idx_c < rng_a)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'd', 'b', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_d * D2_8_SIZE_SLICE_1_D + ll + (blk_idx_b * D2_8_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_8_SIZE_SLICE_1_A + idx_c + 0) * size_b) * size_d) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			if (0 < rng_b && idx_c + 8 < rng_a) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_d * D2_8_SIZE_SLICE_1_D + ll + (blk_idx_b * D2_8_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_8_SIZE_SLICE_1_A + idx_c + 8) * size_b) * size_d) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_c < rng_c && 0 < rng_f)
		for (int ll = 0; ll < rng_e; ll++)
		{
			// ['g', 'c', 'e', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_c < rng_c
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_c * D2_8_SIZE_SLICE_1_C + idx_c + (blk_idx_e * D2_8_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_8_SIZE_SLICE_1_F + 0) * size_e) * size_c) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_8_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 0];
			temp_bv[1] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 8];
			temp_bv[2] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 16];
			temp_bv[3] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 24];
			temp_bv[4] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 32];
			temp_bv[5] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 40];
			temp_bv[6] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 48];
			temp_bv[7] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_b + (idx_a) * D2_8_SIZE_SLICE_1_B + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_b < rng_b && idx_c < rng_c && idx_f < rng_f)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_e && j < rng_d)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_8_kernel__4_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_8_SIZE_SLICE_1_A;
	int idx_b = threadIdx.x / D2_8_SIZE_SLICE_1_A;
	int idx_c = threadIdx.y % D2_8_SIZE_SLICE_1_C;
	int idx_f = threadIdx.y / D2_8_SIZE_SLICE_1_C;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_8_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_8_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_8_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_8_SIZE_SLICE_1_D + (blk_idx_e * D2_8_SIZE_SLICE_1_E + (blk_idx_f * D2_8_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_8_SIZE_SLICE_1_A)) >= D2_8_SIZE_SLICE_1_A)
	{
		rng_a = D2_8_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_8_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_8_SIZE_SLICE_1_B)) >= D2_8_SIZE_SLICE_1_B)
	{
		rng_b = D2_8_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_8_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_8_SIZE_SLICE_1_C)) >= D2_8_SIZE_SLICE_1_C)
	{
		rng_c = D2_8_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_8_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_8_SIZE_SLICE_1_D)) >= D2_8_SIZE_SLICE_1_D)
	{
		rng_d = D2_8_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_8_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_8_SIZE_SLICE_1_E)) >= D2_8_SIZE_SLICE_1_E)
	{
		rng_e = D2_8_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_8_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_8_SIZE_SLICE_1_F)) >= D2_8_SIZE_SLICE_1_F)
	{
		rng_f = D2_8_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_8_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'd', 'b', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'c', 'e', 'f']], '-=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_8_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_8_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (0 < rng_b && idx_c < rng_a && threadIdx.x < D2_8_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'd', 'b', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_d * D2_8_SIZE_SLICE_1_D + ll + (blk_idx_b * D2_8_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_8_SIZE_SLICE_1_A + idx_c + 0) * size_b) * size_d) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_b
			if (threadIdx.x + l < size_internal && 0 < rng_b && idx_c + 8 < rng_a) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_d * D2_8_SIZE_SLICE_1_D + ll + (blk_idx_b * D2_8_SIZE_SLICE_1_B + 0 + (blk_idx_a * D2_8_SIZE_SLICE_1_A + idx_c + 8) * size_b) * size_d) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_c < rng_c && 0 < rng_f && threadIdx.x < D2_8_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_e; ll++)
		{
			// ['g', 'c', 'e', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_c < rng_c
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_c * D2_8_SIZE_SLICE_1_C + idx_c + (blk_idx_e * D2_8_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_8_SIZE_SLICE_1_F + 0) * size_e) * size_c) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_8_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 0];
			temp_bv[1] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 8];
			temp_bv[2] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 16];
			temp_bv[3] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 24];
			temp_bv[4] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 32];
			temp_bv[5] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 40];
			temp_bv[6] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 48];
			temp_bv[7] = sm_b[ll][idx_c + (idx_f) * D2_8_SIZE_SLICE_1_C + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_b + (idx_a) * D2_8_SIZE_SLICE_1_B + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
				reg_tile[4][xx] -= temp_av * temp_bv[4];
				reg_tile[5][xx] -= temp_av * temp_bv[5];
				reg_tile[6][xx] -= temp_av * temp_bv[6];
				reg_tile[7][xx] -= temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_b < rng_b && idx_c < rng_c && idx_f < rng_f)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_e && j < rng_d)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// written by tc_interface.tc_gen_code_interface_Header()
void sd_t_d2_8_cogent(int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, double* t3, double* host_t2, double* host_v2, int cond_kernel_1, int opt_register_transpose)
{
	int num_thread_blocks_kernel_1;

	double* dev_t3;
	double* dev_t2;
	double* dev_v2;

    num_thread_blocks_kernel_1 = CEIL(size_a, D2_8_SIZE_SLICE_1_A) * CEIL(size_b, D2_8_SIZE_SLICE_1_B) * CEIL(size_c, D2_8_SIZE_SLICE_1_C) * CEIL(size_d, D2_8_SIZE_SLICE_1_D) * CEIL(size_e, D2_8_SIZE_SLICE_1_E) * CEIL(size_f, D2_8_SIZE_SLICE_1_F);
    
#ifdef NWCHEM_TCE_CCSD_T
    dev_t3 = t3_d;
#else
	// cudaMalloc()
    cudaMalloc((void**) &dev_t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f);
#endif

	cudaMalloc((void**) &dev_t2, sizeof(double) * size_a * size_b * size_d * size_g);
	cudaMalloc((void**) &dev_v2, sizeof(double) * size_f * size_e * size_c * size_g);

#ifndef NWCHEM_TCE_CCSD_T
	// cudaMemcpy()
	cudaMemcpy(dev_t3, t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f, cudaMemcpyHostToDevice);
#endif

    cudaMemcpy(dev_t2, host_t2, sizeof(double) * size_a * size_b * size_d * size_g, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2, host_v2, sizeof(double) * size_f * size_e * size_c * size_g, cudaMemcpyHostToDevice);

#ifndef NWCHEM_TCE_CCSD_T
	// Related to Kernels
	// There are 1 Basic Kernels
    long long int tmp_operations = 2 * (long long int)(size_a * size_b * size_c * size_d * size_e * size_f) * size_g;

	printf ("========================================= fusedKernels =============================================\n");
	printf ("		Grid Size  : %6d (1D)\n", num_thread_blocks_kernel_1);
	printf ("		Block-size : %2d, %2d (2D)\n", D2_8_SIZE_TB_1_X, D2_8_SIZE_TB_1_Y);
	printf ("		Reg.-size  : %2d, %2d (2D)\n", D2_8_SIZE_REG_1_X, D2_8_SIZE_REG_1_Y);
	printf ("		A thread deals with (%d x %d) elements (basically)\n", D2_8_SIZE_TB_1_X * D2_8_SIZE_REG_1_X, D2_8_SIZE_TB_1_Y * D2_8_SIZE_REG_1_Y);
	printf ("		# of Operations: %lld\n", tmp_operations);
	printf ("====================================================================================================\n");
#endif

    dim3 gridsize_1(num_thread_blocks_kernel_1);
	dim3 blocksize_1(D2_8_SIZE_TB_1_X, D2_8_SIZE_TB_1_Y);

	int stride_output_a = 1;
	int stride_output_b = stride_output_a * size_a;
	int stride_output_c = stride_output_b * size_b;
	int stride_output_d = stride_output_c * size_c;
	int stride_output_e = stride_output_d * size_d;
	int stride_output_f = stride_output_e * size_e;

	int stride_reg_x_1 = stride_output_d;
	int stride_reg_y_1 = stride_output_e;

	int size_internal = size_g;

	int stride_int_t2 = 1;
	int stride_int_v2 = 1;

	// Decision Tree for Kernel Types
	// No Chance to Utilize the Register Transpose
	if (size_a % D2_8_SIZE_SLICE_1_A == 0 && size_b % D2_8_SIZE_SLICE_1_B == 0 && size_c % D2_8_SIZE_SLICE_1_C == 0 && size_d % D2_8_SIZE_SLICE_1_D == 0 && size_e % D2_8_SIZE_SLICE_1_E == 0 && size_f % D2_8_SIZE_SLICE_1_F == 0)
	{
		// [2] Extenral Index: Full
		if (size_g % D2_8_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Full && Internal: Full
			//printf ("External: Full, Internal: Full\n");
			d2_8_kernel__1_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_8_SIZE_SLICE_1_A), CEIL(size_b, D2_8_SIZE_SLICE_1_B), CEIL(size_c, D2_8_SIZE_SLICE_1_C), CEIL(size_d, D2_8_SIZE_SLICE_1_D), CEIL(size_e, D2_8_SIZE_SLICE_1_E), CEIL(size_f, D2_8_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Full && Internal: Partial
			//printf ("External: Full, Internal: Partial\n");
			d2_8_kernel__2_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_8_SIZE_SLICE_1_A), CEIL(size_b, D2_8_SIZE_SLICE_1_B), CEIL(size_c, D2_8_SIZE_SLICE_1_C), CEIL(size_d, D2_8_SIZE_SLICE_1_D), CEIL(size_e, D2_8_SIZE_SLICE_1_E), CEIL(size_f, D2_8_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}
	else
	{
		// [2] Extenral Index: Partial
		if (size_g % D2_8_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Partial && Internal: Full
			//printf ("External: Partial, Internal: Full\n");
			d2_8_kernel__3_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_8_SIZE_SLICE_1_A), CEIL(size_b, D2_8_SIZE_SLICE_1_B), CEIL(size_c, D2_8_SIZE_SLICE_1_C), CEIL(size_d, D2_8_SIZE_SLICE_1_D), CEIL(size_e, D2_8_SIZE_SLICE_1_E), CEIL(size_f, D2_8_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Partial && Internal: Partial
			//printf ("External: Partial, Internal: Partial\n");
			d2_8_kernel__4_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_8_SIZE_SLICE_1_A), CEIL(size_b, D2_8_SIZE_SLICE_1_B), CEIL(size_c, D2_8_SIZE_SLICE_1_C), CEIL(size_d, D2_8_SIZE_SLICE_1_D), CEIL(size_e, D2_8_SIZE_SLICE_1_E), CEIL(size_f, D2_8_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}

#ifndef NWCHEM_TCE_CCSD_T
	// Copy the Result from Device to Host
	cudaMemcpy(t3, dev_t3, sizeof(double) * (size_a * size_b * size_c * size_d * size_e * size_f), cudaMemcpyDeviceToHost);

	// cudaFree()
    cudaFree(dev_t3);
#endif

    cudaFree(dev_t2);	cudaFree(dev_v2);
}

// created by tc_gen_definition_new()
#define D2_9_SIZE_SLICE_1_G 16
#define D2_9_SIZE_SLICE_1_A 16
#define D2_9_SIZE_SLICE_1_D 4
#define D2_9_SIZE_SLICE_1_C 1
#define D2_9_SIZE_SLICE_1_B 8
#define D2_9_SIZE_SLICE_1_E 8
#define D2_9_SIZE_SLICE_1_F 1

#define D2_9_SIZE_INT_UNIT_1 D2_9_SIZE_SLICE_1_G

#define D2_9_SIZE_TB_1_X 	D2_9_SIZE_SLICE_1_A * D2_9_SIZE_SLICE_1_C
#define D2_9_SIZE_TB_1_Y 	D2_9_SIZE_SLICE_1_B * D2_9_SIZE_SLICE_1_F
#define D2_9_SIZE_REG_1_X 	D2_9_SIZE_SLICE_1_D
#define D2_9_SIZE_REG_1_Y 	D2_9_SIZE_SLICE_1_E

// created by tc_gen_code_Kernel()
__global__ void d2_9_kernel__1_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_9_SIZE_SLICE_1_A;
	int idx_c = threadIdx.x / D2_9_SIZE_SLICE_1_A;
	int idx_b = threadIdx.y % D2_9_SIZE_SLICE_1_B;
	int idx_f = threadIdx.y / D2_9_SIZE_SLICE_1_B;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_9_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_9_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_9_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_9_SIZE_SLICE_1_D + (blk_idx_e * D2_9_SIZE_SLICE_1_E + (blk_idx_f * D2_9_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'd', 'c', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'b', 'e', 'f']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_9_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'd', 'c', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_d * D2_9_SIZE_SLICE_1_D + ll + (blk_idx_c * D2_9_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_9_SIZE_SLICE_1_A + idx_b + 0) * size_c) * size_d) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			// Exception: Full-Full
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_d * D2_9_SIZE_SLICE_1_D + ll + (blk_idx_c * D2_9_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_9_SIZE_SLICE_1_A + idx_b + 8) * size_c) * size_d) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		// No Need to Put Boundary-Checks before For-Statement: : 
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'b', 'e', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_b < rng_b
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_b * D2_9_SIZE_SLICE_1_B + idx_b + (blk_idx_e * D2_9_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_9_SIZE_SLICE_1_F + 0) * size_e) * size_b) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_9_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 0];
			temp_bv[1] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 8];
			temp_bv[2] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 16];
			temp_bv[3] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 24];
			temp_bv[4] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 32];
			temp_bv[5] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 40];
			temp_bv[6] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 48];
			temp_bv[7] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_c + (idx_a) * D2_9_SIZE_SLICE_1_C + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_9_kernel__2_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_9_SIZE_SLICE_1_A;
	int idx_c = threadIdx.x / D2_9_SIZE_SLICE_1_A;
	int idx_b = threadIdx.y % D2_9_SIZE_SLICE_1_B;
	int idx_f = threadIdx.y / D2_9_SIZE_SLICE_1_B;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_9_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_9_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_9_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_9_SIZE_SLICE_1_D + (blk_idx_e * D2_9_SIZE_SLICE_1_E + (blk_idx_f * D2_9_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;


	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'd', 'c', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'b', 'e', 'f']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_9_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_9_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_9_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			// ['g', 'd', 'c', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_d * D2_9_SIZE_SLICE_1_D + ll + (blk_idx_c * D2_9_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_9_SIZE_SLICE_1_A + idx_b + 0) * size_c) * size_d) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			if (threadIdx.x + l < size_internal) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_d * D2_9_SIZE_SLICE_1_D + ll + (blk_idx_c * D2_9_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_9_SIZE_SLICE_1_A + idx_b + 8) * size_c) * size_d) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (threadIdx.x < D2_9_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < 8; ll++)
		{
			// ['g', 'b', 'e', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_b < rng_b
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_b * D2_9_SIZE_SLICE_1_B + idx_b + (blk_idx_e * D2_9_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_9_SIZE_SLICE_1_F + 0) * size_e) * size_b) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_9_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 0];
			temp_bv[1] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 8];
			temp_bv[2] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 16];
			temp_bv[3] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 24];
			temp_bv[4] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 32];
			temp_bv[5] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 40];
			temp_bv[6] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 48];
			temp_bv[7] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_c + (idx_a) * D2_9_SIZE_SLICE_1_C + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	#pragma unroll 8
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_9_kernel__3_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_9_SIZE_SLICE_1_A;
	int idx_c = threadIdx.x / D2_9_SIZE_SLICE_1_A;
	int idx_b = threadIdx.y % D2_9_SIZE_SLICE_1_B;
	int idx_f = threadIdx.y / D2_9_SIZE_SLICE_1_B;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_9_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_9_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_9_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_9_SIZE_SLICE_1_D + (blk_idx_e * D2_9_SIZE_SLICE_1_E + (blk_idx_f * D2_9_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_9_SIZE_SLICE_1_A)) >= D2_9_SIZE_SLICE_1_A)
	{
		rng_a = D2_9_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_9_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_9_SIZE_SLICE_1_B)) >= D2_9_SIZE_SLICE_1_B)
	{
		rng_b = D2_9_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_9_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_9_SIZE_SLICE_1_C)) >= D2_9_SIZE_SLICE_1_C)
	{
		rng_c = D2_9_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_9_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_9_SIZE_SLICE_1_D)) >= D2_9_SIZE_SLICE_1_D)
	{
		rng_d = D2_9_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_9_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_9_SIZE_SLICE_1_E)) >= D2_9_SIZE_SLICE_1_E)
	{
		rng_e = D2_9_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_9_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_9_SIZE_SLICE_1_F)) >= D2_9_SIZE_SLICE_1_F)
	{
		rng_f = D2_9_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_9_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'd', 'c', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'b', 'e', 'f']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_9_SIZE_INT_UNIT_1)
	{
		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (0 < rng_c && idx_b < rng_a)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'd', 'c', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_d * D2_9_SIZE_SLICE_1_D + ll + (blk_idx_c * D2_9_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_9_SIZE_SLICE_1_A + idx_b + 0) * size_c) * size_d) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			if (0 < rng_c && idx_b + 8 < rng_a) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_d * D2_9_SIZE_SLICE_1_D + ll + (blk_idx_c * D2_9_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_9_SIZE_SLICE_1_A + idx_b + 8) * size_c) * size_d) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_b < rng_b && 0 < rng_f)
		for (int ll = 0; ll < rng_e; ll++)
		{
			// ['g', 'b', 'e', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_b < rng_b
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_b * D2_9_SIZE_SLICE_1_B + idx_b + (blk_idx_e * D2_9_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_9_SIZE_SLICE_1_F + 0) * size_e) * size_b) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_9_SIZE_INT_UNIT_1; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 0];
			temp_bv[1] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 8];
			temp_bv[2] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 16];
			temp_bv[3] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 24];
			temp_bv[4] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 32];
			temp_bv[5] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 40];
			temp_bv[6] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 48];
			temp_bv[7] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_c + (idx_a) * D2_9_SIZE_SLICE_1_C + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_c < rng_c && idx_b < rng_b && idx_f < rng_f)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_e && j < rng_d)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void d2_9_kernel__4_1(double* dev_t3, double* dev_t2, double* dev_v2, int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, int numBlk_a, int numBlk_b, int numBlk_c, int numBlk_d, int numBlk_e, int numBlk_f, int stride_int_t2, int stride_int_v2, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64];
	__shared__ double sm_b[16][64];


	int internal_upperbound   = 0;
	int internal_offset;

	// when opt_pre_computed == -1, all indices will be calculated manually
	// # of indices mapped on TB_X: 2
	// # of indices mapped on TB_Y: 2
	int idx_a = threadIdx.x % D2_9_SIZE_SLICE_1_A;
	int idx_c = threadIdx.x / D2_9_SIZE_SLICE_1_A;
	int idx_b = threadIdx.y % D2_9_SIZE_SLICE_1_B;
	int idx_f = threadIdx.y / D2_9_SIZE_SLICE_1_B;

	int tmp_blkIdx;
	int blk_idx_f = blockIdx.x / (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = blockIdx.x % (numBlk_e * numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_e = tmp_blkIdx / (numBlk_d * numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_d * numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_d = tmp_blkIdx / (numBlk_c * numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_c * numBlk_b * numBlk_a);

	int blk_idx_c = tmp_blkIdx / (numBlk_b * numBlk_a);
	tmp_blkIdx = tmp_blkIdx % (numBlk_b * numBlk_a);

	int blk_idx_b = tmp_blkIdx / numBlk_a;
	tmp_blkIdx = tmp_blkIdx % (numBlk_a);

	int  blk_idx_a = tmp_blkIdx;

	int t3_base_thread = blk_idx_a * D2_9_SIZE_SLICE_1_A + idx_a + (blk_idx_b * D2_9_SIZE_SLICE_1_B + idx_b + (blk_idx_c * D2_9_SIZE_SLICE_1_C + idx_c + (blk_idx_d * D2_9_SIZE_SLICE_1_D + (blk_idx_e * D2_9_SIZE_SLICE_1_E + (blk_idx_f * D2_9_SIZE_SLICE_1_F + idx_f) * size_e) * size_d) * size_c) * size_b) * size_a;

	// need to support partial tiles
	int rng_a, rng_b, rng_c, rng_d, rng_e, rng_f;
	if ((size_a - (blk_idx_a * D2_9_SIZE_SLICE_1_A)) >= D2_9_SIZE_SLICE_1_A)
	{
		rng_a = D2_9_SIZE_SLICE_1_A;
	}
	else
	{
		rng_a = size_a % D2_9_SIZE_SLICE_1_A;
	}
	if ((size_b - (blk_idx_b * D2_9_SIZE_SLICE_1_B)) >= D2_9_SIZE_SLICE_1_B)
	{
		rng_b = D2_9_SIZE_SLICE_1_B;
	}
	else
	{
		rng_b = size_b % D2_9_SIZE_SLICE_1_B;
	}
	if ((size_c - (blk_idx_c * D2_9_SIZE_SLICE_1_C)) >= D2_9_SIZE_SLICE_1_C)
	{
		rng_c = D2_9_SIZE_SLICE_1_C;
	}
	else
	{
		rng_c = size_c % D2_9_SIZE_SLICE_1_C;
	}
	if ((size_d - (blk_idx_d * D2_9_SIZE_SLICE_1_D)) >= D2_9_SIZE_SLICE_1_D)
	{
		rng_d = D2_9_SIZE_SLICE_1_D;
	}
	else
	{
		rng_d = size_d % D2_9_SIZE_SLICE_1_D;
	}
	if ((size_e - (blk_idx_e * D2_9_SIZE_SLICE_1_E)) >= D2_9_SIZE_SLICE_1_E)
	{
		rng_e = D2_9_SIZE_SLICE_1_E;
	}
	else
	{
		rng_e = size_e % D2_9_SIZE_SLICE_1_E;
	}
	if ((size_f - (blk_idx_f * D2_9_SIZE_SLICE_1_F)) >= D2_9_SIZE_SLICE_1_F)
	{
		rng_f = D2_9_SIZE_SLICE_1_F;
	}
	else
	{
		rng_f = size_f % D2_9_SIZE_SLICE_1_F;
	}

	double temp_av;
	double temp_bv[8];
	double reg_tile[8][4];

	for (int i = 0; i < 8; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction: [[16, 'STR_SD2_T2_H7', 'x', 't2', ['g', 'd', 'c', 'a']], [16, 'STR_SD2_V2_H7', 'y', 'v2', ['g', 'b', 'e', 'f']], '+=']
	#pragma unroll 1
	for (int l = 0; l < size_internal; l += D2_9_SIZE_INT_UNIT_1)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + D2_9_SIZE_INT_UNIT_1) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		//---------------------------------------------------------------------------------------------------
		// This is for the new version
		// This Part is for Loading Input-Left
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (0 < rng_c && idx_b < rng_a && threadIdx.x < D2_9_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_d; ll++)
		{
			// ['g', 'd', 'c', 'a']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			sm_a[threadIdx.x][threadIdx.y + 0 + ll * 16] = dev_t2[(blk_idx_d * D2_9_SIZE_SLICE_1_D + ll + (blk_idx_c * D2_9_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_9_SIZE_SLICE_1_A + idx_b + 0) * size_c) * size_d) * size_g + (threadIdx.x + l)];
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: 0 < rng_c
			if (threadIdx.x + l < size_internal && 0 < rng_c && idx_b + 8 < rng_a) 
			sm_a[threadIdx.x][threadIdx.y + 8 + ll * 16] = dev_t2[(blk_idx_d * D2_9_SIZE_SLICE_1_D + ll + (blk_idx_c * D2_9_SIZE_SLICE_1_C + 0 + (blk_idx_a * D2_9_SIZE_SLICE_1_A + idx_b + 8) * size_c) * size_d) * size_g + (threadIdx.x + l)];
		}
		
		// This Part is for Loading Input-Right
		// tc_gen_code_Kernel_Load_Inputs_Abstracts()
		if (idx_b < rng_b && 0 < rng_f && threadIdx.x < D2_9_SIZE_INT_UNIT_1 - internal_upperbound)
		for (int ll = 0; ll < rng_e; ll++)
		{
			// ['g', 'b', 'e', 'f']
			// Exception: Temp. version!: threadIdx.x + l
			// Exception: Temp. version!: idx_b < rng_b
			sm_b[threadIdx.x][threadIdx.y + ll * 8] = dev_v2[(blk_idx_b * D2_9_SIZE_SLICE_1_B + idx_b + (blk_idx_e * D2_9_SIZE_SLICE_1_E + ll + (blk_idx_f * D2_9_SIZE_SLICE_1_F + 0) * size_e) * size_b) * size_g + (threadIdx.x + l)];
		}
		__syncthreads();
		//---------------------------------------------------------------------------------------------------
		

		// Part: Generalized Threads
		for (int ll = 0; ll < D2_9_SIZE_INT_UNIT_1 - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 0];
			temp_bv[1] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 8];
			temp_bv[2] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 16];
			temp_bv[3] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 24];
			temp_bv[4] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 32];
			temp_bv[5] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 40];
			temp_bv[6] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 48];
			temp_bv[7] = sm_b[ll][idx_b + (idx_f) * D2_9_SIZE_SLICE_1_B + 56];

			for (int xx = 0; xx < 4; xx++) // (1)
			{
				temp_av = sm_a[ll][idx_c + (idx_a) * D2_9_SIZE_SLICE_1_C + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
				reg_tile[4][xx] += temp_av * temp_bv[4];
				reg_tile[5][xx] += temp_av * temp_bv[5];
				reg_tile[6][xx] += temp_av * temp_bv[6];
				reg_tile[7][xx] += temp_av * temp_bv[7];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_a < rng_a && idx_c < rng_c && idx_b < rng_b && idx_f < rng_f)
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_e && j < rng_d)
			{
			    dev_t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
			}
		}
	}
}

// written by tc_interface.tc_gen_code_interface_Header()
void sd_t_d2_9_cogent(int size_a, int size_b, int size_c, int size_d, int size_e, int size_f, int size_g, double* t3, double* host_t2, double* host_v2, int cond_kernel_1, int opt_register_transpose)
{
	int num_thread_blocks_kernel_1;

	double* dev_t3;
	double* dev_t2;
	double* dev_v2;

    num_thread_blocks_kernel_1 = CEIL(size_a, D2_9_SIZE_SLICE_1_A) * CEIL(size_b, D2_9_SIZE_SLICE_1_B) * CEIL(size_c, D2_9_SIZE_SLICE_1_C) * CEIL(size_d, D2_9_SIZE_SLICE_1_D) * CEIL(size_e, D2_9_SIZE_SLICE_1_E) * CEIL(size_f, D2_9_SIZE_SLICE_1_F);
    
#ifdef NWCHEM_TCE_CCSD_T
    dev_t3 = t3_d;
#else
	// cudaMalloc()
    cudaMalloc((void**) &dev_t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f);
#endif    

	cudaMalloc((void**) &dev_t2, sizeof(double) * size_a * size_c * size_d * size_g);
	cudaMalloc((void**) &dev_v2, sizeof(double) * size_f * size_e * size_b * size_g);

#ifndef NWCHEM_TCE_CCSD_T
	// cudaMemcpy()
    cudaMemcpy(dev_t3, t3, sizeof(double) * size_a * size_b * size_c * size_d * size_e * size_f, cudaMemcpyHostToDevice);
#endif

	cudaMemcpy(dev_t2, host_t2, sizeof(double) * size_a * size_c * size_d * size_g, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2, host_v2, sizeof(double) * size_f * size_e * size_b * size_g, cudaMemcpyHostToDevice);

#ifndef NWCHEM_TCE_CCSD_T
	// Related to Kernels
	// There are 1 Basic Kernels
    long long int tmp_operations = 2 * (long long int)(size_a * size_b * size_c * size_d * size_e * size_f) * size_g;

	printf ("========================================= fusedKernels =============================================\n");
	printf ("		Grid Size  : %6d (1D)\n", num_thread_blocks_kernel_1);
	printf ("		Block-size : %2d, %2d (2D)\n", D2_9_SIZE_TB_1_X, D2_9_SIZE_TB_1_Y);
	printf ("		Reg.-size  : %2d, %2d (2D)\n", D2_9_SIZE_REG_1_X, D2_9_SIZE_REG_1_Y);
	printf ("		A thread deals with (%d x %d) elements (basically)\n", D2_9_SIZE_TB_1_X * D2_9_SIZE_REG_1_X, D2_9_SIZE_TB_1_Y * D2_9_SIZE_REG_1_Y);
	printf ("		# of Operations: %lld\n", tmp_operations);
	printf ("====================================================================================================\n");
#endif

    dim3 gridsize_1(num_thread_blocks_kernel_1);
	dim3 blocksize_1(D2_9_SIZE_TB_1_X, D2_9_SIZE_TB_1_Y);

	int stride_output_a = 1;
	int stride_output_b = stride_output_a * size_a;
	int stride_output_c = stride_output_b * size_b;
	int stride_output_d = stride_output_c * size_c;
	int stride_output_e = stride_output_d * size_d;
	int stride_output_f = stride_output_e * size_e;

	int stride_reg_x_1 = stride_output_d;
	int stride_reg_y_1 = stride_output_e;

	int size_internal = size_g;

	int stride_int_t2 = 1;
	int stride_int_v2 = 1;

	// Decision Tree for Kernel Types
	// No Chance to Utilize the Register Transpose
	if (size_a % D2_9_SIZE_SLICE_1_A == 0 && size_b % D2_9_SIZE_SLICE_1_B == 0 && size_c % D2_9_SIZE_SLICE_1_C == 0 && size_d % D2_9_SIZE_SLICE_1_D == 0 && size_e % D2_9_SIZE_SLICE_1_E == 0 && size_f % D2_9_SIZE_SLICE_1_F == 0)
	{
		// [2] Extenral Index: Full
		if (size_g % D2_9_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Full && Internal: Full
			//printf ("External: Full, Internal: Full\n");
			d2_9_kernel__1_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_9_SIZE_SLICE_1_A), CEIL(size_b, D2_9_SIZE_SLICE_1_B), CEIL(size_c, D2_9_SIZE_SLICE_1_C), CEIL(size_d, D2_9_SIZE_SLICE_1_D), CEIL(size_e, D2_9_SIZE_SLICE_1_E), CEIL(size_f, D2_9_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Full && Internal: Partial
			//printf ("External: Full, Internal: Partial\n");
			d2_9_kernel__2_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_9_SIZE_SLICE_1_A), CEIL(size_b, D2_9_SIZE_SLICE_1_B), CEIL(size_c, D2_9_SIZE_SLICE_1_C), CEIL(size_d, D2_9_SIZE_SLICE_1_D), CEIL(size_e, D2_9_SIZE_SLICE_1_E), CEIL(size_f, D2_9_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
	}
	else
	{
		// [2] Extenral Index: Partial
		if (size_g % D2_9_SIZE_SLICE_1_G == 0)
		{
			// [3] Internal Index: Full
			// >>> External: Partial && Internal: Full
			//printf ("External: Partial, Internal: Full\n");
			d2_9_kernel__3_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_9_SIZE_SLICE_1_A), CEIL(size_b, D2_9_SIZE_SLICE_1_B), CEIL(size_c, D2_9_SIZE_SLICE_1_C), CEIL(size_d, D2_9_SIZE_SLICE_1_D), CEIL(size_e, D2_9_SIZE_SLICE_1_E), CEIL(size_f, D2_9_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
		else
		{
			// [4] Internal Index: Partial
			// >>> External: Partial && Internal: Partial
			//printf ("External: Partial, Internal: Partial\n");
			d2_9_kernel__4_1<<<gridsize_1, blocksize_1>>>(dev_t3, dev_t2, dev_v2, size_a, size_b, size_c, size_d, size_e, size_f, size_g, CEIL(size_a, D2_9_SIZE_SLICE_1_A), CEIL(size_b, D2_9_SIZE_SLICE_1_B), CEIL(size_c, D2_9_SIZE_SLICE_1_C), CEIL(size_d, D2_9_SIZE_SLICE_1_D), CEIL(size_e, D2_9_SIZE_SLICE_1_E), CEIL(size_f, D2_9_SIZE_SLICE_1_F), stride_int_t2, stride_int_v2, stride_reg_x_1, stride_reg_y_1, size_internal);
		}
    }
    
#ifndef NWCHEM_TCE_CCSD_T
	// Copy the Result from Device to Host
	cudaMemcpy(t3, dev_t3, sizeof(double) * (size_a * size_b * size_c * size_d * size_e * size_f), cudaMemcpyDeviceToHost);

	// cudaFree()
    cudaFree(dev_t3);	
#endif

    cudaFree(dev_t2);	cudaFree(dev_v2);
}

/*
 *      Partiall Fused Kernels and Fully Fused Kernels
 */
#define SIZE_SLICE_1_H3 4
#define SIZE_SLICE_1_H2 4
#define SIZE_SLICE_1_H1 4
#define SIZE_SLICE_1_P6 4
#define SIZE_SLICE_1_P5 4
#define SIZE_SLICE_1_P4 4
#define SIZE_SLICE_1_P7 16

#define SIZE_SLICE_2_H3 4
#define SIZE_SLICE_2_H2 4
#define SIZE_SLICE_2_H1 4
#define SIZE_SLICE_2_P6 4
#define SIZE_SLICE_2_P5 4
#define SIZE_SLICE_2_P4 4
#define SIZE_SLICE_2_P7 16

#define SIZE_INT_UNIT 	SIZE_SLICE_1_P7

#define SIZE_TB_1_X 	SIZE_SLICE_1_H3 * SIZE_SLICE_1_H2
#define SIZE_TB_1_Y 	SIZE_SLICE_1_P6 * SIZE_SLICE_1_H1
#define SIZE_REG_1_X 	SIZE_SLICE_1_P5
#define SIZE_REG_1_Y 	SIZE_SLICE_1_P4

#define SIZE_TB_2_X 	SIZE_SLICE_2_H3 * SIZE_SLICE_2_H2
#define SIZE_TB_2_Y 	SIZE_SLICE_2_P4 * SIZE_SLICE_2_H1
#define SIZE_REG_2_X 	SIZE_SLICE_2_P5
#define SIZE_REG_2_Y 	SIZE_SLICE_2_P6

#define NUM_INDEX 		6

// created by tc_gen_code_Kernel()
__global__ void kernel_ccsdT_sd2_fully_fused_partial_partial(double* t3, double* d_t2_1, double* d_t2_2, double* d_t2_3, double* d_t2_4, double* d_t2_5, double* d_t2_6, double* d_v2_1, double* d_v2_2, double* d_v2_3, double* d_v2_4, double* d_v2_5, double* d_v2_6, double* d_t2_7, double* d_t2_8, double* d_t2_9, double* d_v2_7, double* d_v2_8, double* d_v2_9, int size_h3, int size_h2, int size_h1, int size_p6, int size_p5, int size_p4, int size_p7, int numBlk_h3, int numBlk_h2, int numBlk_h1, int numBlk_p6, int numBlk_p5, int numBlk_p4, int kernel_1, int kernel_2, int kernel_3, int kernel_4, int kernel_5, int kernel_6, int kernel_7, int kernel_8, int kernel_9, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int internal_upperbound   	= 0;
	int internal_offset;

	// should support for non-full tiles
	int idx_h3 = threadIdx.x % SIZE_SLICE_1_H3;
	int idx_h2 = threadIdx.x / SIZE_SLICE_1_H3;
	int idx_p6 = threadIdx.y % SIZE_SLICE_1_P6;
	int idx_h1 = threadIdx.y / SIZE_SLICE_1_P6;

    // Common for Threads within a Thread Block
	int tmp_blkIdx;        
    int blk_idx_p4  = blockIdx.x / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);
    tmp_blkIdx      = blockIdx.x % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);

    int blk_idx_p5  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);

    int blk_idx_p6  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1);

    int blk_idx_h1 = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2);
    tmp_blkIdx     = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2);

    int blk_idx_h2 = (tmp_blkIdx) / (numBlk_h3);
    int blk_idx_h3 = (tmp_blkIdx) % (numBlk_h3);

    int rng_h3, rng_h2, rng_h1, rng_p6, rng_p5, rng_p4;

    if ((size_h3 - (blk_idx_h3 * SIZE_SLICE_2_H3)) >= SIZE_SLICE_2_H3)
    {
        rng_h3 = SIZE_SLICE_2_H3;
    }
    else
    {
        rng_h3 = size_h3 % SIZE_SLICE_2_H3;
    }
    
    if ((size_h2 - (blk_idx_h2 * SIZE_SLICE_2_H2)) >= SIZE_SLICE_2_H2)
    {
        rng_h2 = SIZE_SLICE_2_H2;
    }
    else
    {
        rng_h2 = size_h2 % SIZE_SLICE_2_H2;
    }

    if ((size_h1 - (blk_idx_h1 * SIZE_SLICE_2_H1)) >= SIZE_SLICE_2_H1)
    {
        rng_h1 = SIZE_SLICE_2_H1;
    }
    else
    {
        rng_h1 = size_h1 % SIZE_SLICE_2_H1;
    }
    
    if ((size_p6 - (blk_idx_p6 * SIZE_SLICE_2_P6)) >= SIZE_SLICE_2_P6)
    {
        rng_p6 = SIZE_SLICE_2_P6;
    }
    else
    {
        rng_p6 = size_p6 % SIZE_SLICE_2_P6;
    }

    if ((size_p5 - (blk_idx_p5 * SIZE_SLICE_2_P5)) >= SIZE_SLICE_2_P5)
    {
        rng_p5 = SIZE_SLICE_2_P5;
    }
    else
    {
        rng_p5 = size_p5 % SIZE_SLICE_2_P5;
    }

    if ((size_p4 - (blk_idx_p4 * SIZE_SLICE_2_P4)) >= SIZE_SLICE_2_P4)
    {
        rng_p4 = SIZE_SLICE_2_P4;
    }
    else
    {
        rng_p4 = size_p4 % SIZE_SLICE_2_P4;
    }

    int t3_base_thread = blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h3 + 
                        (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h2 + 
                        (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_h1 + 
                        (blk_idx_p6 * SIZE_SLICE_2_P6 + idx_p6 + 
                        (blk_idx_p5 * SIZE_SLICE_2_P5 + 
                        (blk_idx_p4 * SIZE_SLICE_2_P4) * size_p5) * size_p6) * size_h1) * size_h2) * size_h3;



	double temp_av;
	double temp_bv[4];
	double reg_tile[4][4];

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	//
	//	sd2_7, 8 and 9
	//
	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_7 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h1 && idx_h1 < rng_h2 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p6; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_7[(blk_idx_p6 *  SIZE_SLICE_2_P6 + ll + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h1) * size_h1) * size_p6) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h3 && idx_h1 < rng_p4 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_7[(blk_idx_h3 * SIZE_SLICE_2_H3 + idx_p6 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_h1) * size_p5) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_2_H3 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_8 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h2 && idx_h1 < rng_h3 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p6; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_8[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h1) * size_h2) * size_p6) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h1 && idx_h1 < rng_p4 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_8[(blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p6 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_1_P4 + idx_h1) * size_p5) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 0];
            temp_bv[1] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 16];
            temp_bv[2] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 32];
            temp_bv[3] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_2_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_9 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > rng_p6) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h1 && idx_h1 < rng_h3 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p6; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_9[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h1) * size_h1) * size_p6) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h2 && idx_h1 < rng_p4 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_9[(blk_idx_h2 * SIZE_SLICE_2_H2 + idx_p6 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_h1) * size_p5) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_2_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	//
	//
	//
	if (threadIdx.y < 4)						// 0, 1, 2, 3
	{
		// sm_a[16][64] <-- (4 x 16) x (4 x 4) = (16 x 64)		  'y''x'
		sm_a[0 + threadIdx.y * 4][threadIdx.x] 			= reg_tile[0][0];
		sm_a[1 + threadIdx.y * 4][threadIdx.x] 			= reg_tile[1][0];
		sm_a[2 + threadIdx.y * 4][threadIdx.x] 			= reg_tile[2][0];
		sm_a[3 + threadIdx.y * 4][threadIdx.x] 			= reg_tile[3][0];
		
		sm_a[0 + threadIdx.y * 4][threadIdx.x + 16] 	= reg_tile[0][1];
		sm_a[1 + threadIdx.y * 4][threadIdx.x + 16] 	= reg_tile[1][1];
		sm_a[2 + threadIdx.y * 4][threadIdx.x + 16] 	= reg_tile[2][1];
		sm_a[3 + threadIdx.y * 4][threadIdx.x + 16] 	= reg_tile[3][1];
		         
		sm_a[0 + threadIdx.y * 4][threadIdx.x + 32] 	= reg_tile[0][2];
		sm_a[1 + threadIdx.y * 4][threadIdx.x + 32] 	= reg_tile[1][2];
		sm_a[2 + threadIdx.y * 4][threadIdx.x + 32] 	= reg_tile[2][2];
		sm_a[3 + threadIdx.y * 4][threadIdx.x + 32] 	= reg_tile[3][2];
		
		sm_a[0 + threadIdx.y * 4][threadIdx.x + 48] 	= reg_tile[0][3];
		sm_a[1 + threadIdx.y * 4][threadIdx.x + 48] 	= reg_tile[1][3];
		sm_a[2 + threadIdx.y * 4][threadIdx.x + 48] 	= reg_tile[2][3];
		sm_a[3 + threadIdx.y * 4][threadIdx.x + 48] 	= reg_tile[3][3];
	}

	if (threadIdx.y >= 4 && threadIdx.y < 8)	// 4, 5, 6, 7
	{
		sm_b[0 + (threadIdx.y - 4) * 4][threadIdx.x] 		= reg_tile[0][0];
        sm_b[1 + (threadIdx.y - 4) * 4][threadIdx.x] 		= reg_tile[1][0];
        sm_b[2 + (threadIdx.y - 4) * 4][threadIdx.x] 		= reg_tile[2][0];
        sm_b[3 + (threadIdx.y - 4) * 4][threadIdx.x] 		= reg_tile[3][0];

        sm_b[0 + (threadIdx.y - 4) * 4][threadIdx.x + 16] 	= reg_tile[0][1];
        sm_b[1 + (threadIdx.y - 4) * 4][threadIdx.x + 16] 	= reg_tile[1][1];
        sm_b[2 + (threadIdx.y - 4) * 4][threadIdx.x + 16] 	= reg_tile[2][1];
        sm_b[3 + (threadIdx.y - 4) * 4][threadIdx.x + 16] 	= reg_tile[3][1];

        sm_b[0 + (threadIdx.y - 4) * 4][threadIdx.x + 32] 	= reg_tile[0][2];
        sm_b[1 + (threadIdx.y - 4) * 4][threadIdx.x + 32] 	= reg_tile[1][2];
        sm_b[2 + (threadIdx.y - 4) * 4][threadIdx.x + 32] 	= reg_tile[2][2];
        sm_b[3 + (threadIdx.y - 4) * 4][threadIdx.x + 32] 	= reg_tile[3][2];

        sm_b[0 + (threadIdx.y - 4) * 4][threadIdx.x + 48] 	= reg_tile[0][3];
        sm_b[1 + (threadIdx.y - 4) * 4][threadIdx.x + 48] 	= reg_tile[1][3];
        sm_b[2 + (threadIdx.y - 4) * 4][threadIdx.x + 48] 	= reg_tile[2][3];
        sm_b[3 + (threadIdx.y - 4) * 4][threadIdx.x + 48] 	= reg_tile[3][3];
	}
	__syncthreads();

	if (threadIdx.y < 4)						// 0, 1, 2, 3
	{
		reg_tile[0][0] = sm_a[threadIdx.y + 0][(threadIdx.x)];
		reg_tile[1][0] = sm_a[threadIdx.y + 4][(threadIdx.x)];
		reg_tile[2][0] = sm_a[threadIdx.y + 8][(threadIdx.x)];
		reg_tile[3][0] = sm_a[threadIdx.y + 12][(threadIdx.x)];

		reg_tile[0][1] = sm_a[threadIdx.y + 0][(threadIdx.x) + 16];
		reg_tile[1][1] = sm_a[threadIdx.y + 4][(threadIdx.x) + 16];
		reg_tile[2][1] = sm_a[threadIdx.y + 8][(threadIdx.x) + 16];
        reg_tile[3][1] = sm_a[threadIdx.y + 12][(threadIdx.x) + 16];
        
		reg_tile[0][2] = sm_a[threadIdx.y + 0][(threadIdx.x) + 32];
		reg_tile[1][2] = sm_a[threadIdx.y + 4][(threadIdx.x) + 32];
		reg_tile[2][2] = sm_a[threadIdx.y + 8][(threadIdx.x) + 32];
		reg_tile[3][2] = sm_a[threadIdx.y + 12][(threadIdx.x) + 32];
        
        reg_tile[0][3] = sm_a[threadIdx.y + 0][(threadIdx.x) + 48];
        reg_tile[1][3] = sm_a[threadIdx.y + 4][(threadIdx.x) + 48];
        reg_tile[2][3] = sm_a[threadIdx.y + 8][(threadIdx.x) + 48];
        reg_tile[3][3] = sm_a[threadIdx.y + 12][(threadIdx.x) + 48];
	}

	if (threadIdx.y >= 4 && threadIdx.y < 8)	// 4, 5, 6, 7
	{
		reg_tile[0][0] = sm_b[(threadIdx.y - 4) + 0][(threadIdx.x)];
		reg_tile[1][0] = sm_b[(threadIdx.y - 4) + 4][(threadIdx.x)];
		reg_tile[2][0] = sm_b[(threadIdx.y - 4) + 8][(threadIdx.x)];
		reg_tile[3][0] = sm_b[(threadIdx.y - 4) + 12][(threadIdx.x)];

		reg_tile[0][1] = sm_b[(threadIdx.y - 4) + 0][(threadIdx.x) + 16];
        reg_tile[1][1] = sm_b[(threadIdx.y - 4) + 4][(threadIdx.x) + 16];
		reg_tile[2][1] = sm_b[(threadIdx.y - 4) + 8][(threadIdx.x) + 16];
		reg_tile[3][1] = sm_b[(threadIdx.y - 4) + 12][(threadIdx.x) + 16];
		
		reg_tile[0][2] = sm_b[(threadIdx.y - 4) + 0][(threadIdx.x) + 32];
		reg_tile[1][2] = sm_b[(threadIdx.y - 4) + 4][(threadIdx.x) + 32];
		reg_tile[2][2] = sm_b[(threadIdx.y - 4) + 8][(threadIdx.x) + 32];
		reg_tile[3][2] = sm_b[(threadIdx.y - 4) + 12][(threadIdx.x) + 32];
        
		reg_tile[0][3] = sm_b[(threadIdx.y - 4) + 0][(threadIdx.x) + 48];	
        reg_tile[1][3] = sm_b[(threadIdx.y - 4) + 4][(threadIdx.x) + 48];
        reg_tile[2][3] = sm_b[(threadIdx.y - 4) + 8][(threadIdx.x) + 48];
        reg_tile[3][3] = sm_b[(threadIdx.y - 4) + 12][(threadIdx.x) + 48];
	}
	__syncthreads();

	if (threadIdx.y >= 8 && threadIdx.y < 12)	// 8, 9, 10, 11
	{
		sm_a[0 + (threadIdx.y - 8) * 4][threadIdx.x] = reg_tile[0][0];
        sm_a[1 + (threadIdx.y - 8) * 4][threadIdx.x] = reg_tile[1][0];
        sm_a[2 + (threadIdx.y - 8) * 4][threadIdx.x] = reg_tile[2][0];
        sm_a[3 + (threadIdx.y - 8) * 4][threadIdx.x] = reg_tile[3][0];

        sm_a[0 + (threadIdx.y - 8) * 4][threadIdx.x + 16] = reg_tile[0][1];
        sm_a[1 + (threadIdx.y - 8) * 4][threadIdx.x + 16] = reg_tile[1][1];
        sm_a[2 + (threadIdx.y - 8) * 4][threadIdx.x + 16] = reg_tile[2][1];
        sm_a[3 + (threadIdx.y - 8) * 4][threadIdx.x + 16] = reg_tile[3][1];

        sm_a[0 + (threadIdx.y - 8) * 4][threadIdx.x + 32] = reg_tile[0][2];
        sm_a[1 + (threadIdx.y - 8) * 4][threadIdx.x + 32] = reg_tile[1][2];
        sm_a[2 + (threadIdx.y - 8) * 4][threadIdx.x + 32] = reg_tile[2][2];
        sm_a[3 + (threadIdx.y - 8) * 4][threadIdx.x + 32] = reg_tile[3][2];

        sm_a[0 + (threadIdx.y - 8) * 4][threadIdx.x + 48] = reg_tile[0][3];
        sm_a[1 + (threadIdx.y - 8) * 4][threadIdx.x + 48] = reg_tile[1][3];
        sm_a[2 + (threadIdx.y - 8) * 4][threadIdx.x + 48] = reg_tile[2][3];
        sm_a[3 + (threadIdx.y - 8) * 4][threadIdx.x + 48] = reg_tile[3][3];
	}

	if (threadIdx.y >= 12)	// 12, 13, 14, 15
	{
		sm_b[0 + (threadIdx.y - 12) * 4][threadIdx.x] = reg_tile[0][0];
        sm_b[1 + (threadIdx.y - 12) * 4][threadIdx.x] = reg_tile[1][0];
        sm_b[2 + (threadIdx.y - 12) * 4][threadIdx.x] = reg_tile[2][0];
		sm_b[3 + (threadIdx.y - 12) * 4][threadIdx.x] = reg_tile[3][0];
		
        sm_b[0 + (threadIdx.y - 12) * 4][threadIdx.x + 16] = reg_tile[0][1];
        sm_b[1 + (threadIdx.y - 12) * 4][threadIdx.x + 16] = reg_tile[1][1];
        sm_b[2 + (threadIdx.y - 12) * 4][threadIdx.x + 16] = reg_tile[2][1];
		sm_b[3 + (threadIdx.y - 12) * 4][threadIdx.x + 16] = reg_tile[3][1];
		
        sm_b[0 + (threadIdx.y - 12) * 4][threadIdx.x + 32] = reg_tile[0][2];
        sm_b[1 + (threadIdx.y - 12) * 4][threadIdx.x + 32] = reg_tile[1][2];
        sm_b[2 + (threadIdx.y - 12) * 4][threadIdx.x + 32] = reg_tile[2][2];
		sm_b[3 + (threadIdx.y - 12) * 4][threadIdx.x + 32] = reg_tile[3][2];
		
        sm_b[0 + (threadIdx.y - 12) * 4][threadIdx.x + 48] = reg_tile[0][3];
        sm_b[1 + (threadIdx.y - 12) * 4][threadIdx.x + 48] = reg_tile[1][3];
        sm_b[2 + (threadIdx.y - 12) * 4][threadIdx.x + 48] = reg_tile[2][3];
        sm_b[3 + (threadIdx.y - 12) * 4][threadIdx.x + 48] = reg_tile[3][3];
	}
	__syncthreads();

	if (threadIdx.y >= 8 && threadIdx.y < 12)	// 8, 9, 10, 11
	{
		reg_tile[0][0] = sm_a[(threadIdx.y - 8) + 0][(threadIdx.x)];
		reg_tile[1][0] = sm_a[(threadIdx.y - 8) + 4][(threadIdx.x)];
		reg_tile[2][0] = sm_a[(threadIdx.y - 8) + 8][(threadIdx.x)];
		reg_tile[3][0] = sm_a[(threadIdx.y - 8) + 12][(threadIdx.x)];

		reg_tile[0][1] = sm_a[(threadIdx.y - 8) + 0][(threadIdx.x) + 16];
		reg_tile[1][1] = sm_a[(threadIdx.y - 8) + 4][(threadIdx.x) + 16];
		reg_tile[2][1] = sm_a[(threadIdx.y - 8) + 8][(threadIdx.x) + 16];
		reg_tile[3][1] = sm_a[(threadIdx.y - 8) + 12][(threadIdx.x) + 16];

		reg_tile[0][2] = sm_a[(threadIdx.y - 8) + 0][(threadIdx.x) + 32];
		reg_tile[1][2] = sm_a[(threadIdx.y - 8) + 4][(threadIdx.x) + 32];		
        reg_tile[2][2] = sm_a[(threadIdx.y - 8) + 8][(threadIdx.x) + 32];
        reg_tile[3][2] = sm_a[(threadIdx.y - 8) + 12][(threadIdx.x) + 32];

        reg_tile[0][3] = sm_a[(threadIdx.y - 8) + 0][(threadIdx.x) + 48];
        reg_tile[1][3] = sm_a[(threadIdx.y - 8) + 4][(threadIdx.x) + 48];
        reg_tile[2][3] = sm_a[(threadIdx.y - 8) + 8][(threadIdx.x) + 48];
		reg_tile[3][3] = sm_a[(threadIdx.y - 8) + 12][(threadIdx.x) + 48];
	}

	if (threadIdx.y >= 12)	// 12, 13, 14, 15
	{
		reg_tile[0][0] = sm_b[(threadIdx.y - 12) + 0][(threadIdx.x)];
		reg_tile[1][0] = sm_b[(threadIdx.y - 12) + 4][(threadIdx.x)];
		reg_tile[2][0] = sm_b[(threadIdx.y - 12) + 8][(threadIdx.x)];
		reg_tile[3][0] = sm_b[(threadIdx.y - 12) + 12][(threadIdx.x)];

		reg_tile[0][1] = sm_b[(threadIdx.y - 12) + 0][(threadIdx.x) + 16];
		reg_tile[1][1] = sm_b[(threadIdx.y - 12) + 4][(threadIdx.x) + 16];
		reg_tile[2][1] = sm_b[(threadIdx.y - 12) + 8][(threadIdx.x) + 16];
		reg_tile[3][1] = sm_b[(threadIdx.y - 12) + 12][(threadIdx.x) + 16];

		reg_tile[0][2] = sm_b[(threadIdx.y - 12) + 0][(threadIdx.x) + 32];
		reg_tile[1][2] = sm_b[(threadIdx.y - 12) + 4][(threadIdx.x) + 32];
		reg_tile[2][2] = sm_b[(threadIdx.y - 12) + 8][(threadIdx.x) + 32];
        reg_tile[3][2] = sm_b[(threadIdx.y - 12) + 12][(threadIdx.x) + 32];

		reg_tile[0][3] = sm_b[(threadIdx.y - 12) + 0][(threadIdx.x) + 48];
		reg_tile[1][3] = sm_b[(threadIdx.y - 12) + 4][(threadIdx.x) + 48];
		reg_tile[2][3] = sm_b[(threadIdx.y - 12) + 8][(threadIdx.x) + 48];
        reg_tile[3][3] = sm_b[(threadIdx.y - 12) + 12][(threadIdx.x) + 48];
	}
	__syncthreads();

	//
	//	sd2_1, 2, 3, 4, 5 and 6.
	//
	// tensor contraction
	internal_upperbound = 0;
	
	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_1 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h1 && idx_h1 < rng_h2 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_1[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h1) * size_h1) * size_p4) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h3 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_1[(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{   
            temp_bv[0] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 48];
            
			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H3 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_2 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h2 && idx_h1 < rng_h3 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_2[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h2) * size_p4) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h1 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_2[(blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 0];
            temp_bv[1] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 16];
            temp_bv[2] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 32];
            temp_bv[3] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_3 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h1 && idx_h1 < rng_h3 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_3[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h1) * size_p4) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h2 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_3[(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_4 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h1 && idx_h1 < rng_h2 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_4[(blk_idx_p5 * SIZE_SLICE_1_P5 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h1) * size_h1) * size_p5) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h3 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_4[(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_5 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h2 && idx_h1 < rng_h3 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_5[(blk_idx_p5 * SIZE_SLICE_1_P5 + ll + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h2) * size_p5) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h1 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_5[(blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_6 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h1 && idx_h1 < rng_h3 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_6[(blk_idx_p5 * SIZE_SLICE_1_P6 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h1) * size_p5) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h2 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_6[(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 0];
            temp_bv[1] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 16];
            temp_bv[2] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 32];
            temp_bv[3] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 48];

			for (int xx = 0; xx < 4; xx++)	// 4 -> rng_p4: Local Transactions...
			{
                temp_av = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_h3 < rng_h3 && idx_h2 < rng_h2 && idx_p6 < rng_p6 && idx_h1 < rng_h1)
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_p4 && j < rng_p5)
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void kernel_ccsdT_sd2_fully_fused_partial_full(double* t3, double* d_t2_1, double* d_t2_2, double* d_t2_3, double* d_t2_4, double* d_t2_5, double* d_t2_6, double* d_v2_1, double* d_v2_2, double* d_v2_3, double* d_v2_4, double* d_v2_5, double* d_v2_6, double* d_t2_7, double* d_t2_8, double* d_t2_9, double* d_v2_7, double* d_v2_8, double* d_v2_9, int size_h3, int size_h2, int size_h1, int size_p6, int size_p5, int size_p4, int size_p7, int numBlk_h3, int numBlk_h2, int numBlk_h1, int numBlk_p6, int numBlk_p5, int numBlk_p4, int kernel_1, int kernel_2, int kernel_3, int kernel_4, int kernel_5, int kernel_6, int kernel_7, int kernel_8, int kernel_9, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int internal_upperbound   	= 0;
	int internal_offset;

    // should support for non-full tiles
	int idx_h3 = threadIdx.x % SIZE_SLICE_1_H3;
	int idx_h2 = threadIdx.x / SIZE_SLICE_1_H3;
	int idx_p6 = threadIdx.y % SIZE_SLICE_1_P6;
	int idx_h1 = threadIdx.y / SIZE_SLICE_1_P6;

    // Common for Threads within a Thread Block
	int tmp_blkIdx;        
    int blk_idx_p4  = blockIdx.x / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);
    tmp_blkIdx      = blockIdx.x % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);

    int blk_idx_p5  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);

    int blk_idx_p6  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1);

    int blk_idx_h1 = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2);
    tmp_blkIdx     = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2);

    int blk_idx_h2 = (tmp_blkIdx) / (numBlk_h3);
    int blk_idx_h3 = (tmp_blkIdx) % (numBlk_h3);

    int t3_base_thread = blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h3 + 
    (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h2 + 
    (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_h1 + 
    (blk_idx_p6 * SIZE_SLICE_2_P6 + idx_p6 + 
    (blk_idx_p5 * SIZE_SLICE_2_P5 + 
    (blk_idx_p4 * SIZE_SLICE_2_P4) * size_p5) * size_p6) * size_h1) * size_h2) * size_h3;


	double temp_av;
	double temp_bv[4];
	double reg_tile[4][4];

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	//
	//	sd2_7, 8 and 9
	//
	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_7 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_7[(blk_idx_p6 *  SIZE_SLICE_2_P6 + ll + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h1) * size_h1) * size_p6) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_7[(blk_idx_h3 * SIZE_SLICE_2_H3 + idx_p6 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_h1) * size_p5) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_2_H3 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_8 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_8[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h1) * size_h2) * size_p6) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_8[(blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p6 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_1_P4 + idx_h1) * size_p5) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 0];
            temp_bv[1] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 16];
            temp_bv[2] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 32];
            temp_bv[3] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_2_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_9 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_9[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h1) * size_h1) * size_p6) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_9[(blk_idx_h2 * SIZE_SLICE_2_H2 + idx_p6 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_h1) * size_p5) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_2_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	//
	//
	//
	if (threadIdx.y < 4)						// 0, 1, 2, 3
	{
		// sm_a[16][64] <-- (4 x 16) x (4 x 4) = (16 x 64)		  'y''x'
		sm_a[0 + threadIdx.y * 4][threadIdx.x] 			= reg_tile[0][0];
		sm_a[1 + threadIdx.y * 4][threadIdx.x] 			= reg_tile[1][0];
		sm_a[2 + threadIdx.y * 4][threadIdx.x] 			= reg_tile[2][0];
		sm_a[3 + threadIdx.y * 4][threadIdx.x] 			= reg_tile[3][0];
		
		sm_a[0 + threadIdx.y * 4][threadIdx.x + 16] 	= reg_tile[0][1];
		sm_a[1 + threadIdx.y * 4][threadIdx.x + 16] 	= reg_tile[1][1];
		sm_a[2 + threadIdx.y * 4][threadIdx.x + 16] 	= reg_tile[2][1];
		sm_a[3 + threadIdx.y * 4][threadIdx.x + 16] 	= reg_tile[3][1];
		         
		sm_a[0 + threadIdx.y * 4][threadIdx.x + 32] 	= reg_tile[0][2];
		sm_a[1 + threadIdx.y * 4][threadIdx.x + 32] 	= reg_tile[1][2];
		sm_a[2 + threadIdx.y * 4][threadIdx.x + 32] 	= reg_tile[2][2];
		sm_a[3 + threadIdx.y * 4][threadIdx.x + 32] 	= reg_tile[3][2];
		
		sm_a[0 + threadIdx.y * 4][threadIdx.x + 48] 	= reg_tile[0][3];
		sm_a[1 + threadIdx.y * 4][threadIdx.x + 48] 	= reg_tile[1][3];
		sm_a[2 + threadIdx.y * 4][threadIdx.x + 48] 	= reg_tile[2][3];
		sm_a[3 + threadIdx.y * 4][threadIdx.x + 48] 	= reg_tile[3][3];
	}

	if (threadIdx.y >= 4 && threadIdx.y < 8)	// 4, 5, 6, 7
	{
		sm_b[0 + (threadIdx.y - 4) * 4][threadIdx.x] 		= reg_tile[0][0];
        sm_b[1 + (threadIdx.y - 4) * 4][threadIdx.x] 		= reg_tile[1][0];
        sm_b[2 + (threadIdx.y - 4) * 4][threadIdx.x] 		= reg_tile[2][0];
        sm_b[3 + (threadIdx.y - 4) * 4][threadIdx.x] 		= reg_tile[3][0];

        sm_b[0 + (threadIdx.y - 4) * 4][threadIdx.x + 16] 	= reg_tile[0][1];
        sm_b[1 + (threadIdx.y - 4) * 4][threadIdx.x + 16] 	= reg_tile[1][1];
        sm_b[2 + (threadIdx.y - 4) * 4][threadIdx.x + 16] 	= reg_tile[2][1];
        sm_b[3 + (threadIdx.y - 4) * 4][threadIdx.x + 16] 	= reg_tile[3][1];

        sm_b[0 + (threadIdx.y - 4) * 4][threadIdx.x + 32] 	= reg_tile[0][2];
        sm_b[1 + (threadIdx.y - 4) * 4][threadIdx.x + 32] 	= reg_tile[1][2];
        sm_b[2 + (threadIdx.y - 4) * 4][threadIdx.x + 32] 	= reg_tile[2][2];
        sm_b[3 + (threadIdx.y - 4) * 4][threadIdx.x + 32] 	= reg_tile[3][2];

        sm_b[0 + (threadIdx.y - 4) * 4][threadIdx.x + 48] 	= reg_tile[0][3];
        sm_b[1 + (threadIdx.y - 4) * 4][threadIdx.x + 48] 	= reg_tile[1][3];
        sm_b[2 + (threadIdx.y - 4) * 4][threadIdx.x + 48] 	= reg_tile[2][3];
        sm_b[3 + (threadIdx.y - 4) * 4][threadIdx.x + 48] 	= reg_tile[3][3];
	}
	__syncthreads();

	if (threadIdx.y < 4)						// 0, 1, 2, 3
	{
		reg_tile[0][0] = sm_a[threadIdx.y + 0][(threadIdx.x)];
		reg_tile[1][0] = sm_a[threadIdx.y + 4][(threadIdx.x)];
		reg_tile[2][0] = sm_a[threadIdx.y + 8][(threadIdx.x)];
		reg_tile[3][0] = sm_a[threadIdx.y + 12][(threadIdx.x)];

		reg_tile[0][1] = sm_a[threadIdx.y + 0][(threadIdx.x) + 16];
		reg_tile[1][1] = sm_a[threadIdx.y + 4][(threadIdx.x) + 16];
		reg_tile[2][1] = sm_a[threadIdx.y + 8][(threadIdx.x) + 16];
        reg_tile[3][1] = sm_a[threadIdx.y + 12][(threadIdx.x) + 16];
        
		reg_tile[0][2] = sm_a[threadIdx.y + 0][(threadIdx.x) + 32];
		reg_tile[1][2] = sm_a[threadIdx.y + 4][(threadIdx.x) + 32];
		reg_tile[2][2] = sm_a[threadIdx.y + 8][(threadIdx.x) + 32];
		reg_tile[3][2] = sm_a[threadIdx.y + 12][(threadIdx.x) + 32];
        
        reg_tile[0][3] = sm_a[threadIdx.y + 0][(threadIdx.x) + 48];
        reg_tile[1][3] = sm_a[threadIdx.y + 4][(threadIdx.x) + 48];
        reg_tile[2][3] = sm_a[threadIdx.y + 8][(threadIdx.x) + 48];
        reg_tile[3][3] = sm_a[threadIdx.y + 12][(threadIdx.x) + 48];
	}

	if (threadIdx.y >= 4 && threadIdx.y < 8)	// 4, 5, 6, 7
	{
		reg_tile[0][0] = sm_b[(threadIdx.y - 4) + 0][(threadIdx.x)];
		reg_tile[1][0] = sm_b[(threadIdx.y - 4) + 4][(threadIdx.x)];
		reg_tile[2][0] = sm_b[(threadIdx.y - 4) + 8][(threadIdx.x)];
		reg_tile[3][0] = sm_b[(threadIdx.y - 4) + 12][(threadIdx.x)];

		reg_tile[0][1] = sm_b[(threadIdx.y - 4) + 0][(threadIdx.x) + 16];
        reg_tile[1][1] = sm_b[(threadIdx.y - 4) + 4][(threadIdx.x) + 16];
		reg_tile[2][1] = sm_b[(threadIdx.y - 4) + 8][(threadIdx.x) + 16];
		reg_tile[3][1] = sm_b[(threadIdx.y - 4) + 12][(threadIdx.x) + 16];
		
		reg_tile[0][2] = sm_b[(threadIdx.y - 4) + 0][(threadIdx.x) + 32];
		reg_tile[1][2] = sm_b[(threadIdx.y - 4) + 4][(threadIdx.x) + 32];
		reg_tile[2][2] = sm_b[(threadIdx.y - 4) + 8][(threadIdx.x) + 32];
		reg_tile[3][2] = sm_b[(threadIdx.y - 4) + 12][(threadIdx.x) + 32];
        
		reg_tile[0][3] = sm_b[(threadIdx.y - 4) + 0][(threadIdx.x) + 48];	
        reg_tile[1][3] = sm_b[(threadIdx.y - 4) + 4][(threadIdx.x) + 48];
        reg_tile[2][3] = sm_b[(threadIdx.y - 4) + 8][(threadIdx.x) + 48];
        reg_tile[3][3] = sm_b[(threadIdx.y - 4) + 12][(threadIdx.x) + 48];
	}
	__syncthreads();

	if (threadIdx.y >= 8 && threadIdx.y < 12)	// 8, 9, 10, 11
	{
		sm_a[0 + (threadIdx.y - 8) * 4][threadIdx.x] = reg_tile[0][0];
        sm_a[1 + (threadIdx.y - 8) * 4][threadIdx.x] = reg_tile[1][0];
        sm_a[2 + (threadIdx.y - 8) * 4][threadIdx.x] = reg_tile[2][0];
        sm_a[3 + (threadIdx.y - 8) * 4][threadIdx.x] = reg_tile[3][0];

        sm_a[0 + (threadIdx.y - 8) * 4][threadIdx.x + 16] = reg_tile[0][1];
        sm_a[1 + (threadIdx.y - 8) * 4][threadIdx.x + 16] = reg_tile[1][1];
        sm_a[2 + (threadIdx.y - 8) * 4][threadIdx.x + 16] = reg_tile[2][1];
        sm_a[3 + (threadIdx.y - 8) * 4][threadIdx.x + 16] = reg_tile[3][1];

        sm_a[0 + (threadIdx.y - 8) * 4][threadIdx.x + 32] = reg_tile[0][2];
        sm_a[1 + (threadIdx.y - 8) * 4][threadIdx.x + 32] = reg_tile[1][2];
        sm_a[2 + (threadIdx.y - 8) * 4][threadIdx.x + 32] = reg_tile[2][2];
        sm_a[3 + (threadIdx.y - 8) * 4][threadIdx.x + 32] = reg_tile[3][2];

        sm_a[0 + (threadIdx.y - 8) * 4][threadIdx.x + 48] = reg_tile[0][3];
        sm_a[1 + (threadIdx.y - 8) * 4][threadIdx.x + 48] = reg_tile[1][3];
        sm_a[2 + (threadIdx.y - 8) * 4][threadIdx.x + 48] = reg_tile[2][3];
        sm_a[3 + (threadIdx.y - 8) * 4][threadIdx.x + 48] = reg_tile[3][3];
	}

	if (threadIdx.y >= 12)	// 12, 13, 14, 15
	{
		sm_b[0 + (threadIdx.y - 12) * 4][threadIdx.x] = reg_tile[0][0];
        sm_b[1 + (threadIdx.y - 12) * 4][threadIdx.x] = reg_tile[1][0];
        sm_b[2 + (threadIdx.y - 12) * 4][threadIdx.x] = reg_tile[2][0];
		sm_b[3 + (threadIdx.y - 12) * 4][threadIdx.x] = reg_tile[3][0];
		
        sm_b[0 + (threadIdx.y - 12) * 4][threadIdx.x + 16] = reg_tile[0][1];
        sm_b[1 + (threadIdx.y - 12) * 4][threadIdx.x + 16] = reg_tile[1][1];
        sm_b[2 + (threadIdx.y - 12) * 4][threadIdx.x + 16] = reg_tile[2][1];
		sm_b[3 + (threadIdx.y - 12) * 4][threadIdx.x + 16] = reg_tile[3][1];
		
        sm_b[0 + (threadIdx.y - 12) * 4][threadIdx.x + 32] = reg_tile[0][2];
        sm_b[1 + (threadIdx.y - 12) * 4][threadIdx.x + 32] = reg_tile[1][2];
        sm_b[2 + (threadIdx.y - 12) * 4][threadIdx.x + 32] = reg_tile[2][2];
		sm_b[3 + (threadIdx.y - 12) * 4][threadIdx.x + 32] = reg_tile[3][2];
		
        sm_b[0 + (threadIdx.y - 12) * 4][threadIdx.x + 48] = reg_tile[0][3];
        sm_b[1 + (threadIdx.y - 12) * 4][threadIdx.x + 48] = reg_tile[1][3];
        sm_b[2 + (threadIdx.y - 12) * 4][threadIdx.x + 48] = reg_tile[2][3];
        sm_b[3 + (threadIdx.y - 12) * 4][threadIdx.x + 48] = reg_tile[3][3];
	}
	__syncthreads();

	if (threadIdx.y >= 8 && threadIdx.y < 12)	// 8, 9, 10, 11
	{
		reg_tile[0][0] = sm_a[(threadIdx.y - 8) + 0][(threadIdx.x)];
		reg_tile[1][0] = sm_a[(threadIdx.y - 8) + 4][(threadIdx.x)];
		reg_tile[2][0] = sm_a[(threadIdx.y - 8) + 8][(threadIdx.x)];
		reg_tile[3][0] = sm_a[(threadIdx.y - 8) + 12][(threadIdx.x)];

		reg_tile[0][1] = sm_a[(threadIdx.y - 8) + 0][(threadIdx.x) + 16];
		reg_tile[1][1] = sm_a[(threadIdx.y - 8) + 4][(threadIdx.x) + 16];
		reg_tile[2][1] = sm_a[(threadIdx.y - 8) + 8][(threadIdx.x) + 16];
		reg_tile[3][1] = sm_a[(threadIdx.y - 8) + 12][(threadIdx.x) + 16];

		reg_tile[0][2] = sm_a[(threadIdx.y - 8) + 0][(threadIdx.x) + 32];
		reg_tile[1][2] = sm_a[(threadIdx.y - 8) + 4][(threadIdx.x) + 32];		
        reg_tile[2][2] = sm_a[(threadIdx.y - 8) + 8][(threadIdx.x) + 32];
        reg_tile[3][2] = sm_a[(threadIdx.y - 8) + 12][(threadIdx.x) + 32];

        reg_tile[0][3] = sm_a[(threadIdx.y - 8) + 0][(threadIdx.x) + 48];
        reg_tile[1][3] = sm_a[(threadIdx.y - 8) + 4][(threadIdx.x) + 48];
        reg_tile[2][3] = sm_a[(threadIdx.y - 8) + 8][(threadIdx.x) + 48];
		reg_tile[3][3] = sm_a[(threadIdx.y - 8) + 12][(threadIdx.x) + 48];
	}

	if (threadIdx.y >= 12)	// 12, 13, 14, 15
	{
		reg_tile[0][0] = sm_b[(threadIdx.y - 12) + 0][(threadIdx.x)];
		reg_tile[1][0] = sm_b[(threadIdx.y - 12) + 4][(threadIdx.x)];
		reg_tile[2][0] = sm_b[(threadIdx.y - 12) + 8][(threadIdx.x)];
		reg_tile[3][0] = sm_b[(threadIdx.y - 12) + 12][(threadIdx.x)];

		reg_tile[0][1] = sm_b[(threadIdx.y - 12) + 0][(threadIdx.x) + 16];
		reg_tile[1][1] = sm_b[(threadIdx.y - 12) + 4][(threadIdx.x) + 16];
		reg_tile[2][1] = sm_b[(threadIdx.y - 12) + 8][(threadIdx.x) + 16];
		reg_tile[3][1] = sm_b[(threadIdx.y - 12) + 12][(threadIdx.x) + 16];

		reg_tile[0][2] = sm_b[(threadIdx.y - 12) + 0][(threadIdx.x) + 32];
		reg_tile[1][2] = sm_b[(threadIdx.y - 12) + 4][(threadIdx.x) + 32];
		reg_tile[2][2] = sm_b[(threadIdx.y - 12) + 8][(threadIdx.x) + 32];
        reg_tile[3][2] = sm_b[(threadIdx.y - 12) + 12][(threadIdx.x) + 32];

		reg_tile[0][3] = sm_b[(threadIdx.y - 12) + 0][(threadIdx.x) + 48];
		reg_tile[1][3] = sm_b[(threadIdx.y - 12) + 4][(threadIdx.x) + 48];
		reg_tile[2][3] = sm_b[(threadIdx.y - 12) + 8][(threadIdx.x) + 48];
        reg_tile[3][3] = sm_b[(threadIdx.y - 12) + 12][(threadIdx.x) + 48];
	}
	__syncthreads();

	//
	//	sd2_1, 2, 3, 4, 5 and 6.
	//
	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_1 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_1[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h1) * size_h1) * size_p4) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_1[(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 48];
            
			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H3 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_2 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_2[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h2) * size_p4) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_2[(blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 0];
            temp_bv[1] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 16];
            temp_bv[2] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 32];
            temp_bv[3] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_3 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_3[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h1) * size_p4) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_3[(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_4 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_4[(blk_idx_p5 * SIZE_SLICE_1_P5 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h1) * size_h1) * size_p5) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_4[(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_5 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_5[(blk_idx_p5 * SIZE_SLICE_1_P5 + ll + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h2) * size_p5) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_5[(blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_6 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_6[(blk_idx_p5 * SIZE_SLICE_1_P6 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h1) * size_p5) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_6[(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 0];
            temp_bv[1] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 16];
            temp_bv[2] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 32];
            temp_bv[3] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 48];

			for (int xx = 0; xx < 4; xx++)	// 4 -> rng_p4: Local Transactions...
			{
                temp_av = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// Store Results (Registers) to Global Memory
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void kernel_ccsdT_sd2_fully_fused_full_full(double* t3, double* d_t2_1, double* d_t2_2, double* d_t2_3, double* d_t2_4, double* d_t2_5, double* d_t2_6, double* d_v2_1, double* d_v2_2, double* d_v2_3, double* d_v2_4, double* d_v2_5, double* d_v2_6, double* d_t2_7, double* d_t2_8, double* d_t2_9, double* d_v2_7, double* d_v2_8, double* d_v2_9, int size_h3,    int size_h2,    int size_h1,    int size_p6,    int size_p5,    int size_p4,    int size_p7, int numBlk_h3,  int numBlk_h2,  int numBlk_h1,  int numBlk_p6,  int numBlk_p5,  int numBlk_p4, int kernel_1, int kernel_2, int kernel_3, int kernel_4, int kernel_5, int kernel_6, int kernel_7, int kernel_8, int kernel_9, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

    // should support for non-full tiles
	int idx_h3 = threadIdx.x % SIZE_SLICE_1_H3;
	int idx_h2 = threadIdx.x / SIZE_SLICE_1_H3;
	int idx_p6 = threadIdx.y % SIZE_SLICE_1_P6;
	int idx_h1 = threadIdx.y / SIZE_SLICE_1_P6;

    // Common for Threads within a Thread Block
	int tmp_blkIdx;        
    int blk_idx_p4  = blockIdx.x / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);
    tmp_blkIdx      = blockIdx.x % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);

    int blk_idx_p5  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);

    int blk_idx_p6  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1);

    int blk_idx_h1 = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2);
    tmp_blkIdx     = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2);

    int blk_idx_h2 = (tmp_blkIdx) / (numBlk_h3);
    int blk_idx_h3 = (tmp_blkIdx) % (numBlk_h3);

    int t3_base_thread = blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h3 + 
    (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h2 + 
    (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_h1 + 
    (blk_idx_p6 * SIZE_SLICE_2_P6 + idx_p6 + 
    (blk_idx_p5 * SIZE_SLICE_2_P5 + 
    (blk_idx_p4 * SIZE_SLICE_2_P4) * size_p5) * size_p6) * size_h1) * size_h2) * size_h3;

	double temp_av;
	double temp_bv[4];
	double reg_tile[4][4];

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	//
	//	sd2_7, 8 and 9
	//
	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_7 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
            sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_7[(blk_idx_p6 *  SIZE_SLICE_2_P6 + ll + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h1) * size_h1) * size_p6) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_7[(blk_idx_h3 * SIZE_SLICE_2_H3 + idx_p6 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_h1) * size_p5) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_2_H3 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_8 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
            sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_8[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h1) * size_h2) * size_p6) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_8[(blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p6 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_1_P4 + idx_h1) * size_p5) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 0];
            temp_bv[1] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 16];
            temp_bv[2] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 32];
            temp_bv[3] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_2_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_9 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
            sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_9[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h1) * size_h1) * size_p6) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_9[(blk_idx_h2 * SIZE_SLICE_2_H2 + idx_p6 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_h1) * size_p5) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_2_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	//
	//
	//
	if (threadIdx.y < 4)						// 0, 1, 2, 3
	{
		// sm_a[16][64] <-- (4 x 16) x (4 x 4) = (16 x 64)		  'y''x'
		sm_a[0 + threadIdx.y * 4][threadIdx.x] 			= reg_tile[0][0];
		sm_a[1 + threadIdx.y * 4][threadIdx.x] 			= reg_tile[1][0];
		sm_a[2 + threadIdx.y * 4][threadIdx.x] 			= reg_tile[2][0];
		sm_a[3 + threadIdx.y * 4][threadIdx.x] 			= reg_tile[3][0];
		
		sm_a[0 + threadIdx.y * 4][threadIdx.x + 16] 	= reg_tile[0][1];
		sm_a[1 + threadIdx.y * 4][threadIdx.x + 16] 	= reg_tile[1][1];
		sm_a[2 + threadIdx.y * 4][threadIdx.x + 16] 	= reg_tile[2][1];
		sm_a[3 + threadIdx.y * 4][threadIdx.x + 16] 	= reg_tile[3][1];
		         
		sm_a[0 + threadIdx.y * 4][threadIdx.x + 32] 	= reg_tile[0][2];
		sm_a[1 + threadIdx.y * 4][threadIdx.x + 32] 	= reg_tile[1][2];
		sm_a[2 + threadIdx.y * 4][threadIdx.x + 32] 	= reg_tile[2][2];
		sm_a[3 + threadIdx.y * 4][threadIdx.x + 32] 	= reg_tile[3][2];
		
		sm_a[0 + threadIdx.y * 4][threadIdx.x + 48] 	= reg_tile[0][3];
		sm_a[1 + threadIdx.y * 4][threadIdx.x + 48] 	= reg_tile[1][3];
		sm_a[2 + threadIdx.y * 4][threadIdx.x + 48] 	= reg_tile[2][3];
		sm_a[3 + threadIdx.y * 4][threadIdx.x + 48] 	= reg_tile[3][3];
	}

	if (threadIdx.y >= 4 && threadIdx.y < 8)	// 4, 5, 6, 7
	{
		sm_b[0 + (threadIdx.y - 4) * 4][threadIdx.x] 		= reg_tile[0][0];
        sm_b[1 + (threadIdx.y - 4) * 4][threadIdx.x] 		= reg_tile[1][0];
        sm_b[2 + (threadIdx.y - 4) * 4][threadIdx.x] 		= reg_tile[2][0];
        sm_b[3 + (threadIdx.y - 4) * 4][threadIdx.x] 		= reg_tile[3][0];

        sm_b[0 + (threadIdx.y - 4) * 4][threadIdx.x + 16] 	= reg_tile[0][1];
        sm_b[1 + (threadIdx.y - 4) * 4][threadIdx.x + 16] 	= reg_tile[1][1];
        sm_b[2 + (threadIdx.y - 4) * 4][threadIdx.x + 16] 	= reg_tile[2][1];
        sm_b[3 + (threadIdx.y - 4) * 4][threadIdx.x + 16] 	= reg_tile[3][1];

        sm_b[0 + (threadIdx.y - 4) * 4][threadIdx.x + 32] 	= reg_tile[0][2];
        sm_b[1 + (threadIdx.y - 4) * 4][threadIdx.x + 32] 	= reg_tile[1][2];
        sm_b[2 + (threadIdx.y - 4) * 4][threadIdx.x + 32] 	= reg_tile[2][2];
        sm_b[3 + (threadIdx.y - 4) * 4][threadIdx.x + 32] 	= reg_tile[3][2];

        sm_b[0 + (threadIdx.y - 4) * 4][threadIdx.x + 48] 	= reg_tile[0][3];
        sm_b[1 + (threadIdx.y - 4) * 4][threadIdx.x + 48] 	= reg_tile[1][3];
        sm_b[2 + (threadIdx.y - 4) * 4][threadIdx.x + 48] 	= reg_tile[2][3];
        sm_b[3 + (threadIdx.y - 4) * 4][threadIdx.x + 48] 	= reg_tile[3][3];
	}
	__syncthreads();

	if (threadIdx.y < 4)						// 0, 1, 2, 3
	{
		reg_tile[0][0] = sm_a[threadIdx.y + 0][(threadIdx.x)];
		reg_tile[1][0] = sm_a[threadIdx.y + 4][(threadIdx.x)];
		reg_tile[2][0] = sm_a[threadIdx.y + 8][(threadIdx.x)];
		reg_tile[3][0] = sm_a[threadIdx.y + 12][(threadIdx.x)];

		reg_tile[0][1] = sm_a[threadIdx.y + 0][(threadIdx.x) + 16];
		reg_tile[1][1] = sm_a[threadIdx.y + 4][(threadIdx.x) + 16];
		reg_tile[2][1] = sm_a[threadIdx.y + 8][(threadIdx.x) + 16];
        reg_tile[3][1] = sm_a[threadIdx.y + 12][(threadIdx.x) + 16];
        
		reg_tile[0][2] = sm_a[threadIdx.y + 0][(threadIdx.x) + 32];
		reg_tile[1][2] = sm_a[threadIdx.y + 4][(threadIdx.x) + 32];
		reg_tile[2][2] = sm_a[threadIdx.y + 8][(threadIdx.x) + 32];
		reg_tile[3][2] = sm_a[threadIdx.y + 12][(threadIdx.x) + 32];
        
        reg_tile[0][3] = sm_a[threadIdx.y + 0][(threadIdx.x) + 48];
        reg_tile[1][3] = sm_a[threadIdx.y + 4][(threadIdx.x) + 48];
        reg_tile[2][3] = sm_a[threadIdx.y + 8][(threadIdx.x) + 48];
        reg_tile[3][3] = sm_a[threadIdx.y + 12][(threadIdx.x) + 48];
	}

	if (threadIdx.y >= 4 && threadIdx.y < 8)	// 4, 5, 6, 7
	{
		reg_tile[0][0] = sm_b[(threadIdx.y - 4) + 0][(threadIdx.x)];
		reg_tile[1][0] = sm_b[(threadIdx.y - 4) + 4][(threadIdx.x)];
		reg_tile[2][0] = sm_b[(threadIdx.y - 4) + 8][(threadIdx.x)];
		reg_tile[3][0] = sm_b[(threadIdx.y - 4) + 12][(threadIdx.x)];

		reg_tile[0][1] = sm_b[(threadIdx.y - 4) + 0][(threadIdx.x) + 16];
        reg_tile[1][1] = sm_b[(threadIdx.y - 4) + 4][(threadIdx.x) + 16];
		reg_tile[2][1] = sm_b[(threadIdx.y - 4) + 8][(threadIdx.x) + 16];
		reg_tile[3][1] = sm_b[(threadIdx.y - 4) + 12][(threadIdx.x) + 16];
		
		reg_tile[0][2] = sm_b[(threadIdx.y - 4) + 0][(threadIdx.x) + 32];
		reg_tile[1][2] = sm_b[(threadIdx.y - 4) + 4][(threadIdx.x) + 32];
		reg_tile[2][2] = sm_b[(threadIdx.y - 4) + 8][(threadIdx.x) + 32];
		reg_tile[3][2] = sm_b[(threadIdx.y - 4) + 12][(threadIdx.x) + 32];
        
		reg_tile[0][3] = sm_b[(threadIdx.y - 4) + 0][(threadIdx.x) + 48];	
        reg_tile[1][3] = sm_b[(threadIdx.y - 4) + 4][(threadIdx.x) + 48];
        reg_tile[2][3] = sm_b[(threadIdx.y - 4) + 8][(threadIdx.x) + 48];
        reg_tile[3][3] = sm_b[(threadIdx.y - 4) + 12][(threadIdx.x) + 48];
	}
	__syncthreads();

	if (threadIdx.y >= 8 && threadIdx.y < 12)	// 8, 9, 10, 11
	{
		sm_a[0 + (threadIdx.y - 8) * 4][threadIdx.x] = reg_tile[0][0];
        sm_a[1 + (threadIdx.y - 8) * 4][threadIdx.x] = reg_tile[1][0];
        sm_a[2 + (threadIdx.y - 8) * 4][threadIdx.x] = reg_tile[2][0];
        sm_a[3 + (threadIdx.y - 8) * 4][threadIdx.x] = reg_tile[3][0];

        sm_a[0 + (threadIdx.y - 8) * 4][threadIdx.x + 16] = reg_tile[0][1];
        sm_a[1 + (threadIdx.y - 8) * 4][threadIdx.x + 16] = reg_tile[1][1];
        sm_a[2 + (threadIdx.y - 8) * 4][threadIdx.x + 16] = reg_tile[2][1];
        sm_a[3 + (threadIdx.y - 8) * 4][threadIdx.x + 16] = reg_tile[3][1];

        sm_a[0 + (threadIdx.y - 8) * 4][threadIdx.x + 32] = reg_tile[0][2];
        sm_a[1 + (threadIdx.y - 8) * 4][threadIdx.x + 32] = reg_tile[1][2];
        sm_a[2 + (threadIdx.y - 8) * 4][threadIdx.x + 32] = reg_tile[2][2];
        sm_a[3 + (threadIdx.y - 8) * 4][threadIdx.x + 32] = reg_tile[3][2];

        sm_a[0 + (threadIdx.y - 8) * 4][threadIdx.x + 48] = reg_tile[0][3];
        sm_a[1 + (threadIdx.y - 8) * 4][threadIdx.x + 48] = reg_tile[1][3];
        sm_a[2 + (threadIdx.y - 8) * 4][threadIdx.x + 48] = reg_tile[2][3];
        sm_a[3 + (threadIdx.y - 8) * 4][threadIdx.x + 48] = reg_tile[3][3];
	}

	if (threadIdx.y >= 12)	// 12, 13, 14, 15
	{
		sm_b[0 + (threadIdx.y - 12) * 4][threadIdx.x] = reg_tile[0][0];
        sm_b[1 + (threadIdx.y - 12) * 4][threadIdx.x] = reg_tile[1][0];
        sm_b[2 + (threadIdx.y - 12) * 4][threadIdx.x] = reg_tile[2][0];
		sm_b[3 + (threadIdx.y - 12) * 4][threadIdx.x] = reg_tile[3][0];
		
        sm_b[0 + (threadIdx.y - 12) * 4][threadIdx.x + 16] = reg_tile[0][1];
        sm_b[1 + (threadIdx.y - 12) * 4][threadIdx.x + 16] = reg_tile[1][1];
        sm_b[2 + (threadIdx.y - 12) * 4][threadIdx.x + 16] = reg_tile[2][1];
		sm_b[3 + (threadIdx.y - 12) * 4][threadIdx.x + 16] = reg_tile[3][1];
		
        sm_b[0 + (threadIdx.y - 12) * 4][threadIdx.x + 32] = reg_tile[0][2];
        sm_b[1 + (threadIdx.y - 12) * 4][threadIdx.x + 32] = reg_tile[1][2];
        sm_b[2 + (threadIdx.y - 12) * 4][threadIdx.x + 32] = reg_tile[2][2];
		sm_b[3 + (threadIdx.y - 12) * 4][threadIdx.x + 32] = reg_tile[3][2];
		
        sm_b[0 + (threadIdx.y - 12) * 4][threadIdx.x + 48] = reg_tile[0][3];
        sm_b[1 + (threadIdx.y - 12) * 4][threadIdx.x + 48] = reg_tile[1][3];
        sm_b[2 + (threadIdx.y - 12) * 4][threadIdx.x + 48] = reg_tile[2][3];
        sm_b[3 + (threadIdx.y - 12) * 4][threadIdx.x + 48] = reg_tile[3][3];
	}
	__syncthreads();

	if (threadIdx.y >= 8 && threadIdx.y < 12)	// 8, 9, 10, 11
	{
		reg_tile[0][0] = sm_a[(threadIdx.y - 8) + 0][(threadIdx.x)];
		reg_tile[1][0] = sm_a[(threadIdx.y - 8) + 4][(threadIdx.x)];
		reg_tile[2][0] = sm_a[(threadIdx.y - 8) + 8][(threadIdx.x)];
		reg_tile[3][0] = sm_a[(threadIdx.y - 8) + 12][(threadIdx.x)];

		reg_tile[0][1] = sm_a[(threadIdx.y - 8) + 0][(threadIdx.x) + 16];
		reg_tile[1][1] = sm_a[(threadIdx.y - 8) + 4][(threadIdx.x) + 16];
		reg_tile[2][1] = sm_a[(threadIdx.y - 8) + 8][(threadIdx.x) + 16];
		reg_tile[3][1] = sm_a[(threadIdx.y - 8) + 12][(threadIdx.x) + 16];

		reg_tile[0][2] = sm_a[(threadIdx.y - 8) + 0][(threadIdx.x) + 32];
		reg_tile[1][2] = sm_a[(threadIdx.y - 8) + 4][(threadIdx.x) + 32];		
        reg_tile[2][2] = sm_a[(threadIdx.y - 8) + 8][(threadIdx.x) + 32];
        reg_tile[3][2] = sm_a[(threadIdx.y - 8) + 12][(threadIdx.x) + 32];

        reg_tile[0][3] = sm_a[(threadIdx.y - 8) + 0][(threadIdx.x) + 48];
        reg_tile[1][3] = sm_a[(threadIdx.y - 8) + 4][(threadIdx.x) + 48];
        reg_tile[2][3] = sm_a[(threadIdx.y - 8) + 8][(threadIdx.x) + 48];
		reg_tile[3][3] = sm_a[(threadIdx.y - 8) + 12][(threadIdx.x) + 48];
	}

	if (threadIdx.y >= 12)	// 12, 13, 14, 15
	{
		reg_tile[0][0] = sm_b[(threadIdx.y - 12) + 0][(threadIdx.x)];
		reg_tile[1][0] = sm_b[(threadIdx.y - 12) + 4][(threadIdx.x)];
		reg_tile[2][0] = sm_b[(threadIdx.y - 12) + 8][(threadIdx.x)];
		reg_tile[3][0] = sm_b[(threadIdx.y - 12) + 12][(threadIdx.x)];

		reg_tile[0][1] = sm_b[(threadIdx.y - 12) + 0][(threadIdx.x) + 16];
		reg_tile[1][1] = sm_b[(threadIdx.y - 12) + 4][(threadIdx.x) + 16];
		reg_tile[2][1] = sm_b[(threadIdx.y - 12) + 8][(threadIdx.x) + 16];
		reg_tile[3][1] = sm_b[(threadIdx.y - 12) + 12][(threadIdx.x) + 16];

		reg_tile[0][2] = sm_b[(threadIdx.y - 12) + 0][(threadIdx.x) + 32];
		reg_tile[1][2] = sm_b[(threadIdx.y - 12) + 4][(threadIdx.x) + 32];
		reg_tile[2][2] = sm_b[(threadIdx.y - 12) + 8][(threadIdx.x) + 32];
        reg_tile[3][2] = sm_b[(threadIdx.y - 12) + 12][(threadIdx.x) + 32];

		reg_tile[0][3] = sm_b[(threadIdx.y - 12) + 0][(threadIdx.x) + 48];
		reg_tile[1][3] = sm_b[(threadIdx.y - 12) + 4][(threadIdx.x) + 48];
		reg_tile[2][3] = sm_b[(threadIdx.y - 12) + 8][(threadIdx.x) + 48];
        reg_tile[3][3] = sm_b[(threadIdx.y - 12) + 12][(threadIdx.x) + 48];
	}
	__syncthreads();

	//
	//	sd2_1, 2, 3, 4, 5 and 6.
	//
	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_1 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
            sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_1[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h1) * size_h1) * size_p4) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_1[(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 48];
            
			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H3 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_2 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
            sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_2[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h2) * size_p4) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_2[(blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 0];
            temp_bv[1] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 16];
            temp_bv[2] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 32];
            temp_bv[3] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_3 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_3[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h1) * size_p4) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_3[(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_4 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_4[(blk_idx_p5 * SIZE_SLICE_1_P5 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h1) * size_h1) * size_p5) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_4[(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_5 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_5[(blk_idx_p5 * SIZE_SLICE_1_P5 + ll + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h2) * size_p5) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_5[(blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_6 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_6[(blk_idx_p5 * SIZE_SLICE_1_P6 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h1) * size_p5) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_6[(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 0];
            temp_bv[1] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 16];
            temp_bv[2] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 32];
            temp_bv[3] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 48];

			for (int xx = 0; xx < 4; xx++)	// 4 -> rng_p4: Local Transactions...
			{
                temp_av = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// Store Results (Registers) to Global Memory
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void kernel_ccsdT_sd2_123456_partial_partial(double* t3, double* d_t2_1, double* d_t2_2, double* d_t2_3, double* d_t2_4, double* d_t2_5, double* d_t2_6, double* d_v2_1, double* d_v2_2, double* d_v2_3, double* d_v2_4, double* d_v2_5, double* d_v2_6, int size_h3, int size_h2, int size_h1, int size_p6, int size_p5, int size_p4, int size_p7, int numBlk_h3, int numBlk_h2, int numBlk_h1, int numBlk_p6, int numBlk_p5, int numBlk_p4, int kernel_1, int kernel_2, int kernel_3, int kernel_4, int kernel_5, int kernel_6, int kernel_7, int kernel_8, int kernel_9, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int internal_upperbound   	= 0;
	int internal_offset;

	// should support for non-full tiles
	int idx_h3 = threadIdx.x % SIZE_SLICE_1_H3;
	int idx_h2 = threadIdx.x / SIZE_SLICE_1_H3;
	int idx_p6 = threadIdx.y % SIZE_SLICE_1_P6;
	int idx_h1 = threadIdx.y / SIZE_SLICE_1_P6;

	// Common for Threads within a Thread Block
	int tmp_blkIdx;        
    int blk_idx_p4  = blockIdx.x / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);
    tmp_blkIdx      = blockIdx.x % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);

    int blk_idx_p5  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);

    int blk_idx_p6  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1);

    int blk_idx_h1 = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2);
    tmp_blkIdx     = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2);

    int blk_idx_h2 = (tmp_blkIdx) / (numBlk_h3);
    int blk_idx_h3 = (tmp_blkIdx) % (numBlk_h3);

    int rng_h3, rng_h2, rng_h1, rng_p6, rng_p5, rng_p4;

    if ((size_h3 - (blk_idx_h3 * SIZE_SLICE_2_H3)) >= SIZE_SLICE_2_H3)
    {
        rng_h3 = SIZE_SLICE_2_H3;
    }
    else
    {
        rng_h3 = size_h3 % SIZE_SLICE_2_H3;
    }
    
    if ((size_h2 - (blk_idx_h2 * SIZE_SLICE_2_H2)) >= SIZE_SLICE_2_H2)
    {
        rng_h2 = SIZE_SLICE_2_H2;
    }
    else
    {
        rng_h2 = size_h2 % SIZE_SLICE_2_H2;
    }

    if ((size_h1 - (blk_idx_h1 * SIZE_SLICE_2_H1)) >= SIZE_SLICE_2_H1)
    {
        rng_h1 = SIZE_SLICE_2_H1;
    }
    else
    {
        rng_h1 = size_h1 % SIZE_SLICE_2_H1;
    }
    
    if ((size_p6 - (blk_idx_p6 * SIZE_SLICE_2_P6)) >= SIZE_SLICE_2_P6)
    {
        rng_p6 = SIZE_SLICE_2_P6;
    }
    else
    {
        rng_p6 = size_p6 % SIZE_SLICE_2_P6;
    }

    if ((size_p5 - (blk_idx_p5 * SIZE_SLICE_2_P5)) >= SIZE_SLICE_2_P5)
    {
        rng_p5 = SIZE_SLICE_2_P5;
    }
    else
    {
        rng_p5 = size_p5 % SIZE_SLICE_2_P5;
    }

    if ((size_p4 - (blk_idx_p4 * SIZE_SLICE_2_P4)) >= SIZE_SLICE_2_P4)
    {
        rng_p4 = SIZE_SLICE_2_P4;
    }
    else
    {
        rng_p4 = size_p4 % SIZE_SLICE_2_P4;
    }

    int t3_base_thread = blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h3 + 
                        (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h2 + 
                        (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_h1 + 
                        (blk_idx_p6 * SIZE_SLICE_2_P6 + idx_p6 + 
                        (blk_idx_p5 * SIZE_SLICE_2_P5 + 
                        (blk_idx_p4 * SIZE_SLICE_2_P4) * size_p5) * size_p6) * size_h1) * size_h2) * size_h3;


	double temp_av;
	double temp_bv[4];
	double reg_tile[4][4];

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_1 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h1 && idx_h1 < rng_h2 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_1[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h1) * size_h1) * size_p4) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h3 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_1[(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{   
            temp_bv[0] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 48];
            
			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H3 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_2 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h2 && idx_h1 < rng_h3 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_2[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h2) * size_p4) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h1 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_2[(blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 0];
            temp_bv[1] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 16];
            temp_bv[2] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 32];
            temp_bv[3] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_3 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h1 && idx_h1 < rng_h3 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_3[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h1) * size_p4) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h2 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_3[(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_4 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h1 && idx_h1 < rng_h2 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_4[(blk_idx_p5 * SIZE_SLICE_1_P5 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h1) * size_h1) * size_p5) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h3 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_4[(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_5 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h2 && idx_h1 < rng_h3 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_5[(blk_idx_p5 * SIZE_SLICE_1_P5 + ll + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h2) * size_p5) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h1 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_5[(blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_6 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p6 < rng_h1 && idx_h1 < rng_h3 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_6[(blk_idx_p5 * SIZE_SLICE_1_P6 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h1) * size_p5) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h2 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_6[(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 0];
            temp_bv[1] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 16];
            temp_bv[2] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 32];
            temp_bv[3] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 48];

			for (int xx = 0; xx < 4; xx++)	// 4 -> rng_p4: Local Transactions...
			{
                temp_av = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_h3 < rng_h3 && idx_h2 < rng_h2 && idx_p6 < rng_p6 && idx_h1 < rng_h1)
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_p4 && j < rng_p5)
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

__global__ void kernel_ccsdT_sd2_123456_partial_full(double* t3, double* d_t2_1, double* d_t2_2, double* d_t2_3, double* d_t2_4, double* d_t2_5, double* d_t2_6, double* d_v2_1, double* d_v2_2, double* d_v2_3, double* d_v2_4, double* d_v2_5, double* d_v2_6, int size_h3, int size_h2, int size_h1, int size_p6, int size_p5, int size_p4, int size_p7, int numBlk_h3, int numBlk_h2, int numBlk_h1, int numBlk_p6, int numBlk_p5, int numBlk_p4, int kernel_1, int kernel_2, int kernel_3, int kernel_4, int kernel_5, int kernel_6, int kernel_7, int kernel_8, int kernel_9, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int internal_upperbound   = 0;
    int internal_offset;
    
    // should support for non-full tiles
	int idx_h3 = threadIdx.x % SIZE_SLICE_1_H3;
	int idx_h2 = threadIdx.x / SIZE_SLICE_1_H3;
	int idx_p6 = threadIdx.y % SIZE_SLICE_1_P6;
	int idx_h1 = threadIdx.y / SIZE_SLICE_1_P6;

	// Common for Threads within a Thread Block
	int tmp_blkIdx;        
    int blk_idx_p4  = blockIdx.x / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);
    tmp_blkIdx      = blockIdx.x % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);

    int blk_idx_p5  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);

    int blk_idx_p6  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1);

    int blk_idx_h1 = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2);
    tmp_blkIdx     = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2);

    int blk_idx_h2 = (tmp_blkIdx) / (numBlk_h3);
    int blk_idx_h3 = (tmp_blkIdx) % (numBlk_h3);

    int t3_base_thread = blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h3 + 
                        (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h2 + 
                        (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_h1 + 
                        (blk_idx_p6 * SIZE_SLICE_2_P6 + idx_p6 + 
                        (blk_idx_p5 * SIZE_SLICE_2_P5 + 
                        (blk_idx_p4 * SIZE_SLICE_2_P4) * size_p5) * size_p6) * size_h1) * size_h2) * size_h3;

	double temp_av;
	double temp_bv[4];
	double reg_tile[4][4];

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_1 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_1[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h1) * size_h1) * size_p4) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_1[(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 48];
            
			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H3 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_2 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_2[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h2) * size_p4) * size_p7 + (threadIdx.x + l)];		
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_2[(blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 0];
            temp_bv[1] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 16];
            temp_bv[2] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 32];
            temp_bv[3] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_3 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_3[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h1) * size_p4) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_3[(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_4 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_4[(blk_idx_p5 * SIZE_SLICE_1_P5 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h1) * size_h1) * size_p5) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_4[(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_5 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_5[(blk_idx_p5 * SIZE_SLICE_1_P5 + ll + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h2) * size_p5) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_5[(blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_6 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_6[(blk_idx_p5 * SIZE_SLICE_1_P6 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h1) * size_p5) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_6[(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 0];
            temp_bv[1] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 16];
            temp_bv[2] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 32];
            temp_bv[3] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 48];

			for (int xx = 0; xx < 4; xx++)	// 4 -> rng_p4: Local Transactions...
			{
                temp_av = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

__global__ void kernel_ccsdT_sd2_123456_full_full(double* t3, double* d_t2_1, double* d_t2_2, double* d_t2_3, double* d_t2_4, double* d_t2_5, double* d_t2_6, double* d_v2_1, double* d_v2_2, double* d_v2_3, double* d_v2_4, double* d_v2_5, double* d_v2_6, int size_h3, int size_h2, int size_h1, int size_p6, int size_p5, int size_p4, int size_p7, int numBlk_h3, int numBlk_h2, int numBlk_h1, int numBlk_p6, int numBlk_p5, int numBlk_p4, int kernel_1, int kernel_2, int kernel_3, int kernel_4, int kernel_5, int kernel_6, int kernel_7, int kernel_8, int kernel_9, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

    // should support for non-full tiles
	int idx_h3 = threadIdx.x % SIZE_SLICE_1_H3;
	int idx_h2 = threadIdx.x / SIZE_SLICE_1_H3;
	int idx_p6 = threadIdx.y % SIZE_SLICE_1_P6;
	int idx_h1 = threadIdx.y / SIZE_SLICE_1_P6;

    // Common for Threads within a Thread Block
	int tmp_blkIdx;        
    int blk_idx_p4  = blockIdx.x / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);
    tmp_blkIdx      = blockIdx.x % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);

    int blk_idx_p5  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);

    int blk_idx_p6  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1);

    int blk_idx_h1 = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2);
    tmp_blkIdx     = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2);

    int blk_idx_h2 = (tmp_blkIdx) / (numBlk_h3);
    int blk_idx_h3 = (tmp_blkIdx) % (numBlk_h3);

    int t3_base_thread = blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h3 + 
                        (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h2 + 
                        (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_h1 + 
                        (blk_idx_p6 * SIZE_SLICE_2_P6 + idx_p6 + 
                        (blk_idx_p5 * SIZE_SLICE_2_P5 + 
                        (blk_idx_p4 * SIZE_SLICE_2_P4) * size_p5) * size_p6) * size_h1) * size_h2) * size_h3;

	double temp_av;
	double temp_bv[4];
	double reg_tile[4][4];

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_1 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_1[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h1) * size_h1) * size_p4) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_1[(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + 48];
            
			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H3 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_2 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
            sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_2[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h2) * size_p4) * size_p7 + (threadIdx.x + l)];		
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_2[(blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 0];
            temp_bv[1] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 16];
            temp_bv[2] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 32];
            temp_bv[3] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_3 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
            sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_3[(blk_idx_p4 * SIZE_SLICE_1_P4 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h1) * size_p4) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_3[(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_1_P5 + ll) * size_p6) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_4 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
            sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_4[(blk_idx_p5 * SIZE_SLICE_1_P5 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h1) * size_h1) * size_p5) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_4[(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_b[ll][idx_h3 + (idx_p6) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_5 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
            sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_5[(blk_idx_p5 * SIZE_SLICE_1_P5 + ll + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h2) * size_p5) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_5[(blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 0];
            temp_bv[1] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 16];
            temp_bv[2] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 32];
            temp_bv[3] = sm_b[ll][idx_h1 + (idx_p6) * SIZE_SLICE_1_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_1_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_6 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_6[(blk_idx_p5 * SIZE_SLICE_1_P6 + ll + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_p6 + (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h1) * size_h1) * size_p5) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_6[(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_p6 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_h1 + (blk_idx_p4 * SIZE_SLICE_1_P4 + ll) * size_p6) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 0];
            temp_bv[1] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 16];
            temp_bv[2] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 32];
            temp_bv[3] = sm_b[ll][idx_h2 + (idx_p6) * SIZE_SLICE_1_H2 + 48];

			for (int xx = 0; xx < 4; xx++)	// 4 -> rng_p4: Local Transactions...
			{
                temp_av = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_1_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] = reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void kernel_ccsdT_sd2_789_partial_partial(double* t3, double* d_t2_7, double* d_t2_8, double* d_t2_9, double* d_v2_7, double* d_v2_8, double* d_v2_9, int size_h3, int size_h2, int size_h1, int size_p6, int size_p5, int size_p4, int size_p7, int numBlk_h3, int numBlk_h2, int numBlk_h1, int numBlk_p6, int numBlk_p5, int numBlk_p4, int kernel_1, int kernel_2, int kernel_3, int kernel_4, int kernel_5, int kernel_6, int kernel_7, int kernel_8, int kernel_9, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int internal_upperbound   = 0;
	int internal_offset;

	// should support for non-full tiles
	int idx_h3 = threadIdx.x % SIZE_SLICE_2_H3;
	int idx_h2 = threadIdx.x / SIZE_SLICE_2_H3;
	int idx_p4 = threadIdx.y % SIZE_SLICE_2_P4;
	int idx_h1 = threadIdx.y / SIZE_SLICE_2_P4;

    int tmp_blkIdx;        
    int blk_idx_p4  = blockIdx.x / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);
    tmp_blkIdx      = blockIdx.x % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);

    int blk_idx_p5  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);

    int blk_idx_p6  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1);

    int blk_idx_h1 = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2);
    tmp_blkIdx     = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2);

    int blk_idx_h2 = (tmp_blkIdx) / (numBlk_h3);
    int blk_idx_h3 = (tmp_blkIdx) % (numBlk_h3);

    int rng_h3, rng_h2, rng_h1, rng_p6, rng_p5, rng_p4;

    if ((size_h3 - (blk_idx_h3 * SIZE_SLICE_2_H3)) >= SIZE_SLICE_2_H3)
    {
        rng_h3 = SIZE_SLICE_2_H3;
    }
    else
    {
        rng_h3 = size_h3 % SIZE_SLICE_2_H3;
    }
    
    if ((size_h2 - (blk_idx_h2 * SIZE_SLICE_2_H2)) >= SIZE_SLICE_2_H2)
    {
        rng_h2 = SIZE_SLICE_2_H2;
    }
    else
    {
        rng_h2 = size_h2 % SIZE_SLICE_2_H2;
    }

    if ((size_h1 - (blk_idx_h1 * SIZE_SLICE_2_H1)) >= SIZE_SLICE_2_H1)
    {
        rng_h1 = SIZE_SLICE_2_H1;
    }
    else
    {
        rng_h1 = size_h1 % SIZE_SLICE_2_H1;
    }
    
    if ((size_p6 - (blk_idx_p6 * SIZE_SLICE_2_P6)) >= SIZE_SLICE_2_P6)
    {
        rng_p6 = SIZE_SLICE_2_P6;
    }
    else
    {
        rng_p6 = size_p6 % SIZE_SLICE_2_P6;
    }

    if ((size_p5 - (blk_idx_p5 * SIZE_SLICE_2_P5)) >= SIZE_SLICE_2_P5)
    {
        rng_p5 = SIZE_SLICE_2_P5;
    }
    else
    {
        rng_p5 = size_p5 % SIZE_SLICE_2_P5;
    }

    if ((size_p4 - (blk_idx_p4 * SIZE_SLICE_2_P4)) >= SIZE_SLICE_2_P4)
    {
        rng_p4 = SIZE_SLICE_2_P4;
    }
    else
    {
        rng_p4 = size_p4 % SIZE_SLICE_2_P4;
    }

    int t3_base_thread = blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h3 + 
                        (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h2 + 
                        (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_h1 + 
                        (blk_idx_p6 * SIZE_SLICE_2_P6 +  
                        (blk_idx_p5 * SIZE_SLICE_2_P5 + 
                        (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_p4) * size_p5) * size_p6) * size_h1) * size_h2) * size_h3;


	double temp_av;
	double temp_bv[4];
	double reg_tile[4][4];

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_7 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		if (idx_p4 < rng_h1 && idx_h1 < rng_h2 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p6; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_7[(blk_idx_p6 *  SIZE_SLICE_2_P6 + ll + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p4 + (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h1) * size_h1) * size_p6) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p4 < rng_h3 && idx_h1 < rng_p4 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_7[(blk_idx_h3 * SIZE_SLICE_2_H3 + idx_p4 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_h1) * size_p5) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h3 + (idx_p4) * SIZE_SLICE_2_H3 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_8 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
        if (idx_p4 < rng_h2 && idx_h1 < rng_h3 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p6; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_8[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_p4 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h1) * size_h2) * size_p6) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
        if (idx_p4 < rng_h1 && idx_h1 < rng_p4 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_8[(blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p4 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_1_P4 + idx_h1) * size_p5) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 0];
            temp_bv[1] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 16];
            temp_bv[2] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 32];
            temp_bv[3] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h1 + (idx_p4) * SIZE_SLICE_2_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_9 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
        if (idx_p4 < rng_h1 && idx_h1 < rng_h3 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p6; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_9[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p4 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h1) * size_h1) * size_p6) * size_p7 + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
        if (idx_p4 < rng_h2 && idx_h1 < rng_p4 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_9[(blk_idx_h2 * SIZE_SLICE_2_H2 + idx_p4 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_h1) * size_p5) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
            temp_bv[0] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h2 + (idx_p4) * SIZE_SLICE_2_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}


	// Store Results (Registers) to Global Memory
	// Part: Generalized Threads
	// Part: Generalized Register-Tiling
	if (idx_h3 < rng_h3 && idx_h2 < rng_h2 && idx_p4 < rng_p4 && idx_h1 < rng_h1)
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if(i < rng_p6 && j < rng_p5)
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] += reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void kernel_ccsdT_sd2_789_partial_full(double* t3, double* d_t2_7, double* d_t2_8, double* d_t2_9, double* d_v2_7, double* d_v2_8, double* d_v2_9, int size_h3, int size_h2, int size_h1, int size_p6, int size_p5, int size_p4, int size_p7, int numBlk_h3, int numBlk_h2, int numBlk_h1, int numBlk_p6, int numBlk_p5, int numBlk_p4, int kernel_1, int kernel_2, int kernel_3, int kernel_4, int kernel_5, int kernel_6, int kernel_7, int kernel_8, int kernel_9, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int internal_upperbound   = 0;
	int internal_offset;

    // should support for non-full tiles
	int idx_h3 = threadIdx.x % SIZE_SLICE_2_H3;
	int idx_h2 = threadIdx.x / SIZE_SLICE_2_H3;
	int idx_p4 = threadIdx.y % SIZE_SLICE_2_P4;
	int idx_h1 = threadIdx.y / SIZE_SLICE_2_P4;

    int tmp_blkIdx;        
    int blk_idx_p4  = blockIdx.x / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);
    tmp_blkIdx      = blockIdx.x % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);

    int blk_idx_p5  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);

    int blk_idx_p6  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1);

    int blk_idx_h1 = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2);
    tmp_blkIdx     = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2);

    int blk_idx_h2 = (tmp_blkIdx) / (numBlk_h3);
    int blk_idx_h3 = (tmp_blkIdx) % (numBlk_h3);

    int t3_base_thread = blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h3 + 
                        (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h2 + 
                        (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_h1 + 
                        (blk_idx_p6 * SIZE_SLICE_2_P6 +  
                        (blk_idx_p5 * SIZE_SLICE_2_P5 + 
                        (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_p4) * size_p5) * size_p6) * size_h1) * size_h2) * size_h3;

	double temp_av;
	double temp_bv[4];
	double reg_tile[4][4];

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_7 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_7[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p4 + (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h1) * size_h1) * size_p6) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_7[(blk_idx_h3 * SIZE_SLICE_2_H3 + idx_p4 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_h1) * size_p5) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h3 + (idx_p4) * SIZE_SLICE_2_H3 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_8 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_8[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_p4 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h1) * size_h2) * size_p6) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_8[(blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p4 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_1_P4 + idx_h1) * size_p5) * size_h1) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 0];
            temp_bv[1] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 16];
            temp_bv[2] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 32];
            temp_bv[3] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h1 + (idx_p4) * SIZE_SLICE_2_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	internal_upperbound = 0;
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_9 == 1; l+= SIZE_INT_UNIT)
	{
		// Part: Generalized Contraction Index (p7b)
		internal_offset = (l + SIZE_INT_UNIT) - size_internal;
		if (internal_offset > 0) internal_upperbound = internal_offset;

		// Load Input Tensor to Shared Memory: 16:16
        // # of Internal Indices: 1
        if (threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_9[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p4 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h1) * size_h1) * size_p6) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_9[(blk_idx_h2 * SIZE_SLICE_2_H2 + idx_p4 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_h1) * size_p5) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h2 + (idx_p4) * SIZE_SLICE_2_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// Store Results (Registers) to Global Memory
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] += reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void kernel_ccsdT_sd2_789_full_full(double* t3, double* d_t2_7, double* d_t2_8, double* d_t2_9, double* d_v2_7, double* d_v2_8, double* d_v2_9, int size_h3, int size_h2, int size_h1, int size_p6, int size_p5, int size_p4, int size_p7, int numBlk_h3, int numBlk_h2, int numBlk_h1, int numBlk_p6, int numBlk_p5, int numBlk_p4, int kernel_1, int kernel_2, int kernel_3, int kernel_4, int kernel_5, int kernel_6, int kernel_7, int kernel_8, int kernel_9, int stride_reg_x, int stride_reg_y, int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

    // should support for non-full tiles
	int idx_h3 = threadIdx.x % SIZE_SLICE_2_H3;
	int idx_h2 = threadIdx.x / SIZE_SLICE_2_H3;
	int idx_p4 = threadIdx.y % SIZE_SLICE_2_P4;
	int idx_h1 = threadIdx.y / SIZE_SLICE_2_P4;

    int tmp_blkIdx;        
    int blk_idx_p4  = blockIdx.x / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);
    tmp_blkIdx      = blockIdx.x % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6 * numBlk_p5);

    int blk_idx_p5  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1 * numBlk_p6);

    int blk_idx_p6  = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2 * numBlk_h1);
    tmp_blkIdx      = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2 * numBlk_h1);

    int blk_idx_h1 = (tmp_blkIdx) / (numBlk_h3 * numBlk_h2);
    tmp_blkIdx     = (tmp_blkIdx) % (numBlk_h3 * numBlk_h2);

    int blk_idx_h2 = (tmp_blkIdx) / (numBlk_h3);
    int blk_idx_h3 = (tmp_blkIdx) % (numBlk_h3);

    int t3_base_thread = blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h3 + 
                        (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h2 + 
                        (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_h1 + 
                        (blk_idx_p6 * SIZE_SLICE_2_P6 +  
                        (blk_idx_p5 * SIZE_SLICE_2_P5 + 
                        (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_p4) * size_p5) * size_p6) * size_h1) * size_h2) * size_h3;


	double temp_av;
	double temp_bv[4];
	double reg_tile[4][4];

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
	reg_tile[i][j] = 0.0;

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_7 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_7[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p4 + (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h1) * size_h1) * size_p6) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_7[(blk_idx_h3 * SIZE_SLICE_2_H3 + idx_p4 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_h1) * size_p5) * size_h3) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h2) * SIZE_SLICE_2_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h3 + (idx_p4) * SIZE_SLICE_2_H3 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_8 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
            sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_8[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_p4 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h1) * size_h2) * size_p6) * size_p7 + (threadIdx.x + l)];
            sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_8[(blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p4 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_1_P4 + idx_h1) * size_p5) * size_h1) * size_p7 + (threadIdx.x + l)];
        }
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 0];
            temp_bv[1] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 16];
            temp_bv[2] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 32];
            temp_bv[3] = sm_a[ll][idx_h2 + (idx_h3) * SIZE_SLICE_2_H2 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h1 + (idx_p4) * SIZE_SLICE_2_H1 + (xx * 16)];

				reg_tile[0][xx] -= temp_av * temp_bv[0];
				reg_tile[1][xx] -= temp_av * temp_bv[1];
				reg_tile[2][xx] -= temp_av * temp_bv[2];
				reg_tile[3][xx] -= temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// tensor contraction
	#pragma unroll 1
	for (int l = 0; l < size_internal && kernel_9 == 1; l+= SIZE_INT_UNIT)
	{
		// Load Input Tensor to Shared Memory: 16:16
		// # of Internal Indices: 1
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_9[(blk_idx_p6 * SIZE_SLICE_2_P6 + ll + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_p4 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h1) * size_h1) * size_p6) * size_p7 + (threadIdx.x + l)];
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_9[(blk_idx_h2 * SIZE_SLICE_2_H2 + idx_p4 + (blk_idx_p5 * SIZE_SLICE_2_P5 + ll + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_h1) * size_p5) * size_h2) * size_p7 + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 0];
            temp_bv[1] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 16];
            temp_bv[2] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 32];
            temp_bv[3] = sm_a[ll][idx_h1 + (idx_h3) * SIZE_SLICE_2_H1 + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
                temp_av = sm_b[ll][idx_h2 + (idx_p4) * SIZE_SLICE_2_H2 + (xx * 16)];

				reg_tile[0][xx] += temp_av * temp_bv[0];
				reg_tile[1][xx] += temp_av * temp_bv[1];
				reg_tile[2][xx] += temp_av * temp_bv[2];
				reg_tile[3][xx] += temp_av * temp_bv[3];
			}
		}
		__syncthreads();
	}

	// Store Results (Registers) to Global Memory
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] += reg_tile[i][j];
		}
	}
}

//
// 	Inputs: Problem Sizes and Tensors 
// 	Temporally, input tensors are host-side memory.
//
extern "C" 
void sd_t_d2_all_cuda(int size_h3, int size_h2, int size_h1, int size_p6, int size_p5, int size_p4, int size_p7, 
	double* t3, 
	double* t2_1, double* v2_1,
	double* t2_2, double* v2_2,
	double* t2_3, double* v2_3,
	double* t2_4, double* v2_4,
	double* t2_5, double* v2_5,
	double* t2_6, double* v2_6,
	double* t2_7, double* v2_7,
	double* t2_8, double* v2_8,
    double* t2_9, double* v2_9, 
    int kernel_1, int kernel_2, int kernel_3, 
    int kernel_4, int kernel_5, int kernel_6, 
    int kernel_7, int kernel_8, int kernel_9,
    int opt_register_transpose)
{
	// # of Blocks for Each Kernel
	int	 num_blocks_kernel_1,		num_blocks_kernel_2;
	int  size_internal = size_p7;

	// Device Memory for Inputs and Output
	double *dev_t3;
	double *dev_t2_1, *dev_t2_2, *dev_t2_3, *dev_t2_4, *dev_t2_5, *dev_t2_6, *dev_t2_7, *dev_t2_8, *dev_t2_9;
	double *dev_v2_1, *dev_v2_2, *dev_v2_3, *dev_v2_4, *dev_v2_5, *dev_v2_6, *dev_v2_7, *dev_v2_8, *dev_v2_9;
    
#ifdef NWCHEM_TCE_CCSD_T
    dev_t3 = t3_d;
#else
    cudaMalloc((void**) &dev_t3, 	sizeof(double) * size_h3 * size_h2 * size_h1 * size_p6 * size_p5 * size_p4);
#endif

	cudaMalloc((void**) &dev_t2_1, 	sizeof(double) * size_h2 * size_h1 * size_p4 * size_p7);
	cudaMalloc((void**) &dev_v2_1, 	sizeof(double) * size_p5 * size_p6 * size_h3 * size_p7);
	cudaMalloc((void**) &dev_t2_2, 	sizeof(double) * size_h3 * size_h2 * size_p4 * size_p7);
	cudaMalloc((void**) &dev_v2_2, 	sizeof(double) * size_p5 * size_p6 * size_h1 * size_p7);
	cudaMalloc((void**) &dev_t2_3, 	sizeof(double) * size_h3 * size_h1 * size_p4 * size_p7);
	cudaMalloc((void**) &dev_v2_3, 	sizeof(double) * size_p5 * size_p6 * size_h2 * size_p7);
	cudaMalloc((void**) &dev_t2_4, 	sizeof(double) * size_h2 * size_h1 * size_p5 * size_p7);
	cudaMalloc((void**) &dev_v2_4, 	sizeof(double) * size_p4 * size_p6 * size_h3 * size_p7);
	cudaMalloc((void**) &dev_t2_5, 	sizeof(double) * size_h3 * size_h2 * size_p5 * size_p7);
	cudaMalloc((void**) &dev_v2_5, 	sizeof(double) * size_p4 * size_p6 * size_h1 * size_p7);
	cudaMalloc((void**) &dev_t2_6, 	sizeof(double) * size_h3 * size_h1 * size_p5 * size_p7);
	cudaMalloc((void**) &dev_v2_6, 	sizeof(double) * size_p4 * size_p6 * size_h2 * size_p7);

	cudaMalloc((void**) &dev_t2_7, 	sizeof(double) * size_h2 * size_h1 * size_p6 * size_p7);
	cudaMalloc((void**) &dev_v2_7, 	sizeof(double) * size_p4 * size_p5 * size_h3 * size_p7);
	cudaMalloc((void**) &dev_t2_8, 	sizeof(double) * size_h3 * size_h2 * size_p6 * size_p7);
	cudaMalloc((void**) &dev_v2_8, 	sizeof(double) * size_p4 * size_p5 * size_h1 * size_p7);
	cudaMalloc((void**) &dev_t2_9, 	sizeof(double) * size_h3 * size_h1 * size_p6 * size_p7);
	cudaMalloc((void**) &dev_v2_9, 	sizeof(double) * size_p4 * size_p5 * size_h2 * size_p7);

#ifndef NWCHEM_TCE_CCSD_T
    cudaMemcpy(dev_t3, 	 t3,   sizeof(double) * size_h3 * size_h2 * size_h1 * size_p6 * size_p5 * size_p4, cudaMemcpyHostToDevice);
#endif

	cudaMemcpy(dev_t2_1, t2_1, sizeof(double) * size_h2 * size_h1 * size_p4 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_1, v2_1, sizeof(double) * size_p5 * size_p6 * size_h3 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_2, t2_2, sizeof(double) * size_h3 * size_h2 * size_p4 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_2, v2_2, sizeof(double) * size_p5 * size_p6 * size_h1 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_3, t2_3, sizeof(double) * size_h3 * size_h1 * size_p4 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_3, v2_3, sizeof(double) * size_p5 * size_p6 * size_h2 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_4, t2_4, sizeof(double) * size_h2 * size_h1 * size_p5 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_4, v2_4, sizeof(double) * size_p4 * size_p6 * size_h3 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_5, t2_5, sizeof(double) * size_h3 * size_h2 * size_p5 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_5, v2_5, sizeof(double) * size_p4 * size_p6 * size_h1 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_6, t2_6, sizeof(double) * size_h3 * size_h1 * size_p5 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_6, v2_6, sizeof(double) * size_p4 * size_p6 * size_h2 * size_p7, cudaMemcpyHostToDevice);

	cudaMemcpy(dev_t2_7, t2_7, sizeof(double) * size_h2 * size_h1 * size_p6 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_7, v2_7, sizeof(double) * size_p4 * size_p5 * size_h3 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_8, t2_8, sizeof(double) * size_h3 * size_h2 * size_p6 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_8, v2_8, sizeof(double) * size_p4 * size_p5 * size_h1 * size_p7, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_9, t2_9, sizeof(double) * size_h3 * size_h1 * size_p6 * size_p7, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_v2_9, v2_9, sizeof(double) * size_p4 * size_p5 * size_h2 * size_p7, cudaMemcpyHostToDevice);
    
    num_blocks_kernel_1 =   CEIL(size_h3, SIZE_SLICE_1_H3) * CEIL(size_h2, SIZE_SLICE_1_H2) * CEIL(size_h1, SIZE_SLICE_1_H1) * 
                            CEIL(size_p6, SIZE_SLICE_1_P6) * CEIL(size_p5, SIZE_SLICE_1_P5) * CEIL(size_p4, SIZE_SLICE_1_P4);

    num_blocks_kernel_2 =   CEIL(size_h3, SIZE_SLICE_2_H3) * CEIL(size_h2, SIZE_SLICE_2_H2) * CEIL(size_h1, SIZE_SLICE_2_H1) * 
                            CEIL(size_p6, SIZE_SLICE_2_P6) * CEIL(size_p5, SIZE_SLICE_2_P5) * CEIL(size_p4, SIZE_SLICE_2_P4);

#ifndef NWCHEM_TCE_CCSD_T
	// (5) launch kernel(s)
	printf ("========================================= fusedKernels =============================================\n");
	printf ("		Grid Size  : %6d (1D)\n", num_blocks_kernel_1);
	printf ("		Block-size : %2d, %2d (2D)\n", SIZE_TB_1_X, SIZE_TB_1_Y);
	printf ("		A thread deals with (%d x %d) elements (basically)\n", 4 * 16, 4 * 16);
	printf ("====================================================================================================\n");
	printf ("========================================= fusedKernels =============================================\n");
	printf ("		Grid Size  : %6d (1D)\n", num_blocks_kernel_2);
	printf ("		Block-size : %2d, %2d (2D)\n", SIZE_TB_2_X, SIZE_TB_2_Y);
	printf ("		A thread deals with (%d x %d) elements (basically)\n", 4 * 16, 4 * 16);
    printf ("====================================================================================================\n");
#endif

	// Depends on # of Fused Kernel
	dim3 gridsize_1(num_blocks_kernel_1);
	dim3 blocksize_1(SIZE_TB_1_X, SIZE_TB_1_Y);

	dim3 gridsize_2(num_blocks_kernel_2);
	dim3 blocksize_2(SIZE_TB_2_X, SIZE_TB_2_Y);

	int	str_sd2_t3_h3 = 1;
	int str_sd2_t3_h2 = str_sd2_t3_h3 * size_h3;
	int str_sd2_t3_h1 = str_sd2_t3_h2 * size_h2;
	int str_sd2_t3_p6 = str_sd2_t3_h1 * size_h1;
	int str_sd2_t3_p5 = str_sd2_t3_p6 * size_p6;
	int str_sd2_t3_p4 = str_sd2_t3_p5 * size_p5;

	int str_reg_x_1 = str_sd2_t3_p5;	// STR_SD2_T3_P5
	int str_reg_y_1 = str_sd2_t3_p4;	// STR_SD2_T3_P4
	int str_reg_x_2 = str_sd2_t3_p5;	// STR_SD2_T3_P5
	int str_reg_y_2 = str_sd2_t3_p6;	// SDT_SD2_T3_P6

    if (kernel_1 || kernel_2 || kernel_3 || kernel_4 || kernel_5 || kernel_6)
    {
        if (kernel_7 || kernel_8 || kernel_9)
        {
            // 1,2,3,4,5,6  [ON]
            // 7,8,9        [ON]
            if (size_h3 % SIZE_SLICE_1_H3 == 0 && size_h2 % SIZE_SLICE_1_H2 == 0 && size_h1 % SIZE_SLICE_1_H1 == 0 && 
                size_p6 % SIZE_SLICE_1_P6 == 0 && size_p5 % SIZE_SLICE_1_P5 == 0 && size_p4 % SIZE_SLICE_1_P4 == 0)
            {
                if (size_p7 % SIZE_INT_UNIT == 0)
                {
                    kernel_ccsdT_sd2_fully_fused_full_full<<<gridsize_1, blocksize_1>>>(dev_t3, 
                    dev_t2_1, dev_t2_2, dev_t2_3, dev_t2_4, dev_t2_5, dev_t2_6, 
                    dev_v2_1, dev_v2_2, dev_v2_3, dev_v2_4, dev_v2_5, dev_v2_6, 
                    dev_t2_7, dev_t2_8, dev_t2_9, 
                    dev_v2_7, dev_v2_8, dev_v2_9, 
                    size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7,
                    CEIL(size_h3, SIZE_SLICE_1_H3), CEIL(size_h2, SIZE_SLICE_1_H2), CEIL(size_h1, SIZE_SLICE_1_H1), 
                    CEIL(size_p6, SIZE_SLICE_1_P6), CEIL(size_p5, SIZE_SLICE_1_P5), CEIL(size_p4, SIZE_SLICE_1_P4), 
                    kernel_1, kernel_2, kernel_3,
                    kernel_4, kernel_5, kernel_6,
                    kernel_7, kernel_8, kernel_9,
                    str_reg_x_1, str_reg_y_1,
                    size_internal);	
                }
                else
                {
                    kernel_ccsdT_sd2_fully_fused_partial_full<<<gridsize_1, blocksize_1>>>(dev_t3, 
                    dev_t2_1, dev_t2_2, dev_t2_3, dev_t2_4, dev_t2_5, dev_t2_6, 
                    dev_v2_1, dev_v2_2, dev_v2_3, dev_v2_4, dev_v2_5, dev_v2_6, 
                    dev_t2_7, dev_t2_8, dev_t2_9, 
                    dev_v2_7, dev_v2_8, dev_v2_9, 
                    size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7,
                    CEIL(size_h3, SIZE_SLICE_1_H3), CEIL(size_h2, SIZE_SLICE_1_H2), CEIL(size_h1, SIZE_SLICE_1_H1), 
                    CEIL(size_p6, SIZE_SLICE_1_P6), CEIL(size_p5, SIZE_SLICE_1_P5), CEIL(size_p4, SIZE_SLICE_1_P4), 
                    kernel_1, kernel_2, kernel_3,
                    kernel_4, kernel_5, kernel_6,
                    kernel_7, kernel_8, kernel_9,
                    str_reg_x_1, str_reg_y_1,
                    size_internal);	
                }
            }
            else
            {
                kernel_ccsdT_sd2_fully_fused_partial_partial<<<gridsize_1, blocksize_1>>>(dev_t3, 
                dev_t2_1, dev_t2_2, dev_t2_3, dev_t2_4, dev_t2_5, dev_t2_6, 
                dev_v2_1, dev_v2_2, dev_v2_3, dev_v2_4, dev_v2_5, dev_v2_6, 
                dev_t2_7, dev_t2_8, dev_t2_9, 
                dev_v2_7, dev_v2_8, dev_v2_9, 
                size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7,
                CEIL(size_h3, SIZE_SLICE_1_H3), CEIL(size_h2, SIZE_SLICE_1_H2), CEIL(size_h1, SIZE_SLICE_1_H1), 
                CEIL(size_p6, SIZE_SLICE_1_P6), CEIL(size_p5, SIZE_SLICE_1_P5), CEIL(size_p4, SIZE_SLICE_1_P4), 
                kernel_1, kernel_2, kernel_3,
                kernel_4, kernel_5, kernel_6,
                kernel_7, kernel_8, kernel_9,
                str_reg_x_1, str_reg_y_1,
                size_internal);	
            }
        }
        else
        {
            // 1,2,3,4,5,6  [ON]
            // 7,8,9        [OFF]
            if (size_h3 % SIZE_SLICE_1_H3 == 0 && size_h2 % SIZE_SLICE_1_H2 == 0 && size_h1 % SIZE_SLICE_1_H1 == 0 && 
                size_p6 % SIZE_SLICE_1_P6 == 0 && size_p5 % SIZE_SLICE_1_P5 == 0 && size_p4 % SIZE_SLICE_1_P4 == 0)
            {
                if (size_p7 % SIZE_INT_UNIT == 0)
                {
                    kernel_ccsdT_sd2_123456_full_full<<<gridsize_1, blocksize_1>>>(dev_t3, 
                    dev_t2_1, dev_t2_2, dev_t2_3, dev_t2_4, dev_t2_5, dev_t2_6, 
                    dev_v2_1, dev_v2_2, dev_v2_3, dev_v2_4, dev_v2_5, dev_v2_6, 
                    size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7,
                    CEIL(size_h3, SIZE_SLICE_1_H3), CEIL(size_h2, SIZE_SLICE_1_H2), CEIL(size_h1, SIZE_SLICE_1_H1), 
                    CEIL(size_p6, SIZE_SLICE_1_P6), CEIL(size_p5, SIZE_SLICE_1_P5), CEIL(size_p4, SIZE_SLICE_1_P4), 
                    kernel_1, kernel_2, kernel_3,
                    kernel_4, kernel_5, kernel_6,
                    kernel_7, kernel_8, kernel_9,
                    str_reg_x_1, str_reg_y_1,
                    size_internal);
                }
                else
                {
                    kernel_ccsdT_sd2_123456_partial_full<<<gridsize_1, blocksize_1>>>(dev_t3, 
                    dev_t2_1, dev_t2_2, dev_t2_3, dev_t2_4, dev_t2_5, dev_t2_6, 
                    dev_v2_1, dev_v2_2, dev_v2_3, dev_v2_4, dev_v2_5, dev_v2_6, 
                    size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7,
                    CEIL(size_h3, SIZE_SLICE_1_H3), CEIL(size_h2, SIZE_SLICE_1_H2), CEIL(size_h1, SIZE_SLICE_1_H1), 
                    CEIL(size_p6, SIZE_SLICE_1_P6), CEIL(size_p5, SIZE_SLICE_1_P5), CEIL(size_p4, SIZE_SLICE_1_P4), 
                    kernel_1, kernel_2, kernel_3,
                    kernel_4, kernel_5, kernel_6,
                    kernel_7, kernel_8, kernel_9,
                    str_reg_x_1, str_reg_y_1,
                    size_internal);
                }
            }
            else
            {
                kernel_ccsdT_sd2_123456_partial_partial<<<gridsize_1, blocksize_1>>>(dev_t3, 
                dev_t2_1, dev_t2_2, dev_t2_3, dev_t2_4, dev_t2_5, dev_t2_6, 
                dev_v2_1, dev_v2_2, dev_v2_3, dev_v2_4, dev_v2_5, dev_v2_6, 
                size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7,
                CEIL(size_h3, SIZE_SLICE_1_H3), CEIL(size_h2, SIZE_SLICE_1_H2), CEIL(size_h1, SIZE_SLICE_1_H1), 
                CEIL(size_p6, SIZE_SLICE_1_P6), CEIL(size_p5, SIZE_SLICE_1_P5), CEIL(size_p4, SIZE_SLICE_1_P4), 
                kernel_1, kernel_2, kernel_3,
                kernel_4, kernel_5, kernel_6,
                kernel_7, kernel_8, kernel_9,
                str_reg_x_1, str_reg_y_1,
                size_internal);
            }
        }
    }
    else
    {
        if (kernel_7 || kernel_8 || kernel_9)
        {
            // 1,2,3,4,5,6  [OFF]
            // 7,8,9        [ON]
            if (size_h3 % SIZE_SLICE_1_H3 == 0 && size_h2 % SIZE_SLICE_1_H2 == 0 && size_h1 % SIZE_SLICE_1_H1 == 0 && 
                size_p6 % SIZE_SLICE_1_P6 == 0 && size_p5 % SIZE_SLICE_1_P5 == 0 && size_p4 % SIZE_SLICE_1_P4 == 0)
            {
                if (size_p7 % SIZE_INT_UNIT == 0)
                {
                    kernel_ccsdT_sd2_789_full_full<<<gridsize_2, blocksize_2>>>(dev_t3, 
                    dev_t2_7, dev_t2_8, dev_t2_9, 
                    dev_v2_7, dev_v2_8, dev_v2_9, 
                    size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7,
                    CEIL(size_h3, SIZE_SLICE_2_H3), CEIL(size_h2, SIZE_SLICE_2_H2), CEIL(size_h1, SIZE_SLICE_2_H1), 
                    CEIL(size_p6, SIZE_SLICE_2_P6), CEIL(size_p5, SIZE_SLICE_2_P5), CEIL(size_p4, SIZE_SLICE_2_P4), 
                    kernel_1, kernel_2, kernel_3,
                    kernel_4, kernel_5, kernel_6,
                    kernel_7, kernel_8, kernel_9,
                    str_reg_x_2, str_reg_y_2,
                    size_internal);
                }
                else
                {
                    kernel_ccsdT_sd2_789_partial_full<<<gridsize_2, blocksize_2>>>(dev_t3, 
                    dev_t2_7, dev_t2_8, dev_t2_9, 
                    dev_v2_7, dev_v2_8, dev_v2_9, 
                    size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7,
                    CEIL(size_h3, SIZE_SLICE_2_H3), CEIL(size_h2, SIZE_SLICE_2_H2), CEIL(size_h1, SIZE_SLICE_2_H1), 
                    CEIL(size_p6, SIZE_SLICE_2_P6), CEIL(size_p5, SIZE_SLICE_2_P5), CEIL(size_p4, SIZE_SLICE_2_P4), 
                    kernel_1, kernel_2, kernel_3,
                    kernel_4, kernel_5, kernel_6,
                    kernel_7, kernel_8, kernel_9,
                    str_reg_x_2, str_reg_y_2,
                    size_internal);
                }
            }
            else
            {
                kernel_ccsdT_sd2_789_partial_partial<<<gridsize_2, blocksize_2>>>(dev_t3, 
                dev_t2_7, dev_t2_8, dev_t2_9, 
                dev_v2_7, dev_v2_8, dev_v2_9, 
                size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7,
                CEIL(size_h3, SIZE_SLICE_2_H3), CEIL(size_h2, SIZE_SLICE_2_H2), CEIL(size_h1, SIZE_SLICE_2_H1), 
                CEIL(size_p6, SIZE_SLICE_2_P6), CEIL(size_p5, SIZE_SLICE_2_P5), CEIL(size_p4, SIZE_SLICE_2_P4), 
                kernel_1, kernel_2, kernel_3,
                kernel_4, kernel_5, kernel_6,
                kernel_7, kernel_8, kernel_9,
                str_reg_x_2, str_reg_y_2,
                size_internal);
            }
        }
        else
        {
            // 1,2,3,4,5,6  [OFF]
            // 7,8,9        [OFF]
            printf (">>> kernel_1,2,3,4,5,6 (OFF) && kernel_7,8,9 (OFF) <<< %d,%d,%d,%d,%d,%d,%d,%d,%d\n", kernel_1, kernel_2, kernel_3, kernel_4, kernel_5, kernel_6, kernel_7, kernel_8, kernel_9);
        }
    }

#ifndef NWCHEM_TCE_CCSD_T
	// Copy the Result from Device to Host
	cudaMemcpy(t3, dev_t3, sizeof(double) * (size_h3 * size_h2 * size_h1 * size_p6 * size_p5 * size_p4), cudaMemcpyDeviceToHost);

	// cudaFree()
    cudaFree(dev_t3);
#endif
	cudaFree(dev_t2_1);	cudaFree(dev_t2_2);	cudaFree(dev_t2_3);	cudaFree(dev_t2_4);	cudaFree(dev_t2_5);	cudaFree(dev_t2_6);	cudaFree(dev_t2_7);	cudaFree(dev_t2_8);	cudaFree(dev_t2_9);	
	cudaFree(dev_v2_1);	cudaFree(dev_v2_2);	cudaFree(dev_v2_3);	cudaFree(dev_v2_4);	cudaFree(dev_v2_5);	cudaFree(dev_v2_6);	cudaFree(dev_v2_7);	cudaFree(dev_v2_8);	cudaFree(dev_v2_9);
}

extern "C" 
void sd_t_d2_all_cuda_(int size_h3, int size_h2, int size_h1, int size_p6, int size_p5, int size_p4, int size_p7, 
                        double* t3, 
                        double* t2_1, double* v2_1,
                        double* t2_2, double* v2_2,
                        double* t2_3, double* v2_3,
                        double* t2_4, double* v2_4,
                        double* t2_5, double* v2_5,
                        double* t2_6, double* v2_6,
                        double* t2_7, double* v2_7,
                        double* t2_8, double* v2_8,
                        double* t2_9, double* v2_9, 
                        int kernel_1, int kernel_2, int kernel_3, 
                        int kernel_4, int kernel_5, int kernel_6, 
                        int kernel_7, int kernel_8, int kernel_9,
                        int opt_register_transpose)
{

    if (kernel_1 && kernel_2 && kernel_3 && kernel_4 && kernel_5 && kernel_6 && kernel_7 && kernel_8 && kernel_8 && kernel_9)
    {
        sd_t_d2_all_cuda(size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7, t3, t2_1, v2_1,
        t2_2, v2_2,
        t2_3, v2_3,
        t2_4, v2_4,
        t2_5, v2_5,
        t2_6, v2_6,
        t2_7, v2_7,
        t2_8, v2_8,
        t2_9, v2_9, 
        kernel_1, kernel_2, kernel_3,
        kernel_4, kernel_5, kernel_6,
        kernel_7, kernel_8, kernel_9,
        opt_register_transpose);
    }
    else
    {
        if (kernel_1)
        sd_t_d2_1_cogent(size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7, t3, t2_1, v2_1, kernel_1, opt_register_transpose);

        if (kernel_2)
        sd_t_d2_2_cogent(size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7, t3, t2_2, v2_2, kernel_2, opt_register_transpose);

        if (kernel_3)
        sd_t_d2_3_cogent(size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7, t3, t2_3, v2_3, kernel_3, opt_register_transpose);

        if (kernel_4)
        sd_t_d2_4_cogent(size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7, t3, t2_4, v2_4, kernel_4, opt_register_transpose);

        if (kernel_5)
        sd_t_d2_5_cogent(size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7, t3, t2_5, v2_5, kernel_5, opt_register_transpose);

        if (kernel_6)
        sd_t_d2_6_cogent(size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7, t3, t2_6, v2_6, kernel_6, opt_register_transpose);

        if (kernel_7)
        sd_t_d2_7_cogent(size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7, t3, t2_7, v2_7, kernel_7, opt_register_transpose);

        if (kernel_8)
        sd_t_d2_8_cogent(size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7, t3, t2_8, v2_8, kernel_8, opt_register_transpose);

        if (kernel_9)
        sd_t_d2_9_cogent(size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7, t3, t2_9, v2_9, kernel_9, opt_register_transpose);
    }
}

//
extern "C"
void sd_t_d2_all_cuda__(Integer* p_size_h3, Integer* p_size_h2, Integer* p_size_h1, Integer* p_size_p6, Integer* p_size_p5, Integer* p_size_p4, Integer* p_size_p7,
                        double* t3, 
                        double* t2_all, Integer* p_size_t2_all,
                        double* v2_all, Integer* p_size_v2_all,
                        Integer* p_kernel_1, Integer* p_kernel_2, Integer* p_kernel_3, 
                        Integer* p_kernel_4, Integer* p_kernel_5, Integer* p_kernel_6, 
                        Integer* p_kernel_7, Integer* p_kernel_8, Integer* p_kernel_9,
			Integer* p_opt_register_transpose)
{
    int size_h3 = *p_size_h3;
    int size_h2 = *p_size_h2;
    int size_h1 = *p_size_h1;
    int size_p6 = *p_size_p6;
    int size_p5 = *p_size_p5;
    int size_p4 = *p_size_p4;
    int size_p7 = *p_size_p7;

    int kernel_1 = *p_kernel_1;
    int kernel_2 = *p_kernel_2;
    int kernel_3 = *p_kernel_3;
    int kernel_4 = *p_kernel_4;
    int kernel_5 = *p_kernel_5;
    int kernel_6 = *p_kernel_6;
    int kernel_7 = *p_kernel_7;
    int kernel_8 = *p_kernel_8;
    int kernel_9 = *p_kernel_9;

    int opt_register_transpose = *p_opt_register_transpose;

    int size_t2_all = *p_size_t2_all;
    int size_v2_all = *p_size_v2_all;

    unsigned int size_t2_each = size_t2_all / 9;
    unsigned int size_v2_each = size_v2_all / 9;

    double* t2_1 = t2_all;
    double* t2_2 = t2_all + (size_t2_each * 1);
    double* t2_3 = t2_all + (size_t2_each * 2);
    double* t2_4 = t2_all + (size_t2_each * 3);
    double* t2_5 = t2_all + (size_t2_each * 4);
    double* t2_6 = t2_all + (size_t2_each * 5);
    double* t2_7 = t2_all + (size_t2_each * 6);
    double* t2_8 = t2_all + (size_t2_each * 7);
    double* t2_9 = t2_all + (size_t2_each * 8);

    double* v2_1 = v2_all;
    double* v2_2 = v2_all + (size_v2_each * 1);
    double* v2_3 = v2_all + (size_v2_each * 2);
    double* v2_4 = v2_all + (size_v2_each * 3);
    double* v2_5 = v2_all + (size_v2_each * 4);
    double* v2_6 = v2_all + (size_v2_each * 5);
    double* v2_7 = v2_all + (size_v2_each * 6);
    double* v2_8 = v2_all + (size_v2_each * 7);
    double* v2_9 = v2_all + (size_v2_each * 8);

    //
    sd_t_d2_all_cuda_(size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7, 
        t3, 
        t2_1, v2_1, t2_2, v2_2, t2_3, v2_3,
        t2_4, v2_4, t2_5, v2_5, t2_6, v2_6,
        t2_7, v2_7, t2_8, v2_8, t2_9, v2_9, 
        kernel_1, kernel_2, kernel_3,
        kernel_4, kernel_5, kernel_6,
        kernel_7, kernel_8, kernel_9,
	    opt_register_transpose);
}
