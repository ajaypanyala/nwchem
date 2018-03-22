#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <locale.h>
#include <algorithm>
#include "header.h"
using namespace std;

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
#define CEIL(a, b) 		(((a) + (b) - 1) / (b))

void Pre_PreComputedArrays_1(int*& h_t3_output_base_addr, 
							int*& h_t2_1_ext_addr, int*& h_v2_1_ext_addr,
							int*& h_t2_2_ext_addr, int*& h_v2_2_ext_addr,
							int*& h_t2_3_ext_addr, int*& h_v2_3_ext_addr,
							int*& h_t2_4_ext_addr, int*& h_v2_4_ext_addr,
							int*& h_t2_5_ext_addr, int*& h_v2_5_ext_addr,
							int*& h_t2_6_ext_addr, int*& h_v2_6_ext_addr,
							int*& h_t3_output_ext_offset,
							int*& h_t2_1_ext_offset, int*& h_v2_1_ext_offset,
							int*& h_t2_2_ext_offset, int*& h_v2_2_ext_offset,
							int*& h_t2_3_ext_offset, int*& h_v2_3_ext_offset,
							int*& h_t2_4_ext_offset, int*& h_v2_4_ext_offset,
							int*& h_t2_5_ext_offset, int*& h_v2_5_ext_offset,
							int*& h_t2_6_ext_offset, int*& h_v2_6_ext_offset,
							int* h_t3_block_index,
							int n_blocks_1, 
							int size_idx_h3, 		int size_idx_h2,		int size_idx_h1,		int size_idx_p6, 		int size_idx_p5, 	int size_idx_p4, 		int size_idx_p7)
{
	//
	int	str_sd2_t3_h3 = 1;
	int str_sd2_t3_h2 = str_sd2_t3_h3 * size_idx_h3;
	int str_sd2_t3_h1 = str_sd2_t3_h2 * size_idx_h2;
	int str_sd2_t3_p6 = str_sd2_t3_h1 * size_idx_h1;
	int str_sd2_t3_p5 = str_sd2_t3_p6 * size_idx_p6;
	int str_sd2_t3_p4 = str_sd2_t3_p5 * size_idx_p5;

	//
	h_t3_output_base_addr = (int*)malloc(sizeof(int) * n_blocks_1);
	h_t2_1_ext_addr = (int*)malloc(sizeof(int) * SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2 * n_blocks_1);
	h_v2_1_ext_addr = (int*)malloc(sizeof(int) * SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5 * n_blocks_1);

	h_t2_2_ext_addr = (int*)malloc(sizeof(int) * SIZE_SLICE_1_P4 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3 * n_blocks_1);
	h_v2_2_ext_addr = (int*)malloc(sizeof(int) * SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5 * n_blocks_1);

	h_t2_3_ext_addr = (int*)malloc(sizeof(int) * SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3 * n_blocks_1);
	h_v2_3_ext_addr = (int*)malloc(sizeof(int) * SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5 * n_blocks_1);

	h_t2_4_ext_addr = (int*)malloc(sizeof(int) * SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2 * n_blocks_1);
	h_v2_4_ext_addr = (int*)malloc(sizeof(int) * SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4 * n_blocks_1);

	h_t2_5_ext_addr = (int*)malloc(sizeof(int) * SIZE_SLICE_1_P5 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3 * n_blocks_1);
	h_v2_5_ext_addr = (int*)malloc(sizeof(int) * SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4 * n_blocks_1);

	h_t2_6_ext_addr = (int*)malloc(sizeof(int) * SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3 * n_blocks_1);
	h_v2_6_ext_addr = (int*)malloc(sizeof(int) * SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4 * n_blocks_1);

	// tc_gen_code_pre_IndirectArray_wo_TB()
	// For Each Thread Block
	for (int i = 0; i < n_blocks_1; i++)
	{
		int blk_idx_h3 = h_t3_block_index[i * NUM_INDEX + 0];
		int blk_idx_h2 = h_t3_block_index[i * NUM_INDEX + 1];
		int blk_idx_h1 = h_t3_block_index[i * NUM_INDEX + 2];
		int blk_idx_p6 = h_t3_block_index[i * NUM_INDEX + 3];
		int blk_idx_p5 = h_t3_block_index[i * NUM_INDEX + 4];
		int blk_idx_p4 = h_t3_block_index[i * NUM_INDEX + 5];

		// calculating t3
		h_t3_output_base_addr[i] = 	blk_idx_h3 * SIZE_SLICE_1_H3 * str_sd2_t3_h3 + blk_idx_h2 * SIZE_SLICE_1_H2 * str_sd2_t3_h2 + blk_idx_h1 * SIZE_SLICE_1_H1 * str_sd2_t3_h1 + 
									blk_idx_p6 * SIZE_SLICE_1_P6 * str_sd2_t3_p6 + blk_idx_p5 * SIZE_SLICE_1_P5 * str_sd2_t3_p5 + blk_idx_p4 * SIZE_SLICE_1_P4 * str_sd2_t3_p4;

		// calculating actual address for t2_1 and v2_1
		for (int idx_p4 = 0; idx_p4 < SIZE_SLICE_1_P4; idx_p4++)
		for (int idx_h2 = 0; idx_h2 < SIZE_SLICE_1_H2; idx_h2++)
		for (int idx_h1 = 0; idx_h1 < SIZE_SLICE_1_H1; idx_h1++)
		{
			int l_offset = idx_h1 + (idx_h2 + (idx_p4) * SIZE_SLICE_1_H2) * SIZE_SLICE_1_H1;
			h_t2_1_ext_addr[l_offset + i * SIZE_SLICE_1_P4 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H1] = (blk_idx_p4 * SIZE_SLICE_1_P4 + idx_p4 + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_h1 + 
																								(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h2) * size_idx_h1) * size_idx_p4) * size_idx_p7;
		}

		for (int idx_p5 = 0; idx_p5 < SIZE_SLICE_1_P5; idx_p5++)
		for (int idx_p6 = 0; idx_p6 < SIZE_SLICE_1_P6; idx_p6++)
		for (int idx_h3 = 0; idx_h3 < SIZE_SLICE_1_H3; idx_h3++)
		{
			int l_offset = idx_h3 + (idx_p6 + (idx_p5) * SIZE_SLICE_1_P6) * SIZE_SLICE_1_H3;
			h_v2_1_ext_addr[l_offset + i * SIZE_SLICE_1_P5 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H3] = (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h3 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_p6 + 
																								(blk_idx_p5 * SIZE_SLICE_1_P5 + idx_p5) * size_idx_p6) * size_idx_h3) * size_idx_p7;
		}
		// calculating actual address for t2_2 and v2_2
		for (int idx_p4 = 0; idx_p4 < SIZE_SLICE_1_P4; idx_p4++)
		for (int idx_h3 = 0; idx_h3 < SIZE_SLICE_1_H3; idx_h3++)
		for (int idx_h2 = 0; idx_h2 < SIZE_SLICE_1_H2; idx_h2++)
		{
			int l_offset = idx_h2 + (idx_h3 + (idx_p4) * SIZE_SLICE_1_H3) * SIZE_SLICE_1_H2;
			h_t2_2_ext_addr[l_offset + i * SIZE_SLICE_1_P4 * SIZE_SLICE_1_H3 * SIZE_SLICE_1_H2] = (blk_idx_p4 * SIZE_SLICE_1_P4 + idx_p4 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h2 + 
																								(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h3) * size_idx_h2) * size_idx_p4) * size_idx_p7;
		}

		for (int idx_p5 = 0; idx_p5 < SIZE_SLICE_1_P5; idx_p5++)
		for (int idx_p6 = 0; idx_p6 < SIZE_SLICE_1_P6; idx_p6++)
		for (int idx_h1 = 0; idx_h1 < SIZE_SLICE_1_H1; idx_h1++)
		{
			int l_offset = idx_h1 + (idx_p6 + (idx_p5) * SIZE_SLICE_1_P6) * SIZE_SLICE_1_H1;
			h_v2_2_ext_addr[l_offset + i * SIZE_SLICE_1_P5 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H1] = (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_h1 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_p6 + 
																								(blk_idx_p5 * SIZE_SLICE_1_P5 + idx_p5) * size_idx_p6) * size_idx_h1) * size_idx_p7;
		}
		// calculating actual address for t2_3 and v2_3
		for (int idx_p4 = 0; idx_p4 < SIZE_SLICE_1_P4; idx_p4++)
		for (int idx_h3 = 0; idx_h3 < SIZE_SLICE_1_H3; idx_h3++)
		for (int idx_h1 = 0; idx_h1 < SIZE_SLICE_1_H1; idx_h1++)
		{
			int l_offset = idx_h1 + (idx_h3 + (idx_p4) * SIZE_SLICE_1_H3) * SIZE_SLICE_1_H1;
			h_t2_3_ext_addr[l_offset + i * SIZE_SLICE_1_P4 * SIZE_SLICE_1_H3 * SIZE_SLICE_1_H1] = (blk_idx_p4 * SIZE_SLICE_1_P4 + idx_p4 + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_h1 + 
																								(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h3) * size_idx_h1) * size_idx_p4) * size_idx_p7;
		}

		for (int idx_p5 = 0; idx_p5 < SIZE_SLICE_1_P5; idx_p5++)
		for (int idx_p6 = 0; idx_p6 < SIZE_SLICE_1_P6; idx_p6++)
		for (int idx_h2 = 0; idx_h2 < SIZE_SLICE_1_H2; idx_h2++)
		{
			int l_offset = idx_h2 + (idx_p6 + (idx_p5) * SIZE_SLICE_1_P6) * SIZE_SLICE_1_H2;
			h_v2_3_ext_addr[l_offset + i * SIZE_SLICE_1_P5 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H2] = (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h2 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_p6 + 
																								(blk_idx_p5 * SIZE_SLICE_1_P5 + idx_p5) * size_idx_p6) * size_idx_h2) * size_idx_p7;
		}
		// calculating actual address for t2_4 and v2_4
		for (int idx_p5 = 0; idx_p5 < SIZE_SLICE_1_P5; idx_p5++)
		for (int idx_h2 = 0; idx_h2 < SIZE_SLICE_1_H2; idx_h2++)
		for (int idx_h1 = 0; idx_h1 < SIZE_SLICE_1_H1; idx_h1++)
		{
			int l_offset = idx_h1 + (idx_h2 + (idx_p5) * SIZE_SLICE_1_H2) * SIZE_SLICE_1_H1;
			h_t2_4_ext_addr[l_offset + i * SIZE_SLICE_1_P5 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H1] = (blk_idx_p5 * SIZE_SLICE_1_P5 + idx_p5 + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_h1 + 
																								(blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h2) * size_idx_h1) * size_idx_p5) * size_idx_p7;
		}

		for (int idx_p4 = 0; idx_p4 < SIZE_SLICE_1_P4; idx_p4++)
		for (int idx_p6 = 0; idx_p6 < SIZE_SLICE_1_P6; idx_p6++)
		for (int idx_h3 = 0; idx_h3 < SIZE_SLICE_1_H3; idx_h3++)
		{
			int l_offset = idx_h3 + (idx_p6 + (idx_p4) * SIZE_SLICE_1_P6) * SIZE_SLICE_1_H3;
			h_v2_4_ext_addr[l_offset + i * SIZE_SLICE_1_P4 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H3] = (blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h3 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_p6 + 
																								(blk_idx_p4 * SIZE_SLICE_1_P4 + idx_p4) * size_idx_p6) * size_idx_h3) * size_idx_p7;
		}
		// calculating actual address for t2_5 and v2_5
		for (int idx_p5 = 0; idx_p5 < SIZE_SLICE_1_P5; idx_p5++)
		for (int idx_h3 = 0; idx_h3 < SIZE_SLICE_1_H3; idx_h3++)
		for (int idx_h2 = 0; idx_h2 < SIZE_SLICE_1_H2; idx_h2++)
		{
			int l_offset = idx_h2 + (idx_h3 + (idx_p5) * SIZE_SLICE_1_H3) * SIZE_SLICE_1_H2;
			h_t2_5_ext_addr[l_offset + i * SIZE_SLICE_1_P5 * SIZE_SLICE_1_H3 * SIZE_SLICE_1_H2] = (blk_idx_p5 * SIZE_SLICE_1_P5 + idx_p5 + (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h2 + 
																								(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h3) * size_idx_h2) * size_idx_p5) * size_idx_p7;
		}

		for (int idx_p4 = 0; idx_p4 < SIZE_SLICE_1_P4; idx_p4++)
		for (int idx_p6 = 0; idx_p6 < SIZE_SLICE_1_P6; idx_p6++)
		for (int idx_h1 = 0; idx_h1 < SIZE_SLICE_1_H1; idx_h1++)
		{
			int l_offset = idx_h1 + (idx_p6 + (idx_p4) * SIZE_SLICE_1_P6) * SIZE_SLICE_1_H1;
			h_v2_5_ext_addr[l_offset + i * SIZE_SLICE_1_P4 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H1] = (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_h1 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_p6 + 
																								(blk_idx_p4 * SIZE_SLICE_1_P4 + idx_p4) * size_idx_p6) * size_idx_h1) * size_idx_p7;
		}
		// calculating actual address for t2_6 and v2_6
		for (int idx_p5 = 0; idx_p5 < SIZE_SLICE_1_P5; idx_p5++)
		for (int idx_h3 = 0; idx_h3 < SIZE_SLICE_1_H3; idx_h3++)
		for (int idx_h1 = 0; idx_h1 < SIZE_SLICE_1_H1; idx_h1++)
		{
			int l_offset = idx_h1 + (idx_h3 + (idx_p5) * SIZE_SLICE_1_H3) * SIZE_SLICE_1_H1;
			h_t2_6_ext_addr[l_offset + i * SIZE_SLICE_1_P5 * SIZE_SLICE_1_H3 * SIZE_SLICE_1_H1] = (blk_idx_p5 * SIZE_SLICE_1_P5 + idx_p5 + (blk_idx_h1 * SIZE_SLICE_1_H1 + idx_h1 + 
																								(blk_idx_h3 * SIZE_SLICE_1_H3 + idx_h3) * size_idx_h1) * size_idx_p5) * size_idx_p7;
		}

		for (int idx_p4 = 0; idx_p4 < SIZE_SLICE_1_P4; idx_p4++)
		for (int idx_p6 = 0; idx_p6 < SIZE_SLICE_1_P6; idx_p6++)
		for (int idx_h2 = 0; idx_h2 < SIZE_SLICE_1_H2; idx_h2++)
		{
			int l_offset = idx_h2 + (idx_p6 + (idx_p4) * SIZE_SLICE_1_P6) * SIZE_SLICE_1_H2;
			h_v2_6_ext_addr[l_offset + i * SIZE_SLICE_1_P4 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H2] = (blk_idx_h2 * SIZE_SLICE_1_H2 + idx_h2 + (blk_idx_p6 * SIZE_SLICE_1_P6 + idx_p6 + 
																								(blk_idx_p4 * SIZE_SLICE_1_P4 + idx_p4) * size_idx_p6) * size_idx_h2) * size_idx_p7;
		}
	}

	// tc_gen_code_pre_IndirectArray_w_TB()
	// Within a Thread Block
	h_t3_output_ext_offset 	= (int*)malloc(sizeof(int) * (SIZE_TB_1_X * SIZE_TB_1_Y));
	h_t2_1_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_1_X * SIZE_TB_1_Y));
	h_v2_1_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_1_X * SIZE_TB_1_Y));
	h_t2_2_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_1_X * SIZE_TB_1_Y));
	h_v2_2_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_1_X * SIZE_TB_1_Y));
	h_t2_3_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_1_X * SIZE_TB_1_Y));
	h_v2_3_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_1_X * SIZE_TB_1_Y));
	h_t2_4_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_1_X * SIZE_TB_1_Y));
	h_v2_4_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_1_X * SIZE_TB_1_Y));
	h_t2_5_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_1_X * SIZE_TB_1_Y));
	h_v2_5_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_1_X * SIZE_TB_1_Y));
	h_t2_6_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_1_X * SIZE_TB_1_Y));
	h_v2_6_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_1_X * SIZE_TB_1_Y));

	for (int idx_h1 = 0; idx_h1 < SIZE_SLICE_1_H1; idx_h1++)
	for (int idx_p6 = 0; idx_p6 < SIZE_SLICE_1_P6; idx_p6++)
	for (int idx_h2 = 0; idx_h2 < SIZE_SLICE_1_H2; idx_h2++)
	for (int idx_h3 = 0; idx_h3 < SIZE_SLICE_1_H3; idx_h3++)
	{
		int t3_index_1 = idx_h3 + (idx_h2 + (idx_p6 + (idx_h1) * SIZE_SLICE_1_P6) * SIZE_SLICE_1_H2) * SIZE_SLICE_1_H3;

		h_t3_output_ext_offset[t3_index_1] = idx_h3 + (idx_h2 + (idx_h1 + (idx_p6) * size_idx_h1) * size_idx_h2) * size_idx_h3;
		h_t2_1_ext_offset[t3_index_1] = idx_h1 + (idx_h2) * SIZE_SLICE_1_H1;
		h_v2_1_ext_offset[t3_index_1] = idx_h3 + (idx_p6) * SIZE_SLICE_1_H3;
		h_t2_2_ext_offset[t3_index_1] = idx_h2 + (idx_h3) * SIZE_SLICE_1_H2;
		h_v2_2_ext_offset[t3_index_1] = idx_h1 + (idx_p6) * SIZE_SLICE_1_H1;
		h_t2_3_ext_offset[t3_index_1] = idx_h1 + (idx_h3) * SIZE_SLICE_1_H1;
		h_v2_3_ext_offset[t3_index_1] = idx_h2 + (idx_p6) * SIZE_SLICE_1_H2;
		h_t2_4_ext_offset[t3_index_1] = idx_h1 + (idx_h2) * SIZE_SLICE_1_H1;
		h_v2_4_ext_offset[t3_index_1] = idx_h3 + (idx_p6) * SIZE_SLICE_1_H3;
		h_t2_5_ext_offset[t3_index_1] = idx_h2 + (idx_h3) * SIZE_SLICE_1_H2;
		h_v2_5_ext_offset[t3_index_1] = idx_h1 + (idx_p6) * SIZE_SLICE_1_H1;
		h_t2_6_ext_offset[t3_index_1] = idx_h1 + (idx_h3) * SIZE_SLICE_1_H1;
		h_v2_6_ext_offset[t3_index_1] = idx_h2 + (idx_p6) * SIZE_SLICE_1_H2;
	}
	// Do not need Indirect Arrays for Internal Indices
}

void Pre_PreComputedArrays_2(int*& h_t3_output_base_addr, 
							int*& h_t2_7_ext_addr, int*& h_v2_7_ext_addr,
							int*& h_t2_8_ext_addr, int*& h_v2_8_ext_addr,
							int*& h_t2_9_ext_addr, int*& h_v2_9_ext_addr,
							int*& h_t3_output_ext_offset,
							int*& h_t2_7_ext_offset, int*& h_v2_7_ext_offset,
							int*& h_t2_8_ext_offset, int*& h_v2_8_ext_offset,
							int*& h_t2_9_ext_offset, int*& h_v2_9_ext_offset,
							int* h_t3_block_index,
							int n_blocks_2, 
							int size_idx_h3, 		int size_idx_h2,		int size_idx_h1,		int size_idx_p6, 		int size_idx_p5, 		int size_idx_p4, 		int size_idx_p7)
{
	//
	int	str_sd2_t3_h3 = 1;
	int str_sd2_t3_h2 = str_sd2_t3_h3 * size_idx_h3;
	int str_sd2_t3_h1 = str_sd2_t3_h2 * size_idx_h2;
	int str_sd2_t3_p6 = str_sd2_t3_h1 * size_idx_h1;
	int str_sd2_t3_p5 = str_sd2_t3_p6 * size_idx_p6;
	int str_sd2_t3_p4 = str_sd2_t3_p5 * size_idx_p5;

	// tc_gen_code_pre_IndirectArray_Init()
	h_t3_output_base_addr 	= (int*)malloc(sizeof(int) * (n_blocks_2));
	h_t2_7_ext_addr 		= (int*)malloc(sizeof(int) * SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H2 * n_blocks_2);
	h_v2_7_ext_addr 		= (int*)malloc(sizeof(int) * SIZE_SLICE_2_H3 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4 * n_blocks_2);
	h_t2_8_ext_addr 		= (int*)malloc(sizeof(int) * SIZE_SLICE_2_P6 * SIZE_SLICE_2_H2 * SIZE_SLICE_2_H3 * n_blocks_2);
	h_v2_8_ext_addr 		= (int*)malloc(sizeof(int) * SIZE_SLICE_2_H1 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4 * n_blocks_2);
	h_t2_9_ext_addr 		= (int*)malloc(sizeof(int) * SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H3 * n_blocks_2);
	h_v2_9_ext_addr 		= (int*)malloc(sizeof(int) * SIZE_SLICE_2_H2 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4 * n_blocks_2);

	// tc_gen_code_pre_IndirectArray_wo_TB()
	// For Each Thread Block
	for (int i = 0; i < n_blocks_2; i++)
	{
		int blk_idx_h3 = h_t3_block_index[i * NUM_INDEX + 0];
		int blk_idx_h2 = h_t3_block_index[i * NUM_INDEX + 1];
		int blk_idx_h1 = h_t3_block_index[i * NUM_INDEX + 2];
		int blk_idx_p6 = h_t3_block_index[i * NUM_INDEX + 3];
		int blk_idx_p5 = h_t3_block_index[i * NUM_INDEX + 4];
		int blk_idx_p4 = h_t3_block_index[i * NUM_INDEX + 5];

		// Calculating t3
		h_t3_output_base_addr[i] = 	blk_idx_h3 * SIZE_SLICE_2_H3 * str_sd2_t3_h3 + blk_idx_h2 * SIZE_SLICE_2_H2 * str_sd2_t3_h2 + blk_idx_h1 * SIZE_SLICE_2_H1 * str_sd2_t3_h1 + 
									blk_idx_p6 * SIZE_SLICE_2_P6 * str_sd2_t3_p6 + blk_idx_p5 * SIZE_SLICE_2_P5 * str_sd2_t3_p5 + blk_idx_p4 * SIZE_SLICE_2_P4 * str_sd2_t3_p4;

		// calculating actual address for t2_7 and v2_7
		for (int idx_p6 = 0; idx_p6 < SIZE_SLICE_2_P6; idx_p6++)
		for (int idx_h2 = 0; idx_h2 < SIZE_SLICE_2_H2; idx_h2++)
		for (int idx_h1 = 0; idx_h1 < SIZE_SLICE_2_H1; idx_h1++)
		{
			int l_offset = idx_h1 + (idx_h2 + (idx_p6) * SIZE_SLICE_2_H2) * SIZE_SLICE_2_H1;
			h_t2_7_ext_addr[l_offset + i * SIZE_SLICE_2_P6 * SIZE_SLICE_2_H2 * SIZE_SLICE_2_H1] = (blk_idx_p6 * SIZE_SLICE_2_P6 + idx_p6 + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_h1 + 
																								(blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h2) * size_idx_h1) * size_idx_p6) * size_idx_p7;
		}

		for (int idx_p5 = 0; idx_p5 < SIZE_SLICE_2_P5; idx_p5++)
		for (int idx_p4 = 0; idx_p4 < SIZE_SLICE_2_P4; idx_p4++)
		for (int idx_h3 = 0; idx_h3 < SIZE_SLICE_2_H3; idx_h3++)
		{
			int l_offset = idx_h3 + (idx_p4 + (idx_p5) * SIZE_SLICE_2_P4) * SIZE_SLICE_2_H3;
			h_v2_7_ext_addr[l_offset + i * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4 * SIZE_SLICE_2_H3] = (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h3 + (blk_idx_p5 * SIZE_SLICE_2_P5 + idx_p5 + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_p4) * size_idx_p5) * size_idx_h3) * size_idx_p7;
		}
		// calculating actual address for t2_8 and v2_8
		for (int idx_p6 = 0; idx_p6 < SIZE_SLICE_2_P6; idx_p6++)
		for (int idx_h3 = 0; idx_h3 < SIZE_SLICE_2_H3; idx_h3++)
		for (int idx_h2 = 0; idx_h2 < SIZE_SLICE_2_H2; idx_h2++)
		{
			int l_offset = idx_h2 + (idx_h3 + (idx_p6) * SIZE_SLICE_2_H3) * SIZE_SLICE_2_H2;
			h_t2_8_ext_addr[l_offset + i * SIZE_SLICE_2_P6 * SIZE_SLICE_2_H3 * SIZE_SLICE_2_H2] = (blk_idx_p6 * SIZE_SLICE_2_P6 + idx_p6 + (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h2 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h3) * size_idx_h2) * size_idx_p6) * size_idx_p7;
		}

		for (int idx_p5 = 0; idx_p5 < SIZE_SLICE_2_P5; idx_p5++)
		for (int idx_p4 = 0; idx_p4 < SIZE_SLICE_2_P4; idx_p4++)
		for (int idx_h1 = 0; idx_h1 < SIZE_SLICE_2_H1; idx_h1++)
		{
			int l_offset = idx_h1 + (idx_p4 + (idx_p5) * SIZE_SLICE_2_P4) * SIZE_SLICE_2_H1;
			h_v2_8_ext_addr[l_offset + i * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4 * SIZE_SLICE_2_H1] = (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_h1 + (blk_idx_p5 * SIZE_SLICE_2_P5 + idx_p5 + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_p4) * size_idx_p5) * size_idx_h1) * size_idx_p7;
		}
		// calculating actual address for t2_9 and v2_9
		for (int idx_p6 = 0; idx_p6 < SIZE_SLICE_2_P6; idx_p6++)
		for (int idx_h3 = 0; idx_h3 < SIZE_SLICE_2_H3; idx_h3++)
		for (int idx_h1 = 0; idx_h1 < SIZE_SLICE_2_H1; idx_h1++)
		{
			int l_offset = idx_h1 + (idx_h3 + (idx_p6) * SIZE_SLICE_2_H3) * SIZE_SLICE_2_H1;
			h_t2_9_ext_addr[l_offset + i * SIZE_SLICE_2_P6 * SIZE_SLICE_2_H3 * SIZE_SLICE_2_H1] = (blk_idx_p6 * SIZE_SLICE_2_P6 + idx_p6 + (blk_idx_h1 * SIZE_SLICE_2_H1 + idx_h1 + (blk_idx_h3 * SIZE_SLICE_2_H3 + idx_h3) * size_idx_h1) * size_idx_p6) * size_idx_p7;
		}

		for (int idx_p5 = 0; idx_p5 < SIZE_SLICE_2_P5; idx_p5++)
		for (int idx_p4 = 0; idx_p4 < SIZE_SLICE_2_P4; idx_p4++)
		for (int idx_h2 = 0; idx_h2 < SIZE_SLICE_2_H2; idx_h2++)
		{
			int l_offset = idx_h2 + (idx_p4 + (idx_p5) * SIZE_SLICE_2_P4) * SIZE_SLICE_2_H2;
			h_v2_9_ext_addr[l_offset + i * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4 * SIZE_SLICE_2_H2] = (blk_idx_h2 * SIZE_SLICE_2_H2 + idx_h2 + (blk_idx_p5 * SIZE_SLICE_2_P5 + idx_p5 + (blk_idx_p4 * SIZE_SLICE_2_P4 + idx_p4) * size_idx_p5) * size_idx_h2) * size_idx_p7;
		}
	}

	// tc_gen_code_pre_IndirectArray_w_TB()
	// Within a Thread Block
	h_t3_output_ext_offset 	= (int*)malloc(sizeof(int) * (SIZE_TB_2_X * SIZE_TB_2_Y));
	h_t2_7_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_2_X * SIZE_TB_2_Y));
	h_v2_7_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_2_X * SIZE_TB_2_Y));
	h_t2_8_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_2_X * SIZE_TB_2_Y));
	h_v2_8_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_2_X * SIZE_TB_2_Y));
	h_t2_9_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_2_X * SIZE_TB_2_Y));
	h_v2_9_ext_offset 		= (int*)malloc(sizeof(int) * (SIZE_TB_2_X * SIZE_TB_2_Y));

	for (int idx_h1 = 0; idx_h1 < SIZE_SLICE_2_H1; idx_h1++)
	for (int idx_p4 = 0; idx_p4 < SIZE_SLICE_2_P4; idx_p4++)
	for (int idx_h2 = 0; idx_h2 < SIZE_SLICE_2_H2; idx_h2++)
	for (int idx_h3 = 0; idx_h3 < SIZE_SLICE_2_H3; idx_h3++)
	{
		int t3_index_1 = idx_h3 + (idx_h2 + (idx_p4 + (idx_h1) * SIZE_SLICE_2_P4) * SIZE_SLICE_2_H2) * SIZE_SLICE_2_H3;

		h_t3_output_ext_offset[t3_index_1] = idx_h3 + (idx_h2 + (idx_h1 + (((idx_p4) * size_idx_p5) * size_idx_p6) * size_idx_h1) * size_idx_h2) * size_idx_h3;
		h_t2_7_ext_offset[t3_index_1] = idx_h1 + (idx_h2) * SIZE_SLICE_2_H1;
		h_v2_7_ext_offset[t3_index_1] = idx_h3 + (idx_p4) * SIZE_SLICE_2_H3;
		h_t2_8_ext_offset[t3_index_1] = idx_h2 + (idx_h3) * SIZE_SLICE_2_H2;
		h_v2_8_ext_offset[t3_index_1] = idx_h1 + (idx_p4) * SIZE_SLICE_2_H1;
		h_t2_9_ext_offset[t3_index_1] = idx_h1 + (idx_h3) * SIZE_SLICE_2_H1;
		h_v2_9_ext_offset[t3_index_1] = idx_h2 + (idx_p4) * SIZE_SLICE_2_H2;
	}
	// Do not need Indirect Arrays for Internal Indices
}

void Pre_TileApproach_1(int*& h_t3_block_index, int*& h_t3_block_range, int* n_blocks_1, int size_h3, int size_h2, int size_h1, int size_p6, int size_p5, int size_p4, int size_p7)
{
	int n_blk_h3 = CEIL(size_h3, 4);
	int n_blk_h2 = CEIL(size_h2, 4);
	int n_blk_h1 = CEIL(size_h1, 4);
	int n_blk_p6 = CEIL(size_p6, 4);
	int n_blk_p5 = CEIL(size_p5, 4);
	int n_blk_p4 = CEIL(size_p4, 4);

	int blk_h3_range[n_blk_h3];
	int blk_h2_range[n_blk_h2];
	int blk_h1_range[n_blk_h1];
	int blk_p6_range[n_blk_p6];
	int blk_p5_range[n_blk_p5];
	int blk_p4_range[n_blk_p4];

	int rng_boundary_h3 = size_h3 % 4;
	int rng_boundary_h2 = size_h2 % 4;
	int rng_boundary_h1 = size_h1 % 4;
	int rng_boundary_p6 = size_p6 % 4;
	int rng_boundary_p5 = size_p5 % 4;
	int rng_boundary_p4 = size_p4 % 4;

	for (int i = 0; i < n_blk_h3; i++)
	{
		blk_h3_range[i] = 4;
		if (rng_boundary_h3 != 0 && i == n_blk_h3 - 1)
		{
			blk_h3_range[i] = rng_boundary_h3;
		}
	}
	for (int i = 0; i < n_blk_h2; i++)
	{
		blk_h2_range[i] = 4;
		if (rng_boundary_h2 != 0 && i == n_blk_h2 - 1)
		{
			blk_h2_range[i] = rng_boundary_h2;
		}
	}
	for (int i = 0; i < n_blk_h1; i++)
	{
		blk_h1_range[i] = 4;
		if (rng_boundary_h1 != 0 && i == n_blk_h1 - 1)
		{
			blk_h1_range[i] = rng_boundary_h1;
		}
	}
	for (int i = 0; i < n_blk_p6; i++)
	{
		blk_p6_range[i] = 4;
		if (rng_boundary_p6 != 0 && i == n_blk_p6 - 1)
		{
			blk_p6_range[i] = rng_boundary_p6;
		}
	}
	for (int i = 0; i < n_blk_p5; i++)
	{
		blk_p5_range[i] = 4;
		if (rng_boundary_p5 != 0 && i == n_blk_p5 - 1)
		{
			blk_p5_range[i] = rng_boundary_p5;
		}
	}
	for (int i = 0; i < n_blk_p4; i++)
	{
		blk_p4_range[i] = 4;
		if (rng_boundary_p4 != 0 && i == n_blk_p4 - 1)
		{
			blk_p4_range[i] = rng_boundary_p4;
		}
	}

	// # of blocks
	*n_blocks_1 = n_blk_h3 * n_blk_h2 * n_blk_h1 * n_blk_p6 * n_blk_p5 * n_blk_p4;
	h_t3_block_index = (int*)malloc(sizeof(int) * (*n_blocks_1) * NUM_INDEX);
	h_t3_block_range = (int*)malloc(sizeof(int) * (*n_blocks_1) * NUM_INDEX);

	for (int blk_h3 = 0; blk_h3 < n_blk_h3; blk_h3++)
	for (int blk_h2 = 0; blk_h2 < n_blk_h2; blk_h2++)
	for (int blk_h1 = 0; blk_h1 < n_blk_h1; blk_h1++)
	for (int blk_p6 = 0; blk_p6 < n_blk_p6; blk_p6++)
	for (int blk_p5 = 0; blk_p5 < n_blk_p5; blk_p5++)
	for (int blk_p4 = 0; blk_p4 < n_blk_p4; blk_p4++)
	{
		int blk_idx_base = (blk_h3 + (blk_h2 + (blk_h1 + (blk_p6 + (blk_p5 + (blk_p4) * n_blk_p5) * n_blk_p6) * n_blk_h1) * n_blk_h2) * n_blk_h3) * NUM_INDEX;
		h_t3_block_index[blk_idx_base + 0] = blk_h3;
		h_t3_block_index[blk_idx_base + 1] = blk_h2;
		h_t3_block_index[blk_idx_base + 2] = blk_h1;
		h_t3_block_index[blk_idx_base + 3] = blk_p6;
		h_t3_block_index[blk_idx_base + 4] = blk_p5;
		h_t3_block_index[blk_idx_base + 5] = blk_p4;

		h_t3_block_range[blk_idx_base + 0] = blk_h3_range[blk_h3];
		h_t3_block_range[blk_idx_base + 1] = blk_h2_range[blk_h2];
		h_t3_block_range[blk_idx_base + 2] = blk_h1_range[blk_h1];
		h_t3_block_range[blk_idx_base + 3] = blk_p6_range[blk_p6];
		h_t3_block_range[blk_idx_base + 4] = blk_p5_range[blk_p5];
		h_t3_block_range[blk_idx_base + 5] = blk_p4_range[blk_p4];
	}
}

//
void Pre_TileApproach_2(int*& h_t3_block_index, int*& h_t3_block_range, int* n_blocks_2, int size_h3, int size_h2, int size_h1, int size_p6, int size_p5, int size_p4, int size_p7)
{
	//
	int n_blk_h3 = CEIL(size_h3, 4);
	int n_blk_h2 = CEIL(size_h2, 4);
	int n_blk_h1 = CEIL(size_h1, 4);
	int n_blk_p6 = CEIL(size_p6, 4);
	int n_blk_p5 = CEIL(size_p5, 4);
	int n_blk_p4 = CEIL(size_p4, 4);

	int blk_h3_range[n_blk_h3];
	int blk_h2_range[n_blk_h2];
	int blk_h1_range[n_blk_h1];
	int blk_p6_range[n_blk_p6];
	int blk_p5_range[n_blk_p5];
	int blk_p4_range[n_blk_p4];

	int rng_boundary_h3 = size_h3 % 4;
	int rng_boundary_h2 = size_h2 % 4;
	int rng_boundary_h1 = size_h1 % 4;
	int rng_boundary_p6 = size_p6 % 4;
	int rng_boundary_p5 = size_p5 % 4;
	int rng_boundary_p4 = size_p4 % 4;

	for (int i = 0; i < n_blk_h3; i++)
	{
		blk_h3_range[i] = 4;
		if (rng_boundary_h3 != 0 && i == n_blk_h3 - 1)
		{
			blk_h3_range[i] = rng_boundary_h3;
		}
	}
	for (int i = 0; i < n_blk_h2; i++)
	{
		blk_h2_range[i] = 4;
		if (rng_boundary_h2 != 0 && i == n_blk_h2 - 1)
		{
			blk_h2_range[i] = rng_boundary_h2;
		}
	}
	for (int i = 0; i < n_blk_h1; i++)
	{
		blk_h1_range[i] = 4;
		if (rng_boundary_h1 != 0 && i == n_blk_h1 - 1)
		{
			blk_h1_range[i] = rng_boundary_h1;
		}
	}
	for (int i = 0; i < n_blk_p6; i++)
	{
		blk_p6_range[i] = 4;
		if (rng_boundary_p6 != 0 && i == n_blk_p6 - 1)
		{
			blk_p6_range[i] = rng_boundary_p6;
		}
	}
	for (int i = 0; i < n_blk_p5; i++)
	{
		blk_p5_range[i] = 4;
		if (rng_boundary_p5 != 0 && i == n_blk_p5 - 1)
		{
			blk_p5_range[i] = rng_boundary_p5;
		}
	}
	for (int i = 0; i < n_blk_p4; i++)
	{
		blk_p4_range[i] = 4;
		if (rng_boundary_p4 != 0 && i == n_blk_p4 - 1)
		{
			blk_p4_range[i] = rng_boundary_p4;
		}
	}

	// # of blocks
	*n_blocks_2 = n_blk_h3 * n_blk_h2 * n_blk_h1 * n_blk_p6 * n_blk_p5 * n_blk_p4;
	h_t3_block_index = (int*)malloc(sizeof(int) * (*n_blocks_2) * NUM_INDEX);
	h_t3_block_range = (int*)malloc(sizeof(int) * (*n_blocks_2) * NUM_INDEX);

	for (int blk_h3 = 0; blk_h3 < n_blk_h3; blk_h3++)
	for (int blk_h2 = 0; blk_h2 < n_blk_h2; blk_h2++)
	for (int blk_h1 = 0; blk_h1 < n_blk_h1; blk_h1++)
	for (int blk_p6 = 0; blk_p6 < n_blk_p6; blk_p6++)
	for (int blk_p5 = 0; blk_p5 < n_blk_p5; blk_p5++)
	for (int blk_p4 = 0; blk_p4 < n_blk_p4; blk_p4++)
	{
		int blk_idx_base = (blk_h3 + (blk_h2 + (blk_h1 + (blk_p6 + (blk_p5 + (blk_p4) * n_blk_p5) * n_blk_p6) * n_blk_h1) * n_blk_h2) * n_blk_h3) * NUM_INDEX;
		h_t3_block_index[blk_idx_base + 0] = blk_h3;
		h_t3_block_index[blk_idx_base + 1] = blk_h2;
		h_t3_block_index[blk_idx_base + 2] = blk_h1;
		h_t3_block_index[blk_idx_base + 3] = blk_p6;
		h_t3_block_index[blk_idx_base + 4] = blk_p5;
		h_t3_block_index[blk_idx_base + 5] = blk_p4;

		h_t3_block_range[blk_idx_base + 0] = blk_h3_range[blk_h3];
		h_t3_block_range[blk_idx_base + 1] = blk_h2_range[blk_h2];
		h_t3_block_range[blk_idx_base + 2] = blk_h1_range[blk_h1];
		h_t3_block_range[blk_idx_base + 3] = blk_p6_range[blk_p6];
		h_t3_block_range[blk_idx_base + 4] = blk_p5_range[blk_p5];
		h_t3_block_range[blk_idx_base + 5] = blk_p4_range[blk_p4];
	}
}

// created by tc_gen_code_Kernel()
__global__ void kernel_ccsdT_0(double* t3, 
									const int* __restrict__ t3_blk_rng_1, const int* __restrict__ t3_output_base_1, const int* __restrict__ t3_output_offset_1, 
									double* d_t2_1, const int* __restrict__ d_t2_1_addr, const int* __restrict__ d_t2_1_offset, 
									double* d_t2_2, const int* __restrict__ d_t2_2_addr, const int* __restrict__ d_t2_2_offset, 
									double* d_t2_3, const int* __restrict__ d_t2_3_addr, const int* __restrict__ d_t2_3_offset, 
									double* d_t2_4, const int* __restrict__ d_t2_4_addr, const int* __restrict__ d_t2_4_offset, 
									double* d_t2_5, const int* __restrict__ d_t2_5_addr, const int* __restrict__ d_t2_5_offset, 
									double* d_t2_6, const int* __restrict__ d_t2_6_addr, const int* __restrict__ d_t2_6_offset, 
									double* d_v2_1, const int* __restrict__ d_v2_1_addr, const int* __restrict__ d_v2_1_offset, 
									double* d_v2_2, const int* __restrict__ d_v2_2_addr, const int* __restrict__ d_v2_2_offset, 
									double* d_v2_3, const int* __restrict__ d_v2_3_addr, const int* __restrict__ d_v2_3_offset, 
									double* d_v2_4, const int* __restrict__ d_v2_4_addr, const int* __restrict__ d_v2_4_offset, 
									double* d_v2_5, const int* __restrict__ d_v2_5_addr, const int* __restrict__ d_v2_5_offset, 
									double* d_v2_6, const int* __restrict__ d_v2_6_addr, const int* __restrict__ d_v2_6_offset, 
									const int* __restrict__ t3_blk_rng_2, const int* __restrict__ t3_output_base_2, const int* __restrict__ t3_output_offset_2, 
									double* d_t2_7, const int* __restrict__ d_t2_7_addr, const int* __restrict__ d_t2_7_offset, 
									double* d_t2_8, const int* __restrict__ d_t2_8_addr, const int* __restrict__ d_t2_8_offset, 
									double* d_t2_9, const int* __restrict__ d_t2_9_addr, const int* __restrict__ d_t2_9_offset, 
									double* d_v2_7, const int* __restrict__ d_v2_7_addr, const int* __restrict__ d_v2_7_offset, 
									double* d_v2_8, const int* __restrict__ d_v2_8_addr, const int* __restrict__ d_v2_8_offset, 
									double* d_v2_9, const int* __restrict__ d_v2_9_addr, const int* __restrict__ d_v2_9_offset, 
									int kernel_1, int kernel_2, int kernel_3, 
									int kernel_4, int kernel_5, int kernel_6, 
									int kernel_7, int kernel_8, int kernel_9, 
									int stride_reg_x, int stride_reg_y,
									int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int l_idx_t3         		= threadIdx.x + threadIdx.y * SIZE_TB_1_X;
	int t3_base_thread   		= t3_output_base_1[blockIdx.x] + t3_output_offset_1[l_idx_t3];
	int internal_upperbound   	= 0;
	int internal_offset;

	// should support for non-full tiles
	int idx_h3 = threadIdx.x % SIZE_SLICE_1_H3;
	int idx_h2 = threadIdx.x / SIZE_SLICE_1_H3;
	int idx_p6 = threadIdx.y % SIZE_SLICE_1_P6;
	int idx_h1 = threadIdx.y / SIZE_SLICE_1_P6;

	int rng_h3 = t3_blk_rng_1[blockIdx.x * NUM_INDEX + 0];
	int rng_h2 = t3_blk_rng_1[blockIdx.x * NUM_INDEX + 1];
	int rng_h1 = t3_blk_rng_1[blockIdx.x * NUM_INDEX + 2];
	int rng_p6 = t3_blk_rng_1[blockIdx.x * NUM_INDEX + 3];
	int rng_p5 = t3_blk_rng_1[blockIdx.x * NUM_INDEX + 4];
	int rng_p4 = t3_blk_rng_1[blockIdx.x * NUM_INDEX + 5];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_7[d_t2_7_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h3 && idx_h1 < rng_p4 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_7[d_v2_7_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H3 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_7_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_8[d_t2_8_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H2 * SIZE_SLICE_2_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h1 && idx_h1 < rng_p4 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_8[d_v2_8_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H1 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_8_offset[l_idx_t3] + (xx * 16)];

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
		if (idx_p6 < rng_h1 && idx_h1 < rng_h3 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p6; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_9[d_t2_9_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h2 && idx_h1 < rng_p4 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_9[d_v2_9_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H2 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_9_offset[l_idx_t3] + (xx * 16)];

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
		if (idx_p6 < rng_h1 && idx_h1 < rng_h2 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_1[d_t2_1_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h3 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_1[d_v2_1_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_1_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_2[d_t2_2_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h1 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_2[d_v2_2_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_2_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_3[d_t2_3_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h2 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_3[d_v2_3_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_3_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_4[d_t2_4_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h3 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_4[d_v2_4_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_4_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_5[d_t2_5_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h1 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_5[d_v2_5_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_5_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_6[d_t2_6_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h2 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_6[d_v2_6_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_6_offset[l_idx_t3] + (xx * 16)];

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
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] += reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void kernel_ccsdT_0_full(double* t3, 
									const int* __restrict__ t3_blk_rng_1, const int* __restrict__ t3_output_base_1, const int* __restrict__ t3_output_offset_1, 
									double* d_t2_1, const int* __restrict__ d_t2_1_addr, const int* __restrict__ d_t2_1_offset, 
									double* d_t2_2, const int* __restrict__ d_t2_2_addr, const int* __restrict__ d_t2_2_offset, 
									double* d_t2_3, const int* __restrict__ d_t2_3_addr, const int* __restrict__ d_t2_3_offset, 
									double* d_t2_4, const int* __restrict__ d_t2_4_addr, const int* __restrict__ d_t2_4_offset, 
									double* d_t2_5, const int* __restrict__ d_t2_5_addr, const int* __restrict__ d_t2_5_offset, 
									double* d_t2_6, const int* __restrict__ d_t2_6_addr, const int* __restrict__ d_t2_6_offset, 
									double* d_v2_1, const int* __restrict__ d_v2_1_addr, const int* __restrict__ d_v2_1_offset, 
									double* d_v2_2, const int* __restrict__ d_v2_2_addr, const int* __restrict__ d_v2_2_offset, 
									double* d_v2_3, const int* __restrict__ d_v2_3_addr, const int* __restrict__ d_v2_3_offset, 
									double* d_v2_4, const int* __restrict__ d_v2_4_addr, const int* __restrict__ d_v2_4_offset, 
									double* d_v2_5, const int* __restrict__ d_v2_5_addr, const int* __restrict__ d_v2_5_offset, 
									double* d_v2_6, const int* __restrict__ d_v2_6_addr, const int* __restrict__ d_v2_6_offset, 
									const int* __restrict__ t3_blk_rng_2, const int* __restrict__ t3_output_base_2, const int* __restrict__ t3_output_offset_2, 
									double* d_t2_7, const int* __restrict__ d_t2_7_addr, const int* __restrict__ d_t2_7_offset, 
									double* d_t2_8, const int* __restrict__ d_t2_8_addr, const int* __restrict__ d_t2_8_offset, 
									double* d_t2_9, const int* __restrict__ d_t2_9_addr, const int* __restrict__ d_t2_9_offset, 
									double* d_v2_7, const int* __restrict__ d_v2_7_addr, const int* __restrict__ d_v2_7_offset, 
									double* d_v2_8, const int* __restrict__ d_v2_8_addr, const int* __restrict__ d_v2_8_offset, 
									double* d_v2_9, const int* __restrict__ d_v2_9_addr, const int* __restrict__ d_v2_9_offset, 
									int kernel_1, int kernel_2, int kernel_3, 
									int kernel_4, int kernel_5, int kernel_6, 
									int kernel_7, int kernel_8, int kernel_9, 
									int stride_reg_x, int stride_reg_y,
									int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int l_idx_t3         		= threadIdx.x + threadIdx.y * SIZE_TB_1_X;
	int t3_base_thread   		= t3_output_base_1[blockIdx.x] + t3_output_offset_1[l_idx_t3];
	int internal_upperbound   	= 0;
	int internal_offset;

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_7[d_t2_7_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_7[d_v2_7_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H3 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_7_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_8[d_t2_8_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H2 * SIZE_SLICE_2_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_8[d_v2_8_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H1 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_8_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_9[d_t2_9_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_9[d_v2_9_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H2 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_9_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_1[d_t2_1_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_1[d_v2_1_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_1_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_2[d_t2_2_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_2[d_v2_2_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_2_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_3[d_t2_3_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_3[d_v2_3_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_3_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_4[d_t2_4_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_4[d_v2_4_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_4_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_5[d_t2_5_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_5[d_v2_5_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_5_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_6[d_t2_6_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_6[d_v2_6_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_6_offset[l_idx_t3] + (xx * 16)];

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
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] += reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void kernel_ccsdT_0_full_full(double* t3, 
									const int* __restrict__ t3_blk_rng_1, const int* __restrict__ t3_output_base_1, const int* __restrict__ t3_output_offset_1, 
									double* d_t2_1, const int* __restrict__ d_t2_1_addr, const int* __restrict__ d_t2_1_offset, 
									double* d_t2_2, const int* __restrict__ d_t2_2_addr, const int* __restrict__ d_t2_2_offset, 
									double* d_t2_3, const int* __restrict__ d_t2_3_addr, const int* __restrict__ d_t2_3_offset, 
									double* d_t2_4, const int* __restrict__ d_t2_4_addr, const int* __restrict__ d_t2_4_offset, 
									double* d_t2_5, const int* __restrict__ d_t2_5_addr, const int* __restrict__ d_t2_5_offset, 
									double* d_t2_6, const int* __restrict__ d_t2_6_addr, const int* __restrict__ d_t2_6_offset, 
									double* d_v2_1, const int* __restrict__ d_v2_1_addr, const int* __restrict__ d_v2_1_offset, 
									double* d_v2_2, const int* __restrict__ d_v2_2_addr, const int* __restrict__ d_v2_2_offset, 
									double* d_v2_3, const int* __restrict__ d_v2_3_addr, const int* __restrict__ d_v2_3_offset, 
									double* d_v2_4, const int* __restrict__ d_v2_4_addr, const int* __restrict__ d_v2_4_offset, 
									double* d_v2_5, const int* __restrict__ d_v2_5_addr, const int* __restrict__ d_v2_5_offset, 
									double* d_v2_6, const int* __restrict__ d_v2_6_addr, const int* __restrict__ d_v2_6_offset, 
									const int* __restrict__ t3_blk_rng_2, const int* __restrict__ t3_output_base_2, const int* __restrict__ t3_output_offset_2, 
									double* d_t2_7, const int* __restrict__ d_t2_7_addr, const int* __restrict__ d_t2_7_offset, 
									double* d_t2_8, const int* __restrict__ d_t2_8_addr, const int* __restrict__ d_t2_8_offset, 
									double* d_t2_9, const int* __restrict__ d_t2_9_addr, const int* __restrict__ d_t2_9_offset, 
									double* d_v2_7, const int* __restrict__ d_v2_7_addr, const int* __restrict__ d_v2_7_offset, 
									double* d_v2_8, const int* __restrict__ d_v2_8_addr, const int* __restrict__ d_v2_8_offset, 
									double* d_v2_9, const int* __restrict__ d_v2_9_addr, const int* __restrict__ d_v2_9_offset, 
									int kernel_1, int kernel_2, int kernel_3, 
									int kernel_4, int kernel_5, int kernel_6, 
									int kernel_7, int kernel_8, int kernel_9, 
									int stride_reg_x, int stride_reg_y,
									int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int l_idx_t3         		= threadIdx.x + threadIdx.y * SIZE_TB_1_X;
	int t3_base_thread   		= t3_output_base_1[blockIdx.x] + t3_output_offset_1[l_idx_t3];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_7[d_t2_7_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_7[d_v2_7_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H3 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_7_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_8[d_t2_8_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H2 * SIZE_SLICE_2_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_8[d_v2_8_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H1 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_8_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_9[d_t2_9_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_9[d_v2_9_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H2 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_9_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_1[d_t2_1_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_1[d_v2_1_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_1_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_2[d_t2_2_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_2[d_v2_2_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_2_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_3[d_t2_3_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_3[d_v2_3_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_3_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_4[d_t2_4_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_4[d_v2_4_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_4_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_5[d_t2_5_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_5[d_v2_5_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_5_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_6[d_t2_6_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_6[d_v2_6_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_6_offset[l_idx_t3] + (xx * 16)];

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
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] += reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void kernel_ccsdT_1(double* t3, const int* __restrict__ t3_blk_rng_1, const int* __restrict__ t3_output_base_1, const int* __restrict__ t3_output_offset_1, 
double* d_t2_1, const int* __restrict__ d_t2_1_addr, const int* __restrict__ d_t2_1_offset, double* d_t2_2, const int* __restrict__ d_t2_2_addr, const int* __restrict__ d_t2_2_offset, double* d_t2_3, const int* __restrict__ d_t2_3_addr, const int* __restrict__ d_t2_3_offset, double* d_t2_4, const int* __restrict__ d_t2_4_addr, const int* __restrict__ d_t2_4_offset, double* d_t2_5, const int* __restrict__ d_t2_5_addr, const int* __restrict__ d_t2_5_offset, double* d_t2_6, const int* __restrict__ d_t2_6_addr, const int* __restrict__ d_t2_6_offset, 
double* d_v2_1, const int* __restrict__ d_v2_1_addr, const int* __restrict__ d_v2_1_offset, double* d_v2_2, const int* __restrict__ d_v2_2_addr, const int* __restrict__ d_v2_2_offset, double* d_v2_3, const int* __restrict__ d_v2_3_addr, const int* __restrict__ d_v2_3_offset, double* d_v2_4, const int* __restrict__ d_v2_4_addr, const int* __restrict__ d_v2_4_offset, double* d_v2_5, const int* __restrict__ d_v2_5_addr, const int* __restrict__ d_v2_5_offset, double* d_v2_6, const int* __restrict__ d_v2_6_addr, const int* __restrict__ d_v2_6_offset, 
int kernel_1, int kernel_2, int kernel_3, 
int kernel_4, int kernel_5, int kernel_6, 
int kernel_7, int kernel_8, int kernel_9, 
int stride_reg_x, int stride_reg_y,
int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int l_idx_t3         		= threadIdx.x + threadIdx.y * SIZE_TB_1_X;
	int t3_base_thread   		= t3_output_base_1[blockIdx.x] + t3_output_offset_1[l_idx_t3];
	int internal_upperbound   	= 0;
	int internal_offset;


	// should support for non-full tiles
	int idx_h3 = threadIdx.x % SIZE_SLICE_1_H3;
	int idx_h2 = threadIdx.x / SIZE_SLICE_1_H3;
	int idx_p6 = threadIdx.y % SIZE_SLICE_1_P6;
	int idx_h1 = threadIdx.y / SIZE_SLICE_1_P6;	// SIZE_SLICE_1_P6 == SIZE_SLICE_1_H3.

	// Common for Threads within a Thread Block
	int rng_h3 = t3_blk_rng_1[blockIdx.x * NUM_INDEX + 0];
	int rng_h2 = t3_blk_rng_1[blockIdx.x * NUM_INDEX + 1];
	int rng_h1 = t3_blk_rng_1[blockIdx.x * NUM_INDEX + 2];
	int rng_p6 = t3_blk_rng_1[blockIdx.x * NUM_INDEX + 3];
	int rng_p5 = t3_blk_rng_1[blockIdx.x * NUM_INDEX + 4];
	int rng_p4 = t3_blk_rng_1[blockIdx.x * NUM_INDEX + 5];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_1[d_t2_1_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h3 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_1[d_v2_1_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_1_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_2[d_t2_2_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h1 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_2[d_v2_2_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_2_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_3[d_t2_3_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h2 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_3[d_v2_3_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_3_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_4[d_t2_4_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h3 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_4[d_v2_4_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_4_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_5[d_t2_5_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h1 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_5[d_v2_5_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_5_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_6[d_t2_6_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p6 < rng_h2 && idx_h1 < rng_p6 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_6[d_v2_6_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 48];

			for (int xx = 0; xx < 4; xx++)	// 4 -> rng_p4: Local Transactions...
			{
				temp_av = sm_a[ll][d_t2_6_offset[l_idx_t3] + (xx * 16)];

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
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] += reg_tile[i][j];
		}
	}
}

__global__ void kernel_ccsdT_1_full(double* t3, const int* __restrict__ t3_blk_rng_1, const int* __restrict__ t3_output_base_1, const int* __restrict__ t3_output_offset_1, 
double* d_t2_1, const int* __restrict__ d_t2_1_addr, const int* __restrict__ d_t2_1_offset, double* d_t2_2, const int* __restrict__ d_t2_2_addr, const int* __restrict__ d_t2_2_offset, double* d_t2_3, const int* __restrict__ d_t2_3_addr, const int* __restrict__ d_t2_3_offset, double* d_t2_4, const int* __restrict__ d_t2_4_addr, const int* __restrict__ d_t2_4_offset, double* d_t2_5, const int* __restrict__ d_t2_5_addr, const int* __restrict__ d_t2_5_offset, double* d_t2_6, const int* __restrict__ d_t2_6_addr, const int* __restrict__ d_t2_6_offset, 
double* d_v2_1, const int* __restrict__ d_v2_1_addr, const int* __restrict__ d_v2_1_offset, double* d_v2_2, const int* __restrict__ d_v2_2_addr, const int* __restrict__ d_v2_2_offset, double* d_v2_3, const int* __restrict__ d_v2_3_addr, const int* __restrict__ d_v2_3_offset, double* d_v2_4, const int* __restrict__ d_v2_4_addr, const int* __restrict__ d_v2_4_offset, double* d_v2_5, const int* __restrict__ d_v2_5_addr, const int* __restrict__ d_v2_5_offset, double* d_v2_6, const int* __restrict__ d_v2_6_addr, const int* __restrict__ d_v2_6_offset, 
int kernel_1, int kernel_2, int kernel_3, 
int kernel_4, int kernel_5, int kernel_6, 
int kernel_7, int kernel_8, int kernel_9, 	// 7, 8, 9 are useless
int stride_reg_x, int stride_reg_y,
int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int l_idx_t3         = threadIdx.x + threadIdx.y * SIZE_TB_1_X;
	int t3_base_thread   = t3_output_base_1[blockIdx.x] + t3_output_offset_1[l_idx_t3];
	int internal_upperbound   = 0;
	int internal_offset;

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_1[d_t2_1_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_1[d_v2_1_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_1_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_2[d_t2_2_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_2[d_v2_2_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_2_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_3[d_t2_3_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_3[d_v2_3_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_3_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_4[d_t2_4_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_4[d_v2_4_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_4_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_5[d_t2_5_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_5[d_v2_5_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_5_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_6[d_t2_6_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_6[d_v2_6_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_6_offset[l_idx_t3] + (xx * 16)];

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
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] += reg_tile[i][j];
		}
	}
}

__global__ void kernel_ccsdT_1_full_full(double* t3, const int* __restrict__ t3_blk_rng_1, const int* __restrict__ t3_output_base_1, const int* __restrict__ t3_output_offset_1, 
double* d_t2_1, const int* __restrict__ d_t2_1_addr, const int* __restrict__ d_t2_1_offset, double* d_t2_2, const int* __restrict__ d_t2_2_addr, const int* __restrict__ d_t2_2_offset, double* d_t2_3, const int* __restrict__ d_t2_3_addr, const int* __restrict__ d_t2_3_offset, double* d_t2_4, const int* __restrict__ d_t2_4_addr, const int* __restrict__ d_t2_4_offset, double* d_t2_5, const int* __restrict__ d_t2_5_addr, const int* __restrict__ d_t2_5_offset, double* d_t2_6, const int* __restrict__ d_t2_6_addr, const int* __restrict__ d_t2_6_offset, 
double* d_v2_1, const int* __restrict__ d_v2_1_addr, const int* __restrict__ d_v2_1_offset, double* d_v2_2, const int* __restrict__ d_v2_2_addr, const int* __restrict__ d_v2_2_offset, double* d_v2_3, const int* __restrict__ d_v2_3_addr, const int* __restrict__ d_v2_3_offset, double* d_v2_4, const int* __restrict__ d_v2_4_addr, const int* __restrict__ d_v2_4_offset, double* d_v2_5, const int* __restrict__ d_v2_5_addr, const int* __restrict__ d_v2_5_offset, double* d_v2_6, const int* __restrict__ d_v2_6_addr, const int* __restrict__ d_v2_6_offset, 
int kernel_1, int kernel_2, int kernel_3, 
int kernel_4, int kernel_5, int kernel_6, 
int kernel_7, int kernel_8, int kernel_9, 	// 7, 8, 9 are useless
int stride_reg_x, int stride_reg_y,
int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int l_idx_t3         = threadIdx.x + threadIdx.y * SIZE_TB_1_X;
	int t3_base_thread   = t3_output_base_1[blockIdx.x] + t3_output_offset_1[l_idx_t3];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_1[d_t2_1_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_1[d_v2_1_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_1_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_1_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_2[d_t2_2_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_2[d_v2_2_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_2_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_2_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_3[d_t2_3_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P4 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_3[d_v2_3_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P5)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_3_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_3_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_4[d_t2_4_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_4[d_v2_4_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H3 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_4_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_4_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_5[d_t2_5_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_5[d_v2_5_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H1 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_5_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_5_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_t2_6[d_t2_6_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_P5 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_1_Y] = d_v2_6[d_v2_6_addr[threadIdx.y + ll * SIZE_TB_1_Y + blockIdx.x * (SIZE_SLICE_1_H2 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_b[ll][d_v2_6_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_a[ll][d_t2_6_offset[l_idx_t3] + (xx * 16)];

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
			t3[t3_base_thread + (i * stride_reg_y) + (j * stride_reg_x)] += reg_tile[i][j];
		}
	}
}

// created by tc_gen_code_Kernel()
__global__ void kernel_ccsdT_2(double* t3, const int* __restrict__ t3_blk_rng_2, const int* __restrict__ t3_output_base_2, const int* __restrict__ t3_output_offset_2, 
double* d_t2_7, const int* __restrict__ d_t2_7_addr, const int* __restrict__ d_t2_7_offset, double* d_t2_8, const int* __restrict__ d_t2_8_addr, const int* __restrict__ d_t2_8_offset, double* d_t2_9, const int* __restrict__ d_t2_9_addr, const int* __restrict__ d_t2_9_offset, 
double* d_v2_7, const int* __restrict__ d_v2_7_addr, const int* __restrict__ d_v2_7_offset, double* d_v2_8, const int* __restrict__ d_v2_8_addr, const int* __restrict__ d_v2_8_offset, double* d_v2_9, const int* __restrict__ d_v2_9_addr, const int* __restrict__ d_v2_9_offset, 
int kernel_1, int kernel_2, int kernel_3, 	// 1, 2, 3, 4, 5, 6 are useless
int kernel_4, int kernel_5, int kernel_6, 
int kernel_7, int kernel_8, int kernel_9, 
int stride_reg_x, int stride_reg_y,
int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int l_idx_t3         = threadIdx.x + threadIdx.y * SIZE_TB_2_X;
	int t3_base_thread   = t3_output_base_2[blockIdx.x] + t3_output_offset_2[l_idx_t3];
	int internal_upperbound   = 0;
	int internal_offset;

	// should support for non-full tiles
	int idx_h3 = threadIdx.x % SIZE_SLICE_2_H3;
	int idx_h2 = threadIdx.x / SIZE_SLICE_2_H3;
	int idx_p4 = threadIdx.y % SIZE_SLICE_2_P4;
	int idx_h1 = threadIdx.y / SIZE_SLICE_2_P4;

	int rng_h3 = t3_blk_rng_2[blockIdx.x * NUM_INDEX + 0];
	int rng_h2 = t3_blk_rng_2[blockIdx.x * NUM_INDEX + 1];
	int rng_h1 = t3_blk_rng_2[blockIdx.x * NUM_INDEX + 2];
	int rng_p6 = t3_blk_rng_2[blockIdx.x * NUM_INDEX + 3];
	int rng_p5 = t3_blk_rng_2[blockIdx.x * NUM_INDEX + 4];
	int rng_p4 = t3_blk_rng_2[blockIdx.x * NUM_INDEX + 5];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_7[d_t2_7_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		if (idx_p4 < rng_h3 && idx_h1 < rng_p4 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_7[d_v2_7_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H3 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_7_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_8[d_t2_8_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H2 * SIZE_SLICE_2_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
        if (idx_p4 < rng_h1 && idx_h1 < rng_p4 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_8[d_v2_8_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H1 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_8_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_9[d_t2_9_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
        if (idx_p4 < rng_h2 && idx_h1 < rng_p4 && threadIdx.x < SIZE_INT_UNIT - internal_upperbound)
		for (int ll = 0; ll < rng_p5; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_9[d_v2_9_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H2 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_9_offset[l_idx_t3] + (xx * 16)];

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
__global__ void kernel_ccsdT_2_full(double* t3, const int* __restrict__ t3_blk_rng_2, const int* __restrict__ t3_output_base_2, const int* __restrict__ t3_output_offset_2, 
double* d_t2_7, const int* __restrict__ d_t2_7_addr, const int* __restrict__ d_t2_7_offset, double* d_t2_8, const int* __restrict__ d_t2_8_addr, const int* __restrict__ d_t2_8_offset, double* d_t2_9, const int* __restrict__ d_t2_9_addr, const int* __restrict__ d_t2_9_offset, 
double* d_v2_7, const int* __restrict__ d_v2_7_addr, const int* __restrict__ d_v2_7_offset, double* d_v2_8, const int* __restrict__ d_v2_8_addr, const int* __restrict__ d_v2_8_offset, double* d_v2_9, const int* __restrict__ d_v2_9_addr, const int* __restrict__ d_v2_9_offset, 
int kernel_1, int kernel_2, int kernel_3, 
int kernel_4, int kernel_5, int kernel_6, 
int kernel_7, int kernel_8, int kernel_9, 
int stride_reg_x, int stride_reg_y,
int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int l_idx_t3         = threadIdx.x + threadIdx.y * SIZE_TB_2_X;
	int t3_base_thread   = t3_output_base_2[blockIdx.x] + t3_output_offset_2[l_idx_t3];
	int internal_upperbound   = 0;
	int internal_offset;

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_7[d_t2_7_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_7[d_v2_7_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H3 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_7_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_8[d_t2_8_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H2 * SIZE_SLICE_2_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_8[d_v2_8_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H1 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_8_offset[l_idx_t3] + (xx * 16)];

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
		for (int ll = 0; ll < 4; ll++)
		{
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_9[d_t2_9_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_9[d_v2_9_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H2 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT - internal_upperbound; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_9_offset[l_idx_t3] + (xx * 16)];

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
__global__ void kernel_ccsdT_2_full_full(double* t3, const int* __restrict__ t3_blk_rng_2, const int* __restrict__ t3_output_base_2, const int* __restrict__ t3_output_offset_2, 
double* d_t2_7, const int* __restrict__ d_t2_7_addr, const int* __restrict__ d_t2_7_offset, double* d_t2_8, const int* __restrict__ d_t2_8_addr, const int* __restrict__ d_t2_8_offset, double* d_t2_9, const int* __restrict__ d_t2_9_addr, const int* __restrict__ d_t2_9_offset, 
double* d_v2_7, const int* __restrict__ d_v2_7_addr, const int* __restrict__ d_v2_7_offset, double* d_v2_8, const int* __restrict__ d_v2_8_addr, const int* __restrict__ d_v2_8_offset, double* d_v2_9, const int* __restrict__ d_v2_9_addr, const int* __restrict__ d_v2_9_offset, 
int kernel_1, int kernel_2, int kernel_3, 
int kernel_4, int kernel_5, int kernel_6, 
int kernel_7, int kernel_8, int kernel_9, 
int stride_reg_x, int stride_reg_y,
int size_internal)
{
	// For Shared Memory,
	__shared__ double sm_a[16][64 + 1];
	__shared__ double sm_b[16][64 + 1];

	int l_idx_t3         = threadIdx.x + threadIdx.y * SIZE_TB_2_X;
	int t3_base_thread   = t3_output_base_2[blockIdx.x] + t3_output_offset_2[l_idx_t3];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_7[d_t2_7_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H2)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_7[d_v2_7_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H3 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_7_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_7_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_8[d_t2_8_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H2 * SIZE_SLICE_2_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_8[d_v2_8_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H1 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		// Part: Generalized Threads
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_8_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_8_offset[l_idx_t3] + (xx * 16)];

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
			sm_a[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_t2_9[d_t2_9_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_P6 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_H3)] + (threadIdx.x + l)];
		}

		// Load Input Tensor to Shared Memory
		for (int ll = 0; ll < 4; ll++)
		{
			sm_b[threadIdx.x][threadIdx.y + ll * SIZE_TB_2_Y] = d_v2_9[d_v2_9_addr[threadIdx.y + ll * SIZE_TB_2_Y + blockIdx.x * (SIZE_SLICE_2_H2 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_P4)] + (threadIdx.x + l)];
		}
		__syncthreads();

		// Cross-Product: 16
		for (int ll = 0; ll < SIZE_INT_UNIT; ll++)
		{
			temp_bv[0] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 0];
			temp_bv[1] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 16];
			temp_bv[2] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 32];
			temp_bv[3] = sm_a[ll][d_t2_9_offset[l_idx_t3] + 48];

			for (int xx = 0 ; xx < 4; xx++)
			{
				temp_av = sm_b[ll][d_v2_9_offset[l_idx_t3] + (xx * 16)];

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

	// Device and Host Meory for Pre-Computed Arrays
	int *h_t3_block_index_1, 		*h_t3_block_range_1;
	int *d_t3_block_range_1;
	int *h_t3_block_index_2, 		*h_t3_block_range_2;
	int *d_t3_block_range_2;
	int *h_t3_output_base_addr_1, 	*h_t3_output_base_addr_2;
	int *d_t3_output_base_addr_1, 	*d_t3_output_base_addr_2;
	int *h_t3_output_ext_offset_1, 	*h_t3_output_ext_offset_2;
	int *d_t3_output_ext_offset_1, 	*d_t3_output_ext_offset_2;

	// Device Memory for Inputs and Output
	double *dev_t3;
	double *dev_t2_1, *dev_t2_2, *dev_t2_3, *dev_t2_4, *dev_t2_5, *dev_t2_6, *dev_t2_7, *dev_t2_8, *dev_t2_9;
	double *dev_v2_1, *dev_v2_2, *dev_v2_3, *dev_v2_4, *dev_v2_5, *dev_v2_6, *dev_v2_7, *dev_v2_8, *dev_v2_9;

	// Device and Host Memory for Pre-Computed Arrays
	int *dev_t2_1_addr, *dev_t2_2_addr, *dev_t2_3_addr, *dev_t2_4_addr, *dev_t2_5_addr, *dev_t2_6_addr, *dev_t2_7_addr, *dev_t2_8_addr, *dev_t2_9_addr;
	int *dev_v2_1_addr, *dev_v2_2_addr, *dev_v2_3_addr, *dev_v2_4_addr, *dev_v2_5_addr, *dev_v2_6_addr, *dev_v2_7_addr, *dev_v2_8_addr, *dev_v2_9_addr;
	int *host_t2_1_addr, *host_t2_2_addr, *host_t2_3_addr, *host_t2_4_addr, *host_t2_5_addr, *host_t2_6_addr, *host_t2_7_addr, *host_t2_8_addr, *host_t2_9_addr;
	int *host_v2_1_addr, *host_v2_2_addr, *host_v2_3_addr, *host_v2_4_addr, *host_v2_5_addr, *host_v2_6_addr, *host_v2_7_addr, *host_v2_8_addr, *host_v2_9_addr;

	int *dev_t2_1_offset,	*dev_t2_2_offset, 	*dev_t2_3_offset, 	*dev_t2_4_offset,	*dev_t2_5_offset, 	*dev_t2_6_offset, 	*dev_t2_7_offset, 	*dev_t2_8_offset, 	*dev_t2_9_offset;
	int *dev_v2_1_offset, 	*dev_v2_2_offset, 	*dev_v2_3_offset, 	*dev_v2_4_offset, 	*dev_v2_5_offset, 	*dev_v2_6_offset, 	*dev_v2_7_offset, 	*dev_v2_8_offset, 	*dev_v2_9_offset;
	int *host_t2_1_offset, 	*host_t2_2_offset, 	*host_t2_3_offset, 	*host_t2_4_offset, 	*host_t2_5_offset, 	*host_t2_6_offset, 	*host_t2_7_offset, 	*host_t2_8_offset, 	*host_t2_9_offset;
	int *host_v2_1_offset, 	*host_v2_2_offset, 	*host_v2_3_offset, 	*host_v2_4_offset, 	*host_v2_5_offset, 	*host_v2_6_offset, 	*host_v2_7_offset, 	*host_v2_8_offset, 	*host_v2_9_offset;

	// (1) Block-Ranges for Partial Tiles (h_t3_block_index is used in (2)) 
	Pre_TileApproach_1(h_t3_block_index_1, h_t3_block_range_1, &num_blocks_kernel_1, size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7);
	Pre_TileApproach_2(h_t3_block_index_2, h_t3_block_range_2, &num_blocks_kernel_2, size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7);

	// (2) Pre-Computed Arrays
	Pre_PreComputedArrays_1(h_t3_output_base_addr_1, 	host_t2_1_addr, host_v2_1_addr, 
														host_t2_2_addr, host_v2_2_addr,
														host_t2_3_addr, host_v2_3_addr,
														host_t2_4_addr, host_v2_4_addr,
														host_t2_5_addr, host_v2_5_addr,
														host_t2_6_addr, host_v2_6_addr,
							h_t3_output_ext_offset_1,	host_t2_1_offset, host_v2_1_offset,
														host_t2_2_offset, host_v2_2_offset,
														host_t2_3_offset, host_v2_3_offset,
														host_t2_4_offset, host_v2_4_offset,
														host_t2_5_offset, host_v2_5_offset,
														host_t2_6_offset, host_v2_6_offset,
							h_t3_block_index_1,
							num_blocks_kernel_1,
							size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7);
	
	Pre_PreComputedArrays_2(h_t3_output_base_addr_2, 	host_t2_7_addr, host_v2_7_addr,
														host_t2_8_addr, host_v2_8_addr,
														host_t2_9_addr, host_v2_9_addr,
							h_t3_output_ext_offset_2,	host_t2_7_offset, host_v2_7_offset,
														host_t2_8_offset, host_v2_8_offset,
														host_t2_9_offset, host_v2_9_offset,
							h_t3_block_index_2,
							num_blocks_kernel_2,
							size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7);
	
	// (3) Malloc CUDA 
	cudaMalloc((void**) &d_t3_block_range_1, 		sizeof(int) * (num_blocks_kernel_1) * NUM_INDEX);
	cudaMalloc((void**) &d_t3_block_range_2, 		sizeof(int) * (num_blocks_kernel_2) * NUM_INDEX);
	cudaMalloc((void**) &d_t3_output_base_addr_1, 	sizeof(int) * (num_blocks_kernel_1));
	cudaMalloc((void**) &d_t3_output_base_addr_2, 	sizeof(int) * (num_blocks_kernel_2));
	cudaMalloc((void**) &d_t3_output_ext_offset_1, 	sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y);
	cudaMalloc((void**) &d_t3_output_ext_offset_2, 	sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y);

	cudaMalloc((void**) &dev_t2_1_addr, 			sizeof(int) * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_P4 * num_blocks_kernel_1);
	cudaMalloc((void**) &dev_t2_2_addr, 			sizeof(int) * SIZE_SLICE_1_H3 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_P4 * num_blocks_kernel_1);
	cudaMalloc((void**) &dev_t2_3_addr, 			sizeof(int) * SIZE_SLICE_1_H3 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_P4 * num_blocks_kernel_1);
	cudaMalloc((void**) &dev_t2_4_addr, 			sizeof(int) * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_P5 * num_blocks_kernel_1);
	cudaMalloc((void**) &dev_t2_5_addr, 			sizeof(int) * SIZE_SLICE_1_H3 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_P5 * num_blocks_kernel_1);
	cudaMalloc((void**) &dev_t2_6_addr, 			sizeof(int) * SIZE_SLICE_1_H3 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_P5 * num_blocks_kernel_1);
	cudaMalloc((void**) &dev_v2_1_addr, 			sizeof(int) * SIZE_SLICE_1_P5 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H3 * num_blocks_kernel_1);
	cudaMalloc((void**) &dev_v2_2_addr, 			sizeof(int) * SIZE_SLICE_1_P5 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H1 * num_blocks_kernel_1);
	cudaMalloc((void**) &dev_v2_3_addr, 			sizeof(int) * SIZE_SLICE_1_P5 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H2 * num_blocks_kernel_1);
	cudaMalloc((void**) &dev_v2_4_addr, 			sizeof(int) * SIZE_SLICE_1_P4 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H3 * num_blocks_kernel_1);
	cudaMalloc((void**) &dev_v2_5_addr, 			sizeof(int) * SIZE_SLICE_1_P4 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H1 * num_blocks_kernel_1);
	cudaMalloc((void**) &dev_v2_6_addr, 			sizeof(int) * SIZE_SLICE_1_P4 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H2 * num_blocks_kernel_1);

	cudaMalloc((void**) &dev_t2_7_addr, 			sizeof(int) * SIZE_SLICE_2_H2 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_P6 * num_blocks_kernel_2);
	cudaMalloc((void**) &dev_t2_8_addr, 			sizeof(int) * SIZE_SLICE_2_H3 * SIZE_SLICE_2_H2 * SIZE_SLICE_2_P6 * num_blocks_kernel_2);
	cudaMalloc((void**) &dev_t2_9_addr, 			sizeof(int) * SIZE_SLICE_2_H3 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_P6 * num_blocks_kernel_2);
	cudaMalloc((void**) &dev_v2_7_addr, 			sizeof(int) * SIZE_SLICE_2_P4 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_H3 * num_blocks_kernel_2);
	cudaMalloc((void**) &dev_v2_8_addr, 			sizeof(int) * SIZE_SLICE_2_P4 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_H1 * num_blocks_kernel_2);
	cudaMalloc((void**) &dev_v2_9_addr, 			sizeof(int) * SIZE_SLICE_2_P4 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_H2 * num_blocks_kernel_2);

	cudaMalloc((void**) &dev_t2_1_offset, 			sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y);
	cudaMalloc((void**) &dev_t2_2_offset, 			sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y);
	cudaMalloc((void**) &dev_t2_3_offset, 			sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y);
	cudaMalloc((void**) &dev_t2_4_offset, 			sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y);
	cudaMalloc((void**) &dev_t2_5_offset, 			sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y);
	cudaMalloc((void**) &dev_t2_6_offset, 			sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y);
	cudaMalloc((void**) &dev_v2_1_offset, 			sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y);
	cudaMalloc((void**) &dev_v2_2_offset, 			sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y);
	cudaMalloc((void**) &dev_v2_3_offset, 			sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y);
	cudaMalloc((void**) &dev_v2_4_offset, 			sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y);
	cudaMalloc((void**) &dev_v2_5_offset, 			sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y);
	cudaMalloc((void**) &dev_v2_6_offset, 			sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y);

	cudaMalloc((void**) &dev_t2_7_offset, 			sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y);
	cudaMalloc((void**) &dev_t2_8_offset, 			sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y);
	cudaMalloc((void**) &dev_t2_9_offset, 			sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y);
	cudaMalloc((void**) &dev_v2_7_offset, 			sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y);
	cudaMalloc((void**) &dev_v2_8_offset, 			sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y);
	cudaMalloc((void**) &dev_v2_9_offset, 			sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y);

	// (4) Memcpy CUDA
	cudaMemcpy(d_t3_block_range_1, 			h_t3_block_range_1, 		sizeof(int) * num_blocks_kernel_1 * NUM_INDEX, 	cudaMemcpyHostToDevice);
	cudaMemcpy(d_t3_block_range_2, 			h_t3_block_range_2, 		sizeof(int) * num_blocks_kernel_2 * NUM_INDEX, 	cudaMemcpyHostToDevice);
	cudaMemcpy(d_t3_output_base_addr_1, 	h_t3_output_base_addr_1, 	sizeof(int) * num_blocks_kernel_1, 				cudaMemcpyHostToDevice);
	cudaMemcpy(d_t3_output_base_addr_2, 	h_t3_output_base_addr_2, 	sizeof(int) * num_blocks_kernel_2, 				cudaMemcpyHostToDevice);
	cudaMemcpy(d_t3_output_ext_offset_1, 	h_t3_output_ext_offset_1, 	sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y, 		cudaMemcpyHostToDevice);
	cudaMemcpy(d_t3_output_ext_offset_2, 	h_t3_output_ext_offset_2, 	sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y, 		cudaMemcpyHostToDevice);

	cudaMemcpy(dev_t2_1_addr, host_t2_1_addr, sizeof(int) * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_P4 * num_blocks_kernel_1, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_1_addr, host_v2_1_addr, sizeof(int) * SIZE_SLICE_1_P5 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H3 * num_blocks_kernel_1, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_2_addr, host_t2_2_addr, sizeof(int) * SIZE_SLICE_1_H3 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_P4 * num_blocks_kernel_1, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_2_addr, host_v2_2_addr, sizeof(int) * SIZE_SLICE_1_P5 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H1 * num_blocks_kernel_1, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_3_addr, host_t2_3_addr, sizeof(int) * SIZE_SLICE_1_H3 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_P4 * num_blocks_kernel_1, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_3_addr, host_v2_3_addr, sizeof(int) * SIZE_SLICE_1_P5 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H2 * num_blocks_kernel_1, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_4_addr, host_t2_4_addr, sizeof(int) * SIZE_SLICE_1_H2 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_P5 * num_blocks_kernel_1, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_4_addr, host_v2_4_addr, sizeof(int) * SIZE_SLICE_1_P4 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H3 * num_blocks_kernel_1, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_5_addr, host_t2_5_addr, sizeof(int) * SIZE_SLICE_1_H3 * SIZE_SLICE_1_H2 * SIZE_SLICE_1_P5 * num_blocks_kernel_1, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_5_addr, host_v2_5_addr, sizeof(int) * SIZE_SLICE_1_P4 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H1 * num_blocks_kernel_1, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_6_addr, host_t2_6_addr, sizeof(int) * SIZE_SLICE_1_H3 * SIZE_SLICE_1_H1 * SIZE_SLICE_1_P5 * num_blocks_kernel_1, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_6_addr, host_v2_6_addr, sizeof(int) * SIZE_SLICE_1_P4 * SIZE_SLICE_1_P6 * SIZE_SLICE_1_H2 * num_blocks_kernel_1, cudaMemcpyHostToDevice);

	cudaMemcpy(dev_t2_7_addr, host_t2_7_addr, sizeof(int) * SIZE_SLICE_2_H2 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_P6 * num_blocks_kernel_2, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_7_addr, host_v2_7_addr, sizeof(int) * SIZE_SLICE_2_P4 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_H3 * num_blocks_kernel_2, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_8_addr, host_t2_8_addr, sizeof(int) * SIZE_SLICE_2_H3 * SIZE_SLICE_2_H2 * SIZE_SLICE_2_P6 * num_blocks_kernel_2, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_8_addr, host_v2_8_addr, sizeof(int) * SIZE_SLICE_2_P4 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_H1 * num_blocks_kernel_2, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_9_addr, host_t2_9_addr, sizeof(int) * SIZE_SLICE_2_H3 * SIZE_SLICE_2_H1 * SIZE_SLICE_2_P6 * num_blocks_kernel_2, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_9_addr, host_v2_9_addr, sizeof(int) * SIZE_SLICE_2_P4 * SIZE_SLICE_2_P5 * SIZE_SLICE_2_H2 * num_blocks_kernel_2, cudaMemcpyHostToDevice);

	cudaMemcpy(dev_t2_1_offset, host_t2_1_offset, sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_1_offset, host_v2_1_offset, sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_2_offset, host_t2_2_offset, sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_2_offset, host_v2_2_offset, sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_3_offset, host_t2_3_offset, sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_3_offset, host_v2_3_offset, sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_4_offset, host_t2_4_offset, sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_4_offset, host_v2_4_offset, sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_5_offset, host_t2_5_offset, sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_5_offset, host_v2_5_offset, sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_6_offset, host_t2_6_offset, sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_6_offset, host_v2_6_offset, sizeof(int) * SIZE_TB_1_X * SIZE_TB_1_Y, cudaMemcpyHostToDevice);

	cudaMemcpy(dev_t2_7_offset, host_t2_7_offset, sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_7_offset, host_v2_7_offset, sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_8_offset, host_t2_8_offset, sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_8_offset, host_v2_8_offset, sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_t2_9_offset, host_t2_9_offset, sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v2_9_offset, host_v2_9_offset, sizeof(int) * SIZE_TB_2_X * SIZE_TB_2_Y, cudaMemcpyHostToDevice);

	// (Temporary)
	//cudaMalloc((void**) &dev_t3, 	sizeof(double) * size_h3 * size_h2 * size_h1 * size_p6 * size_p5 * size_p4);
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

    dev_t3 = t3_d;
	//cudaMemcpy(dev_t3, 	 t3,   sizeof(double) * size_h3 * size_h2 * size_h1 * size_p6 * size_p5 * size_p4, cudaMemcpyHostToDevice);
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

	// (5) launch kernel(s)
	/*
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
	*/
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


	// External Index (Full) : A Multiple of "4"
	if (size_h3 % SIZE_SLICE_1_H3 == 0 && size_h2 % SIZE_SLICE_1_H2 == 0 && size_h1 % SIZE_SLICE_1_H1 == 0 && 
		size_p6 % SIZE_SLICE_1_P6 == 0 && size_p5 % SIZE_SLICE_1_P5 == 0 && size_p4 % SIZE_SLICE_1_P4 == 0)
	{
		// Internal Index (Full) : A Multiple of "16"
		if (size_p7 % SIZE_INT_UNIT == 0)
		{
			if (opt_register_transpose == 1)
			{
				kernel_ccsdT_0_full_full<<<gridsize_1, blocksize_1>>>(dev_t3, d_t3_block_range_1, d_t3_output_base_addr_1, d_t3_output_ext_offset_1, 
																	dev_t2_1, dev_t2_1_addr, dev_t2_1_offset, dev_t2_2, dev_t2_2_addr, dev_t2_2_offset, dev_t2_3, dev_t2_3_addr, dev_t2_3_offset, dev_t2_4, dev_t2_4_addr, dev_t2_4_offset, dev_t2_5, dev_t2_5_addr, dev_t2_5_offset, dev_t2_6, dev_t2_6_addr, dev_t2_6_offset, 
																	dev_v2_1, dev_v2_1_addr, dev_v2_1_offset, dev_v2_2, dev_v2_2_addr, dev_v2_2_offset, dev_v2_3, dev_v2_3_addr, dev_v2_3_offset, dev_v2_4, dev_v2_4_addr, dev_v2_4_offset, dev_v2_5, dev_v2_5_addr, dev_v2_5_offset, dev_v2_6, dev_v2_6_addr, dev_v2_6_offset, 
																	d_t3_block_range_2, d_t3_output_base_addr_2, d_t3_output_ext_offset_2, 
																	dev_t2_7, dev_t2_7_addr, dev_t2_7_offset, dev_t2_8, dev_t2_8_addr, dev_t2_8_offset, dev_t2_9, dev_t2_9_addr, dev_t2_9_offset, 
																	dev_v2_7, dev_v2_7_addr, dev_v2_7_offset, dev_v2_8, dev_v2_8_addr, dev_v2_8_offset, dev_v2_9, dev_v2_9_addr, dev_v2_9_offset, 
																	kernel_1, kernel_2, kernel_3, 
																	kernel_4, kernel_5, kernel_6, 
																	kernel_7, kernel_8, kernel_9, 
																	str_reg_x_1, str_reg_y_1,
																	size_internal);	
			}
			else
			{
				// Depends on # of Fused Kernel
				kernel_ccsdT_1_full_full<<<gridsize_1, blocksize_1>>>(dev_t3, d_t3_block_range_1, d_t3_output_base_addr_1, d_t3_output_ext_offset_1, 
																	dev_t2_1, dev_t2_1_addr, dev_t2_1_offset, dev_t2_2, dev_t2_2_addr, dev_t2_2_offset, dev_t2_3, dev_t2_3_addr, dev_t2_3_offset, dev_t2_4, dev_t2_4_addr, dev_t2_4_offset, dev_t2_5, dev_t2_5_addr, dev_t2_5_offset, dev_t2_6, dev_t2_6_addr, dev_t2_6_offset, 
																	dev_v2_1, dev_v2_1_addr, dev_v2_1_offset, dev_v2_2, dev_v2_2_addr, dev_v2_2_offset, dev_v2_3, dev_v2_3_addr, dev_v2_3_offset, dev_v2_4, dev_v2_4_addr, dev_v2_4_offset, dev_v2_5, dev_v2_5_addr, dev_v2_5_offset, dev_v2_6, dev_v2_6_addr, dev_v2_6_offset, 
																	kernel_1, kernel_2, kernel_3, 
																	kernel_4, kernel_5, kernel_6, 
																	kernel_7, kernel_8, kernel_9, 
																	str_reg_x_1, str_reg_y_1,
																	size_internal);

				kernel_ccsdT_2_full_full<<<gridsize_2, blocksize_2>>>(dev_t3, d_t3_block_range_2, d_t3_output_base_addr_2, d_t3_output_ext_offset_2, 
																	dev_t2_7, dev_t2_7_addr, dev_t2_7_offset, dev_t2_8, dev_t2_8_addr, dev_t2_8_offset, dev_t2_9, dev_t2_9_addr, dev_t2_9_offset, 
																	dev_v2_7, dev_v2_7_addr, dev_v2_7_offset, dev_v2_8, dev_v2_8_addr, dev_v2_8_offset, dev_v2_9, dev_v2_9_addr, dev_v2_9_offset, 
																	kernel_1, kernel_2, kernel_3, 
																	kernel_4, kernel_5, kernel_6, 
																	kernel_7, kernel_8, kernel_9, 
																	str_reg_x_2, str_reg_y_2,
																	size_internal);
			}
		}
		else
		{
			if (opt_register_transpose == 1)
			{
				kernel_ccsdT_0_full<<<gridsize_1, blocksize_1>>>(dev_t3, d_t3_block_range_1, d_t3_output_base_addr_1, d_t3_output_ext_offset_1, 
																dev_t2_1, dev_t2_1_addr, dev_t2_1_offset, dev_t2_2, dev_t2_2_addr, dev_t2_2_offset, dev_t2_3, dev_t2_3_addr, dev_t2_3_offset, 
																dev_t2_4, dev_t2_4_addr, dev_t2_4_offset, dev_t2_5, dev_t2_5_addr, dev_t2_5_offset, dev_t2_6, dev_t2_6_addr, dev_t2_6_offset, 
																dev_v2_1, dev_v2_1_addr, dev_v2_1_offset, dev_v2_2, dev_v2_2_addr, dev_v2_2_offset, dev_v2_3, dev_v2_3_addr, dev_v2_3_offset, 
																dev_v2_4, dev_v2_4_addr, dev_v2_4_offset, dev_v2_5, dev_v2_5_addr, dev_v2_5_offset, dev_v2_6, dev_v2_6_addr, dev_v2_6_offset, 
																d_t3_block_range_2, d_t3_output_base_addr_2, d_t3_output_ext_offset_2, 
																dev_t2_7, dev_t2_7_addr, dev_t2_7_offset, dev_t2_8, dev_t2_8_addr, dev_t2_8_offset, dev_t2_9, dev_t2_9_addr, dev_t2_9_offset, 
																dev_v2_7, dev_v2_7_addr, dev_v2_7_offset, dev_v2_8, dev_v2_8_addr, dev_v2_8_offset, dev_v2_9, dev_v2_9_addr, dev_v2_9_offset, 
																kernel_1, kernel_2, kernel_3, 
																kernel_4, kernel_5, kernel_6, 
																kernel_7, kernel_8, kernel_9, 
																str_reg_x_1, str_reg_y_1,
																size_internal);	
			}
			else
			{
				kernel_ccsdT_1_full<<<gridsize_1, blocksize_1>>>(dev_t3, d_t3_block_range_1, d_t3_output_base_addr_1, d_t3_output_ext_offset_1, 
																dev_t2_1, dev_t2_1_addr, dev_t2_1_offset, dev_t2_2, dev_t2_2_addr, dev_t2_2_offset, dev_t2_3, dev_t2_3_addr, dev_t2_3_offset, 
																dev_t2_4, dev_t2_4_addr, dev_t2_4_offset, dev_t2_5, dev_t2_5_addr, dev_t2_5_offset, dev_t2_6, dev_t2_6_addr, dev_t2_6_offset, 
																dev_v2_1, dev_v2_1_addr, dev_v2_1_offset, dev_v2_2, dev_v2_2_addr, dev_v2_2_offset, dev_v2_3, dev_v2_3_addr, dev_v2_3_offset, 
																dev_v2_4, dev_v2_4_addr, dev_v2_4_offset, dev_v2_5, dev_v2_5_addr, dev_v2_5_offset, dev_v2_6, dev_v2_6_addr, dev_v2_6_offset, 
																kernel_1, kernel_2, kernel_3, 
																kernel_4, kernel_5, kernel_6, 
																kernel_7, kernel_8, kernel_9, 
																str_reg_x_1, str_reg_y_1,
																size_internal);

				kernel_ccsdT_2_full<<<gridsize_2, blocksize_2>>>(dev_t3, d_t3_block_range_2, d_t3_output_base_addr_2, d_t3_output_ext_offset_2, 
																dev_t2_7, dev_t2_7_addr, dev_t2_7_offset, dev_t2_8, dev_t2_8_addr, dev_t2_8_offset, dev_t2_9, dev_t2_9_addr, dev_t2_9_offset, 
																dev_v2_7, dev_v2_7_addr, dev_v2_7_offset, dev_v2_8, dev_v2_8_addr, dev_v2_8_offset, dev_v2_9, dev_v2_9_addr, dev_v2_9_offset, 
																kernel_1, kernel_2, kernel_3, 
																kernel_4, kernel_5, kernel_6, 
																kernel_7, kernel_8, kernel_9, 
																str_reg_x_2, str_reg_y_2,
																size_internal);
			}
		}
	}
	else
	{
		if (opt_register_transpose == 1)
		{
			kernel_ccsdT_0<<<gridsize_1, blocksize_1>>>(dev_t3, d_t3_block_range_1, d_t3_output_base_addr_1, d_t3_output_ext_offset_1, 
			dev_t2_1, dev_t2_1_addr, dev_t2_1_offset, dev_t2_2, dev_t2_2_addr, dev_t2_2_offset, dev_t2_3, dev_t2_3_addr, dev_t2_3_offset, 
			dev_t2_4, dev_t2_4_addr, dev_t2_4_offset, dev_t2_5, dev_t2_5_addr, dev_t2_5_offset, dev_t2_6, dev_t2_6_addr, dev_t2_6_offset, 
			dev_v2_1, dev_v2_1_addr, dev_v2_1_offset, dev_v2_2, dev_v2_2_addr, dev_v2_2_offset, dev_v2_3, dev_v2_3_addr, dev_v2_3_offset, 
			dev_v2_4, dev_v2_4_addr, dev_v2_4_offset, dev_v2_5, dev_v2_5_addr, dev_v2_5_offset, dev_v2_6, dev_v2_6_addr, dev_v2_6_offset, 
			d_t3_block_range_2, d_t3_output_base_addr_2, d_t3_output_ext_offset_2, 
			dev_t2_7, dev_t2_7_addr, dev_t2_7_offset, dev_t2_8, dev_t2_8_addr, dev_t2_8_offset, dev_t2_9, dev_t2_9_addr, dev_t2_9_offset, 
			dev_v2_7, dev_v2_7_addr, dev_v2_7_offset, dev_v2_8, dev_v2_8_addr, dev_v2_8_offset, dev_v2_9, dev_v2_9_addr, dev_v2_9_offset, 
			kernel_1, kernel_2, kernel_3, 
			kernel_4, kernel_5, kernel_6, 
			kernel_7, kernel_8, kernel_9, 
			str_reg_x_1, str_reg_y_1,
			size_internal);	
		}
		else
		{
			kernel_ccsdT_1<<<gridsize_1, blocksize_1>>>(dev_t3, d_t3_block_range_1, d_t3_output_base_addr_1, d_t3_output_ext_offset_1, 
			dev_t2_1, dev_t2_1_addr, dev_t2_1_offset, dev_t2_2, dev_t2_2_addr, dev_t2_2_offset, dev_t2_3, dev_t2_3_addr, dev_t2_3_offset, 
			dev_t2_4, dev_t2_4_addr, dev_t2_4_offset, dev_t2_5, dev_t2_5_addr, dev_t2_5_offset, dev_t2_6, dev_t2_6_addr, dev_t2_6_offset, 
			dev_v2_1, dev_v2_1_addr, dev_v2_1_offset, dev_v2_2, dev_v2_2_addr, dev_v2_2_offset, dev_v2_3, dev_v2_3_addr, dev_v2_3_offset, 
			dev_v2_4, dev_v2_4_addr, dev_v2_4_offset, dev_v2_5, dev_v2_5_addr, dev_v2_5_offset, dev_v2_6, dev_v2_6_addr, dev_v2_6_offset, 
			kernel_1, kernel_2, kernel_3, 
			kernel_4, kernel_5, kernel_6, 
			kernel_7, kernel_8, kernel_9, 
			str_reg_x_1, str_reg_y_1,
			size_internal);

			kernel_ccsdT_2<<<gridsize_2, blocksize_2>>>(dev_t3, d_t3_block_range_2, d_t3_output_base_addr_2, d_t3_output_ext_offset_2, 
			dev_t2_7, dev_t2_7_addr, dev_t2_7_offset, dev_t2_8, dev_t2_8_addr, dev_t2_8_offset, dev_t2_9, dev_t2_9_addr, dev_t2_9_offset, 
			dev_v2_7, dev_v2_7_addr, dev_v2_7_offset, dev_v2_8, dev_v2_8_addr, dev_v2_8_offset, dev_v2_9, dev_v2_9_addr, dev_v2_9_offset, 
			kernel_1, kernel_2, kernel_3, 
			kernel_4, kernel_5, kernel_6, 
			kernel_7, kernel_8, kernel_9, 
			str_reg_x_2, str_reg_y_2,
			size_internal);
		}
	}

	// (Temporary)
	//cudaMemcpy(t3, dev_t3, sizeof(double) * size_h3 * size_h2 * size_h1 * size_p6 * size_p5 * size_p4, cudaMemcpyDeviceToHost);

	// (Temporary)
	//cudaFree(dev_t3);
	cudaFree(dev_t2_1);	cudaFree(dev_t2_2);	cudaFree(dev_t2_3);	cudaFree(dev_t2_4);	cudaFree(dev_t2_5);	cudaFree(dev_t2_6);	cudaFree(dev_t2_7);	cudaFree(dev_t2_8);	cudaFree(dev_t2_9);	
	cudaFree(dev_v2_1);	cudaFree(dev_v2_2);	cudaFree(dev_v2_3);	cudaFree(dev_v2_4);	cudaFree(dev_v2_5);	cudaFree(dev_v2_6);	cudaFree(dev_v2_7);	cudaFree(dev_v2_8);	cudaFree(dev_v2_9);

	// (6) Free Pre-Computed Arrays (in Device Memory)
	cudaFree(d_t3_block_range_1);		cudaFree(d_t3_block_range_2);
	cudaFree(d_t3_output_base_addr_1);	cudaFree(d_t3_output_base_addr_2);
	cudaFree(d_t3_output_ext_offset_1);	cudaFree(d_t3_output_ext_offset_2);
	cudaFree(dev_t2_1_addr);	cudaFree(dev_t2_2_addr);	cudaFree(dev_t2_3_addr);	cudaFree(dev_t2_4_addr);	cudaFree(dev_t2_5_addr);	cudaFree(dev_t2_6_addr);	cudaFree(dev_t2_7_addr);	cudaFree(dev_t2_8_addr);	cudaFree(dev_t2_9_addr);
	cudaFree(dev_v2_1_addr);	cudaFree(dev_v2_2_addr);	cudaFree(dev_v2_3_addr);	cudaFree(dev_v2_4_addr);	cudaFree(dev_v2_5_addr);	cudaFree(dev_v2_6_addr);	cudaFree(dev_v2_7_addr);	cudaFree(dev_v2_8_addr);	cudaFree(dev_v2_9_addr);
	cudaFree(dev_t2_1_offset);	cudaFree(dev_t2_2_offset);	cudaFree(dev_t2_3_offset);	cudaFree(dev_t2_4_offset);	cudaFree(dev_t2_5_offset);	cudaFree(dev_t2_6_offset);	cudaFree(dev_t2_7_offset);	cudaFree(dev_t2_8_offset);	cudaFree(dev_t2_9_offset);	
	cudaFree(dev_v2_1_offset);	cudaFree(dev_v2_2_offset);	cudaFree(dev_v2_3_offset);	cudaFree(dev_v2_4_offset);	cudaFree(dev_v2_5_offset);	cudaFree(dev_v2_6_offset);	cudaFree(dev_v2_7_offset);	cudaFree(dev_v2_8_offset);	cudaFree(dev_v2_9_offset);


	// (7) Free Pre-Computed Arrays (in Host Memory)
	free(h_t3_block_index_1);		free(h_t3_block_index_2);
	free(h_t3_block_range_1);		free(h_t3_block_range_2);
	free(h_t3_output_base_addr_1);	free(h_t3_output_base_addr_2);
    free(h_t3_output_ext_offset_1);	free(h_t3_output_ext_offset_2);
	free(host_t2_1_addr);	free(host_t2_2_addr); 	free(host_t2_3_addr);	free(host_t2_4_addr);	free(host_t2_5_addr);	free(host_t2_6_addr);	free(host_t2_7_addr);	free(host_t2_8_addr);	free(host_t2_9_addr);
	free(host_v2_1_addr);	free(host_v2_2_addr);	free(host_v2_3_addr);	free(host_v2_4_addr);	free(host_v2_5_addr);	free(host_v2_6_addr);	free(host_v2_7_addr);	free(host_v2_8_addr);	free(host_v2_9_addr);
	free(host_t2_1_offset);	free(host_t2_2_offset); free(host_t2_3_offset);	free(host_t2_4_offset);	free(host_t2_5_offset);	free(host_t2_6_offset);	free(host_t2_7_offset);	free(host_t2_8_offset);	free(host_t2_9_offset);
	free(host_v2_1_offset);	free(host_v2_2_offset);	free(host_v2_3_offset);	free(host_v2_4_offset);	free(host_v2_5_offset);	free(host_v2_6_offset);	free(host_v2_7_offset);	free(host_v2_8_offset);	free(host_v2_9_offset);
}

//
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
	sd_t_d2_all_cuda(size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7, 
					t3, 
					t2_1, v2_1,	t2_2, v2_2, t2_3, v2_3,
					t2_4, v2_4,	t2_5, v2_5,	t2_6, v2_6,
					t2_7, v2_7,	t2_8, v2_8,	t2_9, v2_9, 
					kernel_1, kernel_2, kernel_3,
					kernel_4, kernel_5, kernel_6,
					kernel_7, kernel_8, kernel_9,
					opt_register_transpose);
}

//
//	>>> Interface <<<
//	Index' Sizes: h3, h2, h1, p6, p5, p4 and p7.
//	Tensors: t3 (output), t2_all (all inputs) and v2_all (all inputs).
//	Conditions: sd2_1, sd2_2, sd2_3, sd2_4, sd2_5, sd2_6, sd2_7, sd2_8 and sd2_9.
// 	Option: Register-Transpose 	(ON: a single kernel for all 9 functions), 
//								(OFF: two kernels for 6 and 3 functions, respectivley).
//
extern "C"
void sd_t_d2_all_cuda__(int 	size_h3, 	int 	size_h2, 	int 	size_h1, 
						int 	size_p6, 	int 	size_p5, 	int 	size_p4, 	int 	size_p7, 
						double* t3,   		double* t2_all,		double* v2_all,
						int 	kernel_1, 	int 	kernel_2, 	int 	kernel_3, 
						int 	kernel_4, 	int 	kernel_5, 	int 	kernel_6, 
						int 	kernel_7, 	int 	kernel_8, 	int 	kernel_9, 
						int 	opt_register_transpose)
{
	unsigned int size_T2_1, size_T2_2, size_T2_3, size_T2_4, size_T2_5, size_T2_6, size_T2_7, size_T2_8;//, size_T2_9;
	unsigned int size_V2_1, size_V2_2, size_V2_3, size_V2_4, size_V2_5, size_V2_6, size_V2_7, size_V2_8;//, size_V2_9;

	double* t2_1;	double* v2_1;
	double* t2_2;	double* v2_2;
	double* t2_3;	double* v2_3;
	double* t2_4;	double* v2_4;
	double* t2_5;	double* v2_5;
	double* t2_6;	double* v2_6;
	double* t2_7;	double* v2_7;
	double* t2_8;	double* v2_8;
	double* t2_9;	double* v2_9;

	// Each Input Tensor's Size
	size_T2_1 		= size_p7 * size_p4 * size_h1 * size_h2;
	size_T2_2 		= size_p7 * size_p4 * size_h2 * size_h3;
	size_T2_3 		= size_p7 * size_p4 * size_h1 * size_h3;
	size_T2_4 		= size_p7 * size_p5 * size_h1 * size_h2;
	size_T2_5 		= size_p7 * size_p5 * size_h2 * size_h3;
	size_T2_6 		= size_p7 * size_p5 * size_h1 * size_h3;
	size_T2_7 		= size_p7 * size_p6 * size_h1 * size_h2;
	size_T2_8 		= size_p7 * size_p6 * size_h2 * size_h3;

	size_V2_1 		= size_p7 * size_h3 * size_p6 * size_p5;
	size_V2_2 		= size_p7 * size_h1 * size_p6 * size_p5;
	size_V2_3 		= size_p7 * size_h2 * size_p6 * size_p5;
	size_V2_4 		= size_p7 * size_h3 * size_p6 * size_p4;
	size_V2_5 		= size_p7 * size_h1 * size_p6 * size_p4;
	size_V2_6 		= size_p7 * size_h2 * size_p6 * size_p4;
    size_V2_7 		= size_p7 * size_h3 * size_p5 * size_p4;
    size_V2_8 		= size_p7 * size_h1 * size_p5 * size_p4;

	//
	sd_t_d2_all_cuda(size_h3, size_h2, size_h1, size_p6, size_p5, size_p4, size_p7,
					t3, 
					t2_1, v2_1, t2_2, v2_2, t2_3, v2_3, 
					t2_4, v2_4, t2_5, v2_5, t2_6, v2_6,	
					t2_7, v2_7, t2_8, v2_8, t2_9, v2_9, 
					kernel_1, kernel_2, kernel_3,
					kernel_4, kernel_5, kernel_6,
					kernel_7, kernel_8, kernel_9,
					opt_register_transpose);

}
