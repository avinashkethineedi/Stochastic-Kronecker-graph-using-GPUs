#pragma once
#include<math.h>
#include<algorithm>
typedef struct block_probablity
{
	float a, b, c, d;
	long nnz;
}block;
typedef struct EDGE
{
	//edge from v to u with weight w
	int v;
	long u;
	float w;
}edge;
typedef struct CSR_MATRIX
{
	int *row_ptr;
	long *col_ptr;
	float *val_ptr;
	long nnz, rows;
}csr_data;
typedef struct TIME
{
	double t, avg_t, min_t, max_t;
}time_stats;