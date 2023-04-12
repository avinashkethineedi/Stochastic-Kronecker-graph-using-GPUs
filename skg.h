#pragma once
#include<math.h>
#include<algorithm>
typedef struct block_probablity
{
	float a, b, c, d;
	long edges;
}block;
typedef struct EDGE
{
	//edge from v to u with weight w
	long v;
	long u;
	float w;
}edge;
typedef struct CSR_MATRIX
{
	long *row_ptr;
	long *col_ptr;
	float *val_ptr;
	long edges, nodes;
}csr_data;
typedef struct TIME
{
	double t, avg_t, min_t, max_t;
}time_stats;
long** block_allocation(int n, int m)
{
	size_t size = (size_t)n*m;
	long **pe_blocks = (long**)malloc(sizeof(long*)*n);
	pe_blocks[0] = (long*)malloc(sizeof(long)*size);
	for(int i=1;i<n;i++)
		pe_blocks[i] = pe_blocks[0]+(size_t)i*m;
	return pe_blocks;
}
long* calculate_edge_distribution(int rank, int npes, block *mat_prop)
{
	long **pe_blocks = block_allocation(npes, npes);
	int i, j, k, l, m, x, y;
	int n=log2(npes);
	for(i=0;i<npes;i++)
		for(j=0;j<npes;j++)
			pe_blocks[i][j] = mat_prop->edges;
	int num_blocks, grid_size, block_row, block_col, grid_row, grid_col;
	float prob;
	for(i=0;i<n;i++)
	{
		num_blocks = 1<<i;
		grid_size = 1<<(n-1-i);
		for(j=0;j<num_blocks;j++)
		{
			block_row = j*grid_size*2;
			for(k=0;k<num_blocks;k++)
			{
				block_col = k*grid_size*2;
				for(l=0;l<2;l++)
				{
					grid_row = block_row+grid_size*l;
					for(m=0;m<2;m++)
					{
						grid_col = block_col+grid_size*m;
						prob = l?(m?mat_prop->d:mat_prop->c):(m?mat_prop->b:mat_prop->a);
						for(x=0;x<grid_size;x++)
							for(y=0;y<grid_size;y++)
								pe_blocks[grid_row+x][grid_col+y] = (long)round(pe_blocks[grid_row+x][grid_col+y]*prob);
					}
				}
			}
		}
	}
	long *edges_dist = (long*)malloc(sizeof(long)*npes);
	memcpy(edges_dist, pe_blocks[rank], (size_t)npes*sizeof(long));
	free(pe_blocks[0]);
	free(pe_blocks);
	return edges_dist;
}
long calculate_edges(long *edges_dist, int n)
{
	int i;
	long edges=0;
	for(i=0;i<n;i++)
		edges += edges_dist[i];
	return edges;
}
void stochastic_Kronecker_grpah(edge *edge_list, long pe_edges, long dim, long start_idx, float *prob, float *c_prob, int *node_edge_count)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	srand(start_idx+rank);
	long mat_dim, row, col, edge_count=0;
	int p_idx, prow, pcol, i, k;
	float p;
	k = log2(dim);
	for(edge_count=0;edge_count<pe_edges;edge_count++)
	{
		mat_dim = dim;
		row=0;
		col=0;
		for(i=0;i<k;i++)
		{
			p = (float)rand()/RAND_MAX;
			p_idx = p<c_prob[0]?0:(p<c_prob[1]?1:(p<c_prob[2]?2:3));
			prow = p_idx/2;
			pcol = p_idx%2;
			mat_dim/=2;
			row = row + mat_dim*prow;
			col = col + mat_dim*pcol;
		}
		edge_list[edge_count].v = row;
		edge_list[edge_count].u = start_idx+col;
		edge_list[edge_count].w = p;
		node_edge_count[row]++;
	}
}
edge* create_edge_list(long *edges_dist, long pe_edges, long nodes_per_pe, int npes, block *mat_prob, int **node_edge_count)
{
	edge *edge_list = (edge*)malloc((size_t)pe_edges*sizeof(edge));
	int i;
	*node_edge_count = (int*)calloc(nodes_per_pe, sizeof(int));
	long block_edges, start_idx, idx = 0;
	float prob[] = {mat_prob->a, mat_prob->b, mat_prob->c, mat_prob->d}, c_prob[4];
	c_prob[0] = prob[0];
	for(i=1;i<4;i++) c_prob[i] = prob[i] + c_prob[i-1];
	for(i=0;i<npes;i++)
	{
		block_edges = edges_dist[npes-1-i];
		start_idx = nodes_per_pe*(npes-1-i);
		stochastic_Kronecker_grpah(&edge_list[idx], block_edges, nodes_per_pe, start_idx, prob, c_prob, *node_edge_count);
		idx += block_edges;
	}
	return edge_list;
}
csr_data* create_csr_data(edge *edge_list, int *node_edge_count, long pe_edges, long nodes_per_pe)
{
	size_t n = (size_t)sizeof(csr_data)+(nodes_per_pe+1)*sizeof(long) + pe_edges*(sizeof(long) + sizeof(float));
	csr_data *csr_mat = (csr_data*)malloc(n);
	csr_mat->edges = pe_edges;
	csr_mat->nodes = nodes_per_pe;
	long i, idx;
	if(csr_mat == NULL)
	{
		printf("csr memory allocation failed (%ld B)\n", n);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
	//printf("edges: %ld\tnodes: %ld\n", csr_mat->edges, csr_mat->nodes);
	n = sizeof(csr_data);
	csr_mat->row_ptr = (long*)((char*)csr_mat + n);
	n+=((nodes_per_pe+1)*sizeof(long));
	csr_mat->col_ptr = (long*)((char*)csr_mat + n);
	n+=(pe_edges*sizeof(long));
	csr_mat->val_ptr = (float*)((char*)csr_mat + n);
	//row_ptr
	csr_mat->row_ptr[0] = 0;
	for(i=1;i<=nodes_per_pe;i++)
		csr_mat->row_ptr[i] = csr_mat->row_ptr[i-1] + node_edge_count[i-1];
	for(i=0;i<pe_edges;i++)
	{
		idx = edge_list[i].v;
		node_edge_count[idx]--;
		idx = csr_mat->row_ptr[idx] + node_edge_count[idx];
		csr_mat->col_ptr[idx] = edge_list[i].u;
		csr_mat->val_ptr[idx] = edge_list[i].w;
	}
	for(i=0;i<nodes_per_pe;i++)
		std::sort(&csr_mat->col_ptr[csr_mat->row_ptr[i]], &csr_mat->col_ptr[csr_mat->row_ptr[i+1]]);

	free(node_edge_count);
	free(edge_list);
	return csr_mat;
}