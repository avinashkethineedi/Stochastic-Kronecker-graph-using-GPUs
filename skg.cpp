#include<iostream>
#include<mpi.h>
#include<unistd.h>
#include"skg.h"
void print_help(char);
void set_parameters(int, char**, int*, int*, block*, long*);
void print_parameters(int *scaling_factor, long *mat_size, block *mat_prob, long *nodes_per_pe, int *mat_blocks, int edge_factor)
{
	printf("scaling_factor: %d\n", *scaling_factor);
	printf("matrix_size: %ld\n", *mat_size);
	printf("Probabity matrix (a: %0.3f, b: %0.3f, c: %0.3f, d: %0.3f)\n", mat_prob->a, mat_prob->b, mat_prob->c, mat_prob->d);
	printf("Total edges: %ld (edge_factor: %d)\n", mat_prob->edges, edge_factor);
	printf("nodes_per_pe: %ld\n", *nodes_per_pe);
	printf("mat_blocks: %d\n", *mat_blocks);
}
int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int scaling_factor = 15, edge_factor = 20, rank, npes, mat_blocks, *node_edge_count;
	long mat_size, nodes_per_pe, pe_edges, edges=0, *edges_dist, offset;
	edge* edge_list;
	csr_data *csr_mat;
	time_stats graph_time = {0, 0, 0, 0};
	block mat_prob = {0.25, 0.25, 0.25, 0.25, 1};
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	set_parameters(argc, argv, &scaling_factor, &edge_factor, &mat_prob, &mat_size);
	nodes_per_pe = mat_size/npes;
	mat_blocks = mat_size/nodes_per_pe;
	if(!rank)print_parameters(&scaling_factor, &mat_size, &mat_prob, &nodes_per_pe, &mat_blocks, edge_factor);
	//Edges distribution using given probabilty
	edges_dist = calculate_edge_distribution(rank, npes, &mat_prob);
	pe_edges = calculate_edges(edges_dist, npes);
	edge_list = create_edge_list(edges_dist, pe_edges, nodes_per_pe, npes, &mat_prob, &node_edge_count);
	csr_mat = create_csr_data(edge_list, node_edge_count, pe_edges, nodes_per_pe);
	MPI_Barrier(MPI_COMM_WORLD);
	//File write
	MPI_Exscan(&csr_mat->edges, &offset, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	if(!rank) offset=0;
	for(int i=0;i<=csr_mat->nodes;i++) csr_mat->row_ptr[i]+=offset; //update row_ptr
	offset = (nodes_per_pe+1)*rank*sizeof(long);
	file_write("row_ptr.bin", MPI_LONG, offset, csr_mat->row_ptr, nodes_per_pe+1, rank); //row_ptr file
	offset = csr_mat->row_ptr[0];
	file_write("col_ptr.bin", MPI_LONG, offset*sizeof(long), csr_mat->col_ptr, csr_mat->edges, rank); //col_ptr file
	file_write("val_ptr.bin", MPI_FLOAT, offset*sizeof(float), csr_mat->val_ptr, csr_mat->edges, rank); //val_ptr file
	printf("rank: %d\tnpes: %d, edges: %ld\n", rank, npes, pe_edges);
	MPI_Barrier(MPI_COMM_WORLD);
	free(csr_mat);
	MPI_Finalize();
	return EXIT_SUCCESS;
}
void print_help(char ch)
{
	if(ch == 's') printf("invalid value for -s (it requires positive integer)\n");
	if(ch == 'e') printf("invalid value for -e (it requires positive integer)\n");
	if(ch == 'a') printf("invalid value for -a (it requires positive value [0,1])\n");
	if(ch == 'b') printf("invalid value for -b (it requires positive value [0,1])\n");
	if(ch == 'c') printf("invalid value for -c (it requires positive value [0,1])\n");
	printf("scaling factor -s (15)\n");
	printf("edge factor -e (20)\n");
	printf("probability values: -a(0.57) -b(0.19) -c(0.19)\n");
	printf("help -h\n");
	MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}
void set_parameters(int argc, char **argv, int *scaling_factor, int *edge_factor, block *mat_prob, long *mat_size)
{
	int opt;
	while((opt = getopt(argc, argv, ":s:e:a:b:c:h:")) != -1)
	{
		switch(opt)
		{
			case 's':
				*scaling_factor=atoi(optarg);
				if(*scaling_factor < 1) print_help(opt);
				break;
			case 'e':
				*edge_factor=atoi(optarg);
				if(*edge_factor < 1) print_help(opt);
				break;
			case 'a':
				mat_prob->a = atof(optarg);
				if(mat_prob->a > 1 || mat_prob->a < 0) print_help(opt);
				break;
			case 'b':
				mat_prob->b = atof(optarg);
				if(mat_prob->b > 1 || mat_prob->b < 0) print_help(opt);
				break;
			case 'c':
				mat_prob->c = atof(optarg);
				if(mat_prob->c > 1 || mat_prob->c < 0) print_help(opt);
				break;
			case 'h':
				print_help(opt);
				break;
			default:
				if(optopt == 's') *scaling_factor = 15;
				else if(optopt == 'e') *edge_factor = 14;
				else if(optopt == 'a') mat_prob->a = 0.25;
				else if(optopt == 'b') mat_prob->b = 0.25;
				else if(optopt == 'c') mat_prob->c = 0.25;
				else print_help('h');
		}
	}
	mat_prob->d = 1 - (mat_prob->a + mat_prob->b + mat_prob->c);
	if(mat_prob->d > 1 || mat_prob->d < 0)
	{
		printf("sum of probabilities (a: %0.3f, b: %0.3f, c: %0.3f, d: %0.3f) != 1\n", mat_prob->a, mat_prob->b, mat_prob->c, mat_prob->d);
		print_help('h');
	}
	*mat_size = (long)1<<(*scaling_factor);
	mat_prob->edges = (*mat_size)*(*edge_factor);
}