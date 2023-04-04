#include<iostream>
#include<mpi.h>
#include<unistd.h>
#include"skg.h"
void print_help(char);
void set_parameters(int, char**, int*, int*, block*, long*);
void print_parameters(int *scaling_factor, long *mat_size, block *mat_prob, long *nodes_per_pe, int *mat_blocks)
{
	printf("scaling_factor: %d\n", *scaling_factor);
	printf("matrix_size: %ld\n", *mat_size);
	printf("Probabity matrix (a: %0.3f, b: %0.3f, c: %0.3f, d: %0.3f)\n", mat_prob->a, mat_prob->b, mat_prob->c, mat_prob->d);
	printf("nodes_per_pe: %ld\n", *nodes_per_pe);
	printf("mat_blocks: %d\n", *mat_blocks);
}
int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int scaling_factor = 15, edge_factor = 20, rank, npes, mat_blocks;
	long mat_size, nodes_per_pe, pe_edges, edges=0, *edges_dist;
	time_stats graph_time = {0, 0, 0, 0};
	block mat_prob = {0.25, 0.25, 0.25, 0.25, 1};
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	set_parameters(argc, argv, &scaling_factor, &edge_factor, &mat_prob, &mat_size);
	nodes_per_pe = mat_size/npes;
	mat_blocks = mat_size/nodes_per_pe;
	printf("rank: %d\tnpes: %d\n", rank, npes);
	if(!rank)print_parameters(&scaling_factor, &mat_size, &mat_prob, &nodes_per_pe, &mat_blocks);
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
	mat_prob->nnz = (*mat_size)*(*edge_factor);
}