#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<sys/sysinfo.h>
#include<math.h>
#include<stdbool.h>
// --------------- Notation ----------------
// A: n x m, B: m x n, C: n x n
// p: number of processors available in pc
// -----------------------------------------

// ----------- helper fxns -------------------------------------
void init_matrix(double *m, int rows, int cols, bool rand_init);
void print_matrix(double *m, int rows, int cols);
// -------------------------------------------------------------

int main()
{
	// initializing the MPI env
	MPI_Init(NULL, NULL);
	// getting my rank
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	if(my_rank == 0)
	{
		// reading inputs
		int n, m, p;
		printf("size of the matrix A:\n");
		scanf("%d %d", &n, &m);
		// getting the # proc on this pc
		p = get_nprocs();
		int sqrt_p = sqrt(p);
		// initializing the matrix A and B
		double *A, *B, *C;
		A = (double *)malloc(n*m*sizeof(double));
		B = (double *)malloc(m*n*sizeof(double));
		C = (double *)malloc(n*n*sizeof(double));

		time_t t;
		srand((unsigned)time(&t));
		init_matrix(A, n, m, true);
		init_matrix(B, m, n, true);
		init_matrix(C, n, n, false);
	}

	// finalizing MPI env
	MPI_Finalize();
	return 0;
}

// initialized a matrix m with random values if rand_init is 1
// else with 0
void init_matrix(double *m, int rows, int cols, bool rand_init)
{
	if(!rand_init)	for(int i=0; i<rows; i++)	for(int j=0; j<cols; j++)	m[i*cols + j] = 0;
	else	for(int i=0; i<rows; i++)	for(int j=0; j<cols; j++)	m[i*cols + j] = (rand()/(double)RAND_MAX);
}

void print_matrix(double *m, int rows, int cols)
{
	for(int i=0; i<rows; i++)
	{
		for(int j=0; j<cols; j++)
			printf("%lf ", m[i*cols + j]);
		printf("\n");
	}
	printf("\n");
}
