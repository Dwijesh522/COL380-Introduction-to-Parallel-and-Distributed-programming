#include<stdio.h>
#include<mpi.h>

// initialiing array by paralleliing it.
int main()
{
	MPI_Init(NULL, NULL);
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	long long int arr_size, arr_chunk_size;
	// read input from stdin
	if(my_rank == 0)
	{
		scanf("%lld", &arr_size);
		arr_chunk_size = arr_size/4;
	}
	// bcasting arr_chunk_size
	MPI_Bcast(&arr_chunk_size, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	int arr[arr_chunk_size];
	for(int i=0; i<arr_chunk_size; i++)	arr[i]	= my_rank*arr_chunk_size + i+1;

	int complete_arr[4*arr_chunk_size];
	MPI_Allgather(arr, arr_chunk_size, MPI_INT, complete_arr, arr_chunk_size, MPI_INT, MPI_COMM_WORLD);
	for(int i=0; i<arr_chunk_size*4; i++)	printf("%d ", complete_arr[i]);
	printf("\n");
	MPI_Finalize();
}
