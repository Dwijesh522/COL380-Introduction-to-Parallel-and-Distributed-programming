#include<stdio.h>
#include<mpi.h>
int main()
{
	// initializing mpi environment
	MPI_Init(NULL, NULL);
	// getting the rank of my process
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	// creating array of variable size
	int arr_size, arr_size_others;
	if(my_rank == 0)
	{
		scanf("%d", &arr_size);
		arr_size_others = arr_size/4;
	}
	MPI_Bcast(&arr_size_others, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	// master process
	if(my_rank == 0)
	{
		int arr[arr_size];
		int arr_o[arr_size_others];
		for(int i=0; i<arr_size; i++)	arr[i] = i+1;
		MPI_Scatter(arr, arr_size_others, MPI_INT, arr_o, arr_size_others, MPI_INT, 0, MPI_COMM_WORLD);

		int sum=0;
		for(int i=0; i<arr_size_others; i++)	sum += arr_o[i];

		int sum_p[4];
		MPI_Gather(&sum, 1, MPI_INT, sum_p, 1, MPI_INT, 0, MPI_COMM_WORLD);

		sum += sum_p[1] + sum_p[2] + sum_p[3];
		printf("sum of first %d elements: %d\n", arr_size, sum);
	}
	else
	{
		int arr_o[arr_size_others];
		MPI_Scatter(NULL, 0, MPI_INT, arr_o, arr_size_others, MPI_INT, 0, MPI_COMM_WORLD);
		
		int sum=0;
		for(int i=0; i<arr_size_others; i++)	sum += arr_o[i];

		MPI_Gather(&sum, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
	}

	// finalizing mpi environment
	MPI_Finalize();
}
