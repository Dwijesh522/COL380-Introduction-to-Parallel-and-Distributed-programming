#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

int main()
{
	// initializing the mpi environment
	MPI_Init(NULL, NULL);
	// getting the process rank
	int my_rank, arr_size=100;
	int arr[arr_size];
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	if(my_rank == 0)
		for(int i=0; i<arr_size; i++)	arr[i] = i+1;
	MPI_Bcast(arr, arr_size, MPI_INT, 0, MPI_COMM_WORLD);
	int sum=0;
	if(my_rank == 0)
		for(int i=0; i<arr_size/4; i++)	sum += arr[i];
	if(my_rank == 1)
		for(int i=arr_size/4; i<arr_size/2; i++)	sum += arr[i];
	if(my_rank == 2)
		for(int i=arr_size/2; i<3*arr_size/4; i++)	sum += arr[i];
	if(my_rank == 3)
		for(int i=3*arr_size/4; i<arr_size; i++)	sum += arr[i];
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank != 0)
		MPI_Send(&sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	else
	{
		int sum_p[3];
		for(int i=0; i<3; i++)
		{
			MPI_Recv(&sum_p[i], 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		sum += sum_p[0] + sum_p[1] + sum_p[2];
		printf("sum of 100 natural numbers: %d\n", sum);
	}
	MPI_Finalize();
}
