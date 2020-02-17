#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

int main()
{
	// initializing the mpi
	MPI_Init(NULL, NULL);
	// getting the rank of the process
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	if(my_rank == 0)
	{
		time_t t;
		srand((unsigned) time(&t));
		int arr_size = (rand()/(float)RAND_MAX)*10;
		int arr[arr_size];
		for(int i=0; i<arr_size; i++)	arr[i] = i;
		printf("array created of size: %d\n", arr_size);
		MPI_Send(arr, arr_size, MPI_INT, 1, 0, MPI_COMM_WORLD);
	}
	else if(my_rank == 1)
	{
		MPI_Status status;
		// storing the size into the MPI_Status structure
		MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
		// getting the size into a variable
		int arr_size;
		MPI_Get_count(&status, MPI_INT, &arr_size);
		// receiving the message
		int arr[arr_size];
		MPI_Recv(arr, arr_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("got the message\n");
		for(int i=0; i<arr_size; i++)	printf("%d ",arr[i]);
		printf("\n");
	}
	MPI_Finalize();
	return 0;
}
