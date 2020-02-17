#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>

const int MAX_COUNT = 100;

int main()
{
	// initializing the mpi environemnt
	MPI_Init(NULL, NULL);
	int arr[MAX_COUNT], my_rank, message_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	if(my_rank == 0)
	{
		time_t t;
		srand((unsigned) time(&t));
		message_size = (rand()/(float)RAND_MAX) * MAX_COUNT;
		printf("sending message of size %d\n", message_size);
		MPI_Send(arr, message_size, MPI_INT, 1, 0, MPI_COMM_WORLD);
	}
	else if(my_rank == 1)
	{
		MPI_Status status;
		MPI_Recv(arr, MAX_COUNT, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_INT, &message_size);
		printf("received message of size: %d\n", message_size);
		printf("message was from source rank: %d\n", status.MPI_SOURCE);
		printf("message has the tag: %d", status.MPI_TAG);
	}
	return 0;
}
