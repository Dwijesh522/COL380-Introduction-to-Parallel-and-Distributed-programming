#include<mpi.h>
#include<stdio.h>

int main()
{
	// initializing mpi
	MPI_Init(NULL, NULL);
	int number;
	// getting the current process's id
	int process_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
	if(process_rank == 0)
	{
		number = -1;
		MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	}
	else if(process_rank == 1)
	{
		MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("message received: %d\n", number);
	}
	MPI_Finalize();
}
