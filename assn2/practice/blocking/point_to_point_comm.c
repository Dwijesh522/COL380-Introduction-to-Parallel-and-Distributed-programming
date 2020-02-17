#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

const int STEPS_LIMIT = 100;

int main()
{
	// initializing the mpi environment
	MPI_Init(NULL, NULL);
	// getting the rank of my process
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	// generating a random number for random walk size
	time_t t;
	srand((unsigned) time(&t));
	int walk_steps = (rand()/(float)RAND_MAX)*STEPS_LIMIT;
	int steps_division = STEPS_LIMIT/4;
	// temp = steps_limit/4
	// process 1: 0 	- temp-1
	// process 2: temp	- 2*temp-1 
	// process 3: 2*temp 	- 3*temp-1
	// process 4: 3*temp 	- 4*temp-1
	int steps = 0;
	// rank of the currently active process for random walk
	int curr_rank=0;
	if(my_rank == 0)	printf("steps to be taken: %d\n", walk_steps);
	while(steps < walk_steps)
	{
		if(my_rank == curr_rank)
		{
			steps++;
			if(steps == (my_rank+1)*steps_division-1)
			{
				// time to change the duty
				curr_rank++;
				printf("process %d completed %d steps.\n", curr_rank, steps+1);
				
				MPI_Send(&steps, 1, MPI_INT, (my_rank+1)%4, 0, MPI_COMM_WORLD);
				MPI_Send(&steps, 1, MPI_INT, (my_rank+2)%4, 0, MPI_COMM_WORLD);
				MPI_Send(&steps, 1, MPI_INT, (my_rank+3)%4, 0, MPI_COMM_WORLD);
				
				MPI_Send(&curr_rank, 1, MPI_INT, (my_rank+1)%4, 0, MPI_COMM_WORLD);
				MPI_Send(&curr_rank, 1, MPI_INT, (my_rank+2)%4, 0, MPI_COMM_WORLD);
				MPI_Send(&curr_rank, 1, MPI_INT, (my_rank+3)%4, 0, MPI_COMM_WORLD);
			}
		}
		else
		{
			MPI_Recv(&steps, 1, MPI_INT, curr_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if(steps == walk_steps)	break;
			MPI_Recv(&curr_rank, 1, MPI_INT, curr_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	if(my_rank == curr_rank)
	{
		printf("process %d completed %d stpes.\n", curr_rank, steps);
		MPI_Send(&steps, 1, MPI_INT, (my_rank+1)%4, 0, MPI_COMM_WORLD);
		MPI_Send(&steps, 1, MPI_INT, (my_rank+2)%4, 0, MPI_COMM_WORLD);
		MPI_Send(&steps, 1, MPI_INT, (my_rank+3)%4, 0, MPI_COMM_WORLD);
	}	
	MPI_Finalize();
}
