#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<time.h>
#include<unistd.h>

// give illusion of doing long work
int do_work(int sec)
{
	sleep(sec);
	return sec;
}

int main()
{
	// initializing mpi environment
	MPI_Init(NULL, NULL);
	// getting my rank
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	// some common variables
	int val;
	MPI_Request request;
	MPI_Status status;
	int request_accepted=0;
	// if my_rank == 0 then send, test and wait request
	if(my_rank == 0)
	{
		val = do_work(3);
		MPI_Isend(&val, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &request);
		printf("0:	request sent, moving on...\n");
		for(int i=0; i<5; i++)
		{
			do_work(1);
			printf("0:	testing after %d seconds.\n", 3 + i+1);
			printf("%d++++++\n", request_accepted);
			if(!request_accepted)	MPI_Test(&request, &request_accepted, &status);
			printf("%d------\n", request_accepted);
			if(request_accepted)	printf("0:	request completed.\n");
		}
		printf("0:	waiting for the request to complete...\n");
		if(!request_accepted)	MPI_Wait(&request, &status);
		if(request_accepted)	printf("0:	request accepted...\n");
	}
	else
	{
		time_t t;
		srand((unsigned) time(&t));
		int random = rand()%14;
		printf("1:	doing work for %d sec\n", random);
		// doing work for random amount of time
		do_work(random);
		MPI_Irecv(&val, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
		MPI_Wait(&request, &status);
		printf("1:	request accepted\n");
	}

	// finializing mpi environment
	MPI_Finalize();
}
