#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<sys/sysinfo.h>
#include<math.h>
#include<stdbool.h>
#include<time.h>
// -------------- Assumption ----------------------
// if A : (n x m) then
// B : (m x n) and so C:(n x n)
// ------------------------------------------------

// --------------- Notation ----------------
// A: n x m, B: m x n, C: n x n
// p: number of processors available in pc
// -----------------------------------------

// ----------- helper fxns -----------------------------------------------------------------------------------------------------------
void init_matrix(double *m, int rows, int cols, bool rand_init);
void print_matrix(double *m, int rows, int cols);
void seq_matrix_multiplication_request_check(double *a, double *b, double *c, int m, int n, MPI_Request request[], int request_count);
void seq_matrix_multiplication(double *a, double *b, double *c, int m, int n);
bool isEqual(double *a, double *b, int r, int c);
// -----------------------------------------------------------------------------------------------------------------------------------

int main()
{
	int n, m, p, sqrt_p;
	double *A, *B, *C, *A_dash, *B_dash, *C_dash;
	double *ans_seq;
	clock_t time_par, time_seq;
	// initializing the MPI env
	MPI_Init(NULL, NULL);
	// getting my rank
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	// initializing parameters for creating the grid
	if(my_rank == 0)
	{
		// reading inputs
		printf("size of the matrix A:\n");
		scanf("%d %d", &n, &m);
		printf("number of processes:\n");
		scanf("%d", &p);
		sqrt_p = sqrt(p);
		// getting the # proc on this pc
		// initializing the matrix A and B
		double *big_A, *big_B, *big_C;
		big_A = (double *)malloc(n*m*sizeof(double));
		big_B = (double *)malloc(m*n*sizeof(double));
		A = (double *)malloc((n/sqrt_p)*(m/sqrt_p)*sizeof(double));
		B = (double *)malloc((n/sqrt_p)*(m/sqrt_p)*sizeof(double));
		C = (double *)malloc((n/sqrt_p)*(n/sqrt_p)*sizeof(double));
		A_dash = (double *)malloc((n/sqrt_p)*(m/sqrt_p)*sizeof(double));
		B_dash = (double *)malloc((n/sqrt_p)*(m/sqrt_p)*sizeof(double));
		C_dash = (double *)malloc((n/sqrt_p)*(n/sqrt_p)*sizeof(double));

		time_t t;
		srand((unsigned)time(&t));
		init_matrix(big_A, n, m, true);
		init_matrix(big_B, m, n, true);
		init_matrix(C, n/sqrt_p, n/sqrt_p, false);
		// calculating the sequential answer
		ans_seq = (double *)malloc(n*n*sizeof(double));
		printf("Starting seq multiplication...\n");
		time_seq = clock();
		seq_matrix_multiplication(big_A, big_B, ans_seq, n, m);
		time_seq = clock() - time_seq;
		printf("Starting par multiplication...\n");
		
		// temp matrix of matrix for sending info to others
		double *temp_A[p-1], *temp_B[p-1];
		MPI_Request request[3*(p-1)];
		MPI_Status status[3*(p-1)];
		
		int temp_par[3];
		temp_par[0] = n;
		temp_par[1] = m;
		temp_par[2] = p;
		int request_counter = 0;
		int index_count;
		int indices[p-1];
		// let other process know n, m, p
		time_par = clock();
		for(int k=1; k<p; k++)
		{
			temp_A[k] = (double *)malloc((n/sqrt_p)*(m/sqrt_p)*sizeof(double));
			temp_B[k] = (double *)malloc((n/sqrt_p)*(m/sqrt_p)*sizeof(double));

			MPI_Isend(temp_par, 3, MPI_INT, k, 0, MPI_COMM_WORLD, &request[request_counter++]);
			// cart coordinates of the porcess i
			int proc_cart_coord[2];
			proc_cart_coord[0] = k/sqrt_p;
			proc_cart_coord[1] = k%sqrt_p;
			// corresponding matrix coordinates for process i
			int A_coord[2], B_coord[2];
			A_coord[0] = proc_cart_coord[0]*n/sqrt_p;
			A_coord[1] = proc_cart_coord[1]*m/sqrt_p;
			B_coord[0] = proc_cart_coord[0]*m/sqrt_p;
			B_coord[1] = proc_cart_coord[1]*n/sqrt_p;
			// coping corresponding elements for sending to process i
			for(int i=0; i<n/sqrt_p; i++)
				for(int j=0; j<m/sqrt_p; j++)
					temp_A[k][i*(m/sqrt_p) + j] = big_A[(A_coord[0]+i)*m + A_coord[1]+j];
			for(int i=0; i<m/sqrt_p; i++)
				for(int j=0; j<n/sqrt_p; j++)
					temp_B[k][i*(n/sqrt_p) + j] = big_B[(B_coord[0]+i)*n + B_coord[1]+j];
			MPI_Isend(temp_A[k], m*n/p, MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &request[request_counter++]);
			MPI_Isend(temp_B[k], m*n/p, MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &request[request_counter++]);
			MPI_Testsome(3*(k), request, &index_count, indices, MPI_STATUS_IGNORE);
		}
		// coping my own A, B.
		for(int i=0; i<n/sqrt_p; i++)
			for(int j=0; j<m/sqrt_p; j++)
				A[i*(m/sqrt_p) + j] = big_A[i*m + j];
		for(int i=0; i<m/sqrt_p; i++)
			for(int j=0; j<n/sqrt_p; j++)
				B[i*(n/sqrt_p) + j] = big_B[i*n + j];
		printf("waiting for request...\n");
		// waiting for the requests to be accepted
		MPI_Waitall(3*(p-1), request, status);
		// freeing extra memory
		free(big_A);
		free(big_B);
	}
	else
	{
		int temp_par[3];
		MPI_Recv(temp_par, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		n = temp_par[0];
		m = temp_par[1];
		p = temp_par[2];
		// size of the blocked matrices
		sqrt_p = sqrt(p);
		// getting the block matrices
		A = (double *)malloc((n/sqrt_p)*(m/sqrt_p)*sizeof(double));
		B = (double *)malloc((n/sqrt_p)*(m/sqrt_p)*sizeof(double));
		C = (double *)malloc((n/sqrt_p)*(n/sqrt_p)*sizeof(double));
		A_dash = (double *)malloc((n/sqrt_p)*(n/sqrt_p)*sizeof(double));
		B_dash = (double *)malloc((n/sqrt_p)*(n/sqrt_p)*sizeof(double));
		MPI_Recv(A, n*m/p, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(B, n*m/p, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("receiver: message received\n");
		init_matrix(C, n/sqrt_p, n/sqrt_p, false);
	}
	// +----------------------------------------------------------------------------------------------------+
	// |			At this point processes have their own local block of matrix 			|
	// +----------------------------------------------------------------------------------------------------+
	
	
	
	// ----------------- creating 2d matrix of mpi processes ---------------
	MPI_Comm MPI_COMM_WORLD_2D;
	int dims[2];
	dims[0] = sqrt_p;
	dims[1] = sqrt_p;
	int periods[2];
	periods[0] = 1;
	periods[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &MPI_COMM_WORLD_2D);
	// cartesian coord
	int cart_coord[2];
	MPI_Cart_coords(MPI_COMM_WORLD_2D, my_rank, 2, cart_coord);
	bool is_even = 1 - (cart_coord[0] + cart_coord[1])%2;
	//-----------------------------------------------------------------------
	
	
	
	// -------------------------------------- initial shift -----------------------------------------------------------------------
	int rank_source, rank_dest;
	bool A_dash_active = false;
	bool B_dash_active = false;
	MPI_Request shift_request[2];
	MPI_Status shift_status[2];
	int request_count = 0;
	// !!!!!   direction=1 for row shift,	direction=0 for column shift !!!!!!!!!!
	// row shifting
	MPI_Cart_shift(MPI_COMM_WORLD_2D, 1, -cart_coord[0], &rank_source, &rank_dest);
	if(is_even && my_rank != rank_dest)
	{
		MPI_Send(A, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_dest, 0, MPI_COMM_WORLD);
		MPI_Recv(A_dash, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		A_dash_active = true;
	}
	else if(my_rank != rank_dest)
	{	
		MPI_Recv(A_dash, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Isend(A, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_dest, 0, MPI_COMM_WORLD, &shift_request[request_count++]);
		A_dash_active = true;
	}
	// col shifting
	MPI_Cart_shift(MPI_COMM_WORLD_2D, 0, -cart_coord[1], &rank_source, &rank_dest);
	if(is_even && my_rank != rank_dest)
	{
		MPI_Send(B, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_dest, 0, MPI_COMM_WORLD);
		MPI_Recv(B_dash, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		B_dash_active = true;
	}
	else if(my_rank != rank_dest)
	{
		MPI_Recv(B_dash, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Isend(B, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_dest, 0, MPI_COMM_WORLD, &shift_request[request_count++]);
		B_dash_active = true;
	}
	//-----------------------------------------------------------------------------------------------------------------------------
	
	
	
	// --------------------------------------------------------------------- block multiplication ----------------------------------------------------------------------
	for(int i=0; i<sqrt_p; i++)
	{
		if(A_dash_active && B_dash_active)		
			seq_matrix_multiplication_request_check(A_dash, B_dash, C, (n/sqrt_p), (m/sqrt_p), shift_request, request_count);
		else if(A_dash_active && !B_dash_active)	
			seq_matrix_multiplication_request_check(A_dash, B, C, (n/sqrt_p), (m/sqrt_p), shift_request, request_count);
		else if(!A_dash_active && B_dash_active)	seq_matrix_multiplication_request_check(A, B_dash, C, (n/sqrt_p), (m/sqrt_p), shift_request, request_count);
		else						seq_matrix_multiplication_request_check(A, B, C, (n/sqrt_p), (m/sqrt_p), shift_request, request_count);
		// send local copy of matrx A to left by 1 unit
		if(!is_even)	MPI_Waitall(request_count, shift_request, shift_status);
		request_count = 0;
		MPI_Cart_shift(MPI_COMM_WORLD_2D, 1, -1, &rank_source, &rank_dest);
		if(is_even)
		{
			if(A_dash_active)
			{
				MPI_Send(A_dash, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_dest, 0, MPI_COMM_WORLD);
				MPI_Recv(A, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				A_dash_active = false;
			}
			else
			{
				MPI_Send(A, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_dest, 0, MPI_COMM_WORLD);
				MPI_Recv(A_dash, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				A_dash_active = true;
			}
		}
		else
		{
			if(A_dash_active)
			{
				MPI_Recv(A, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Isend(A_dash, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_dest, 0, MPI_COMM_WORLD, &shift_request[request_count++]);
				A_dash_active = false;
			}
			else
			{
				MPI_Recv(A_dash, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Isend(A, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_dest, 0, MPI_COMM_WORLD, &shift_request[request_count++]);
				A_dash_active = true;
			}
		}
		// send local copy of matrix B to up by 1 unit
		MPI_Cart_shift(MPI_COMM_WORLD_2D, 0, -1, &rank_source, &rank_dest);
		if(is_even)
		{
			if(B_dash_active)
			{
				MPI_Send(B_dash, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_dest, 0, MPI_COMM_WORLD);
				MPI_Recv(B, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				B_dash_active = false;
			}
			else
			{
				MPI_Send(B, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_dest, 0, MPI_COMM_WORLD);
				MPI_Recv(B_dash, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				B_dash_active = true;
			}
		}
		else
		{
			if(B_dash_active)
			{
				MPI_Recv(B, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Isend(B_dash, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_dest, 0, MPI_COMM_WORLD, &shift_request[request_count++]);
				B_dash_active = false;
			}
			else
			{
				MPI_Recv(B_dash, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Isend(B, (n/sqrt_p)*(m/sqrt_p), MPI_DOUBLE, rank_dest, 0, MPI_COMM_WORLD, &shift_request[request_count++]);
				B_dash_active = true;
			}
		}
	}
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// +----------------------------------------------------------------------------------------------------+
	// |			At this point processes have calculated their own blocked ans 			|
	// +----------------------------------------------------------------------------------------------------+
	
	

	// -------------------------------- getting the complete ans in process 0 -------------------------------------
	double *ans;
	if(my_rank == 0)
	{
		ans = (double *)malloc(n*n*sizeof(double));
		for(int i=1; i<p; i++)
		{
			// reading block ans
			double *block_ans;
			block_ans = (double *)malloc((n/sqrt_p)*(n/sqrt_p)*sizeof(double));
			MPI_Recv(block_ans, (n/sqrt_p)*(n/sqrt_p), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// cart coordinates of the porcess i
			int proc_cart_coord[2];
			proc_cart_coord[0] = i/sqrt_p;
			proc_cart_coord[1] = i%sqrt_p;
			// corresponding matrix coordinates for process i
			int C_coord[2];
			C_coord[0] = proc_cart_coord[0]*n/sqrt_p;
			C_coord[1] = proc_cart_coord[1]*n/sqrt_p;
			// coping corresponding elements for sending to process i
			for(int i=0; i<n/sqrt_p; i++)
				for(int j=0; j<n/sqrt_p; j++)
					ans[ (C_coord[0]+i)*n + (C_coord[1]+j) ] = block_ans[i*(n/sqrt_p) + j];
		}
		for(int i=0; i<n/sqrt_p; i++)
			for(int j=0; j<n/sqrt_p; j++)
				ans[ i*n + j ] = C[i*(n/sqrt_p) + j];
		time_par = clock() - time_par;
		free(C);
	}
	else
	{
		MPI_Send(C, (n/sqrt_p)*(n/sqrt_p), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		free(A);
		free(B);
		free(C);
		free(A_dash);
		free(B_dash);
	}
	// ---------------------------------------------------------------------------------------------------------------
	// +----------------------------------------------------------------------------------------------------+
	// |				At this point process 0 has the complete ans in 'ans' 			|
	// +----------------------------------------------------------------------------------------------------+



	// -------------------------------------- check with sequential output--------------------------------------------
	if(my_rank == 0)
	{
		printf("checking the correctness...\n");
		if(isEqual(ans, ans_seq, n, n))		printf("Sequential and parallel multiplication gives same answer.\n");
		else					printf("Sequential and parallel multiplication gives different answer.\n");

		double time_seq_sec = ((double)time_seq)/CLOCKS_PER_SEC;
		double time_par_sec = ((double)time_par)/CLOCKS_PER_SEC;
		printf("Sequntial time: %lf\n", time_seq_sec);
		printf("Parallel time: %lf\n", time_par_sec);
	}
	//----------------------------------------------------------------------------------------------------------------
	


	// finalizing MPI env
	MPI_Finalize();
	return 0;
}

// comaper if both the matrices are same or not
bool isEqual(double *a, double *b, int r, int c)
{
	for(int i=0; i<r; i++)
		for(int j=0; j<c; j++)
			if(abs(a[i*c + j] - b[i*c + j]) >= 0.000001)
			{
				printf("wrong at: (%d, %d).\n", i, j);
				printf("values are: %lf\n", a[i*c+j]-b[i*c+j]);
				return false;
			}
	return true;
}

// does c = a*b, a: mxn, b: nxm, c: mxm
void seq_matrix_multiplication_request_check(double *a, double *b, double *c, int m, int n, MPI_Request request[], int request_count)
{
	int accepted1 = 0, accepted2 = 0;
	for(int i=0; i<m; i++)
		for(int j=0; j<m; j++)
		{
			for(int k=0; k<n; k++)				c[i*m + j] += a[i*n + k]*b[k*m + j];
			if(!accepted1 && request_count >= 1)		MPI_Test(&request[0], &accepted1, MPI_STATUS_IGNORE);
			if(!accepted2 && request_count >= 2)		MPI_Test(&request[1], &accepted2, MPI_STATUS_IGNORE);
		}
}

// does c = a*b, a: mxn, b: nxm, c: mxm
void seq_matrix_multiplication(double *a, double *b, double *c, int m, int n)
{
	for(int i=0; i<m; i++)
		for(int j=0; j<m; j++)
			for(int k=0; k<n; k++)
				c[i*m + j] += a[i*n + k]*b[k*m + j];
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
