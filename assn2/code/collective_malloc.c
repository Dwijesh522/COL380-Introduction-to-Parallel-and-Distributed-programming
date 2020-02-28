#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<mpi.h>
#include<stdbool.h>

//A and C are row of columns. B is column of rows
void Matrix_Multiply_Seq(float *A, float *B, float *C, int m, int n, int p);
//Both the matrices are row of columns
int IsEqual(float *A, float *B, int m, int n);
float* create_matrix_space(int x, int y);
void create_matrix(bool random, float* mat, int X, int Y);
void print_matrix(float* mat, int X, int Y);
int main(int argc, char** argv)
{
  //Initialize MPI environment
  MPI_Init(NULL, NULL);
  // -------------- Assumption ----------------------
  // if A : (n x m) then
  // B : (m x n) and so C:(n x n)
  // ------------------------------------------------
  int n,m,p;
  //GET rank and size of World
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  float* A, *C, *C_serial;
  float* B;
  if(world_rank == 0)
  {
    //Take values of m and n from standard input. P is number of processes
    p = world_size;
    printf("Give the values of n and m\n");
    scanf("%d %d",&n, &m);
  }
  //ALL COLLECTIVE COMMUNICATION ROUTINES IMPLY SYNCHRONIZATION POINTS
  //broadcasting n,m,p
  MPI_Bcast((void*)&n, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast((void*)&m, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast((void*)&p, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  // // Checking consistency of n,m,p
  // if (n%p != 0)
  // {
  //   if(world_rank == 0)
  //     printf("Give valid input size N divisible by %d in process %d\n", p, world_rank);
  //   MPI_Finalize();
  //   return 0;
  // }
  //Serial Multiplication
  if(world_rank == 0)
  {
    //time_seq is sequential multiplication time time
    double seq_time = 0;
    clock_t time_seq_start, time_seq_end;
    //Initialize matrices A,B and C, C_serial
    A = create_matrix_space(n,m);
    B= create_matrix_space(m, n);
    C = create_matrix_space(n,n);
    C_serial = create_matrix_space(n,n);
    //Generate random matrices A and B between 0 and 1
    srand(time(NULL));
    create_matrix(true, A, n, m);
    create_matrix(true, B, m, n);
    create_matrix(false, C_serial, n, n);
    create_matrix(false, C, n, n);
    //Time the serial multiplication function implementation
    time_seq_start = clock();
    Matrix_Multiply_Seq(A, B, C_serial, n, n, m);
    time_seq_end = clock();
    // print_matrix(C_serial, n, n);
    seq_time = ((double)(time_seq_end - time_seq_start))/CLOCKS_PER_SEC;
    printf("The time for sequential multiplication with N=%d is %lf\n", n, seq_time);
  }
  else
  {
    B = create_matrix_space(m,n);
  }


  int sub_matrix_size;
  float* sub_matrix_A, *sub_matrix_C;
  //Parallel computation Begins
  if(world_rank==0)
  {
     //time_parallel is parallel implementation
    double parallel_time = 0;
    clock_t time_parallel_start, time_parallel_end;
    time_parallel_start = clock();
    //Broadcasting matrix B to all processors
    MPI_Bcast((void*)B, m*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    //sub_matrix_A
    sub_matrix_size = n/p;
    sub_matrix_A = create_matrix_space(sub_matrix_size, m);
    //Scattering A to sub_matrix_A
      MPI_Scatter((void*)A, sub_matrix_size*m, MPI_FLOAT, (void*)sub_matrix_A, sub_matrix_size*m,
      MPI_FLOAT, 0, MPI_COMM_WORLD);
    //submatrix of C
    sub_matrix_C = create_matrix_space(sub_matrix_size, n);
    create_matrix(false, sub_matrix_C, sub_matrix_size, n);
    //Obtaining sub matrix of final answer C
    Matrix_Multiply_Seq(sub_matrix_A, B, sub_matrix_C, sub_matrix_size, n, m);
    //Gather answer to rank 0
      MPI_Gather((void*)sub_matrix_C, sub_matrix_size*n, MPI_FLOAT, (void*)C, sub_matrix_size*n,
      MPI_FLOAT, 0, MPI_COMM_WORLD);
    //Computation for remainder number of rows
    int rem = n%p;
    if (rem != 0)
    {
      int index_last_row = n - rem;
      MPI_Group new_gp, world_gp;
      MPI_Comm_group(MPI_COMM_WORLD, &world_gp);
      MPI_Comm new_comm;
      int new_comm_ranks[rem];
      for(int i=0; i<rem; i++) new_comm_ranks[i] = i;
      MPI_Group_incl(world_gp, rem, new_comm_ranks, &new_gp);
      MPI_Comm_create(MPI_COMM_WORLD, new_gp, &new_comm);
      if(new_comm != MPI_COMM_NULL)
      {
        //Scattering each of the row to each processor
        MPI_Scatter((void*)(A+index_last_row*m), m, MPI_FLOAT, (void*)sub_matrix_A, m, MPI_FLOAT, 0, new_comm);
        //multiplying the matrices
        Matrix_Multiply_Seq(sub_matrix_A, B, sub_matrix_C, 1, n, m);
        //Gathering back the answer
        MPI_Gather((void*)sub_matrix_C, n, MPI_FLOAT, (void*)(C+index_last_row*n), n, MPI_FLOAT, 0, new_comm);
      }
    }
    time_parallel_end = clock();
    parallel_time = ((double)(time_parallel_end - time_parallel_start))/CLOCKS_PER_SEC;
    printf("The time for parallel multiplication with N=%d is %lf\n", n, parallel_time);
    //Check if both the impleentations are correct
    if(IsEqual(C, C_serial, n, n))
    {
      printf("The two matrices are equal\n");
    }
  }
  else
  {

    //Broadcasting matrix B to all processors
      MPI_Bcast((void*)B, m*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    //sub_matrix_A
    sub_matrix_size = n/p;
    sub_matrix_A = create_matrix_space(sub_matrix_size, m);
    //Scattering A to sub_matrix_A
      MPI_Scatter((void*)A, sub_matrix_size*m, MPI_FLOAT, (void*)sub_matrix_A, sub_matrix_size*m,
      MPI_FLOAT, 0, MPI_COMM_WORLD);
    //submatrix of C
    sub_matrix_C = create_matrix_space(sub_matrix_size, n);
    create_matrix(false, sub_matrix_C, sub_matrix_size, n);
    //Obtaining sub matrix of final answer C
    Matrix_Multiply_Seq(sub_matrix_A, B, sub_matrix_C, sub_matrix_size, n, m);
    //Gather answer to rank 0
      MPI_Gather((void*)sub_matrix_C, sub_matrix_size*n, MPI_FLOAT, (void*)C, sub_matrix_size*n,
      MPI_FLOAT, 0, MPI_COMM_WORLD);

    //Computation for remainder number of rows
    int rem = n%p;
    if (rem != 0)
    {
      int index_last_row = n - rem;
      MPI_Group new_gp, world_gp;
      MPI_Comm_group(MPI_COMM_WORLD, &world_gp);
      MPI_Comm new_comm;
      int new_comm_ranks[rem];
      for(int i=0; i<rem; i++) new_comm_ranks[i] = i;
      MPI_Group_incl(world_gp, rem, new_comm_ranks, &new_gp);
      MPI_Comm_create(MPI_COMM_WORLD, new_gp, &new_comm);
      if(new_comm != MPI_COMM_NULL)
      {
        //Scattering each of the row to each processor
        MPI_Scatter((void*)(A+index_last_row*m), m, MPI_FLOAT, (void*)sub_matrix_A, m, MPI_FLOAT, 0, new_comm);
        //multiplying the matrices
        Matrix_Multiply_Seq(sub_matrix_A, B, sub_matrix_C, 1, n, m);
        //Gathering back the answer
        MPI_Gather((void*)sub_matrix_C, n, MPI_FLOAT, (void*)(C+index_last_row*n), n, MPI_FLOAT, 0, new_comm);
      }
    }
  }

  //End MPI
  MPI_Finalize();
  return 0;
}

//A and C are row of columns. B is column of rows
void Matrix_Multiply_Seq(float *A, float *B, float *C, int m, int n, int p){
	int i, j, k;
	for (i = 0; i < m; i++){
		for (j = 0; j < n; j++){
      float temp = 0;
			for (k = 0; k < p; k++)
				temp += A[k+i*p]*B[j+k*n];
      C[j+i*n] = temp;
		}
	}
  return;
}

//Both the matrices are row of columns
int IsEqual(float *A, float *B, int X, int Y)
{
  for(int i=0; i<X; i++)
  {
    for(int j=0; j<Y; j++)
      {
        if(abs(A[j+i*Y] - B[j+Y*i]) > 0.0000000001) return 0;
      }
  }
  return 1;
}

float* create_matrix_space(int X, int Y)
{
  float* mat = (float*)malloc(sizeof(float)*X*Y);
  return mat;
}

void create_matrix(bool random, float* mat, int X, int Y)
{
  if(random)
  {
    for(int i=0; i<X; i++)
    {
      for(int j=0; j<Y; j++)
      {
        mat[j+i*Y] = rand()/(double)RAND_MAX;
      }
    }
  }
  else
  {
    for(int i=0; i<X; i++)
    {
      for(int j=0; j<Y; j++)
      {
        mat[j+i*Y]  = 0;
      }
    }
  }
}

void print_matrix(float* mat, int X, int Y)
{
  for (int i=0; i<X; i++)
  {
    for (int j=0; j<Y; j++)
      printf("%f ", mat[j+i*Y]);
    printf("\n");
  }
  return;
}
