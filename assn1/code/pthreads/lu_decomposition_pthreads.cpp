#include<iostream>
#include<vector>
#include<chrono>
#include<time.h>
#include<algorithm>
#include<cmath>
#include<pthread.h>
#include<cmath>
#include <unistd.h>
#include <sys/syscall.h>
#define gettid() syscall(SYS_gettid)

using namespace std;
using namespace std::chrono;

// default number of threads = 4
int NUM_THREADS = 4;
int MATRIX_SIZE = 4;
// to implement static scheduling with chunk_size = 4
const int CHUNK_SIZE = 4;

// -------------- thread function declaration --------------//
void *init(void *data);

// -------------- Helper fuction declaration ---------------//
void print_matrix(double **mat, int size);
void print_vector(int *vec, int size);

// Through out the code we will work with this structure 
// so that we do not need to create datastructure again and again
// for parallel blocks
struct data_bin
{
	double **mat, **lower, **upper;
	int *permutation;
	data_bin()
	{
		// all matrices will point to corresponding array of MATRIX_SIZE size.
		// later each element of these array will point to corresponding array of double elements
		mat = 	(double **)(new double[MATRIX_SIZE]);
		lower = (double **)(new double[MATRIX_SIZE]);
		upper = (double **)(new double[MATRIX_SIZE]);
		permutation = new int[MATRIX_SIZE];
		// allocate memory to each array elements to create two dimensional matrices
		for(int i=0; i<MATRIX_SIZE; i++)
		{
			mat[i] = new double[MATRIX_SIZE];
			lower[i] = new double[MATRIX_SIZE];
			upper[i] = new double[MATRIX_SIZE];
			lower[i][i] = 1;
		}
	}
};

int main(int argc, char *argv[])
{
	// getting the current time
	auto start_timer = high_resolution_clock::now();

	// checking for correct command line input
	if(argc < 3)
	{
		cout << "Command line input expected in one of the following format.\n./a.out <size_of_matrix> <num_threads>\nmake A=<size_of_matrix> B=<num_threads>\n";
		return 0;
	}

	sscanf(argv[1], "%d", &MATRIX_SIZE);
	sscanf(argv[2], "%d", &NUM_THREADS);

	// data_bin has all matrices within it.
	data_bin data = data_bin();
	data_bin *data_pointer = &data;
	
	// implementing the following for intializing mat elements
	// #pragma omp parallel for schedule(static) num_threads(NUM_THREADS)
	pthread_t threads[NUM_THREADS];
	for(int i=0; i<NUM_THREADS; i++)	pthread_create(&threads[i], NULL, init, (void *)data_pointer);
	// adding barrier
	for(int i=0; i<NUM_THREADS; i++)	pthread_join(threads[i], NULL);

	// getting the end time
	auto end_time = high_resolution_clock::now();
	// get the duration
	auto duration_time = duration_cast<microseconds>(end_time - start_timer);
	cout << "Time taken: " << duration_time.count()/1000000.0 << " seconds" << endl;	

//	print_matrix(data.mat, MATRIX_SIZE);
	
	exit(0);
}

// implementing the following for intializing mat elements
// #pragma omp parallel for schedule(static) num_threads(NUM_THREADS)
// initialize mat matrix element
void *init(void *arg)
{
	// retrieving important data structures
	data_bin *data_pointer = (data_bin *)arg;
	double **mat = data_pointer->mat;
	// getting current thread's thread_id
	int tid = syscall(SYS_gettid) % NUM_THREADS;
	// figuring out the which iterations are allotted to this thread
	// ceiling is important. consider matrix_size=21, NUM_THREADS=4 then 1 iteration may be left to be executed without taking ceiling
	// it is used to incorporate ramainder iterations
	int start_iter_index = tid * ceil( ( (MATRIX_SIZE*(1.0)) / NUM_THREADS) );
	int end_iter_index = start_iter_index + (MATRIX_SIZE-1);
	end_iter_index = end_iter_index >= MATRIX_SIZE ? (MATRIX_SIZE-1) : end_iter_index;

	for(int i=start_iter_index; i<= end_iter_index; i++)
		for(int j=0; j<MATRIX_SIZE; j++)	mat[i][j] = i;
}

void print_matrix(double **mat, int size)
{
	for(int i=0; i<size; i++)
	{
		for(int j=0; j<size; j++)	cout << mat[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}
void print_vector(int *vec, int size)
{
	for(int i=0; i<size; i++)	cout << vec[i] << " ";
	cout << endl;
	cout << endl;
}
