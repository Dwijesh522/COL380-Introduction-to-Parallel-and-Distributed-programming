#include<iostream>
#include<omp.h>
#include<vector>
#include<chrono>
#include<time.h>

using namespace std;
using namespace std::chrono;

// default number of threads = 4
int NUM_THREADS = 4;
int MATRIX_SIZE = 4;

// Helper fuctions
void print_matrix(vector<vector<float*>> mat, int r, int c);
void print_vector(vector<int*> permutation, int size);
int get_maximum_element_index(vector<vector<float*>> &mat, int c, int r1, int r2);

int main(int argc, char *argv[])
{
	// getting the current time
	auto start = high_resolution_clock::now();
	
	// checking for correct command line input
	if(argc < 3)	
	{
		cout << "Command line input expected in one of the following format.\n./a.out <size_of_matrix> <num_threads>\nmake A=<size_of_matrix> B=<num_threads>\n";
		return 0;
	}
	
	sscanf(argv[1], "%d", &MATRIX_SIZE);
	sscanf(argv[2], "%d", &NUM_THREADS);
	
	// below format will help us parallelize initialization of matrices.
	// declaring matrices as two dimensional vector of pointer to floats
	vector<vector<float *>> mat(MATRIX_SIZE, vector<float*>(MATRIX_SIZE));
	vector<vector<float *>> upper(MATRIX_SIZE, vector<float*>(MATRIX_SIZE, new float(0)));
	vector<vector<float *>> lower(MATRIX_SIZE, vector<float*>(MATRIX_SIZE, new float(0)));
	// declaring permutation matrix as vector of pointers to integers
	vector<int *> permutation(MATRIX_SIZE);

	// seed the random number generator
	srand48(time(0));

	// -----------------------------------------------------------------------
	// ------------- sequential timing	: 3.5 sec
	// ------------- parallel timing	: 2.9 sec(optimal with 4 threads)
	// -----------------------------------------------------------------------
	#pragma omp parallel num_threads(NUM_THREADS)
	{
		// Initializing mat: 64000000 elements with random floating point values
		#pragma omp for schedule(static, 4) collapse(2) nowait
		for(int i=0; i<MATRIX_SIZE; i++)
			for(int j=0; j<MATRIX_SIZE; j++)
				mat[i][j] = new float(drand48());	// pointer to a random float value
		// Initializing permutation matrix and diagonal entries of lower triangular matrix: 8000 elements
		// No need to parallelize it
		#pragma omp for schedule(static, 4)
		for(int i=0; i<MATRIX_SIZE; i++)
		{
			permutation[i] = new int(i+1);
			lower[i][i] = new float(1);
		}
	}
	
	// outermost loop: non parallelizable due to data dependecies across iterations
	for(int k=0; k< MATRIX_SIZE; k++)
	{
		int k_dash = get_maximum_element_index(mat, k, k, MATRIX_SIZE-1);
		// if singular matrix then stop
		if(k_dash == -1)	{ cout << "Singular matrix.\n"; return 0;}
	}
	
	// getting the end time
	auto end = high_resolution_clock::now();
	// get the duration
	auto duration = duration_cast<microseconds>(end - start);
	cout << "Time taken: " << duration.count() << " micro seconds" << endl;

	return 0;
}

// This function finds the maximum element from a fix column and between raw index r1 nad r2
// If such non-zero element exists then it returns its index
// Else it returns -1.
// While comparing absolute values is compared against each other
// -----------------------------------------------------------------------
// ------------- sequential time:	4.9 seconds
// ------------- parallel time:		4.9 seconds
// -----------------------------------------------------------------------
int get_maximum_element_index(vector<vector<float*>> &mat, int c, int r1, int r2)
{
	float max_element=0;
	int index=-1;
	// parallel if enough number of elements
	if(r2-r1 >= 7700)
	{
		#pragma omp parallel for schedule(static, 4) num_threads(NUM_THREADS) reduction(max: index)
		for(int i= r1; i<= r2; i++)
		{
			index = max_element > abs(*mat[i][c]) ? index : i;
			max_element = max_element > (*mat[i][c]) ? max_element : (*mat[i][c]);
		}
	}
	// sequential for smaller number elements
	else
	{
		for(int i= r1; i<= r2; i++)
		{
			if(abs(*mat[i][c]) > max_element)
			{
				max_element = abs(*mat[i][c]);
				index = i;
			}
		}	
	}
	return index;
}

void print_matrix(vector<vector<float*>> mat, int r, int c)
{
	for(int i=0; i<r; i++)
	{
		for(int j=0; j<c; j++)	cout << *mat[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}
void print_vector(vector<int*> permutation, int size)
{
	for(int i=0; i<size; i++)	cout << *permutation[i] << " ";
	cout << endl << endl;
}
