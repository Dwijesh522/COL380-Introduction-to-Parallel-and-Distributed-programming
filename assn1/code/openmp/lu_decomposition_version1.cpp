#include<iostream>
#include<omp.h>
#include<vector>
#include<chrono>
#include<time.h>
#include<algorithm>

using namespace std;
using namespace std::chrono;

// default number of threads = 4
int NUM_THREADS = 4;
int MATRIX_SIZE = 4;

// Helper fuctions
void print_matrix(vector<vector<double>* > mat, int r, int c);
void print_vector(vector<int> permutation, int size);
pair<int, double> get_maximum_element_index(vector<vector<double>* >& mat, int c, int r1, int r2);
void print_error(vector<vector<double>*> &mat_dup, vector<int> &permutation, vector<vector<double>*> &lower, vector<vector<double>*> &upper);

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
	// declaring matrices as two dimensional vector of pointer to doubles
	//vector<vector<double>* > mat(MATRIX_SIZE, new vector<double>(MATRIX_SIZE, 0));
	vector<vector<double>* > mat, upper, lower, mat_dup;
	// declaring permutation matrix as vector of pointers to integers
	vector<int> permutation(MATRIX_SIZE);

	// seed the random number generator
	srand48(time(0));											/// seed differently for each thread ###################

	// ------------------------------------------------------------------------------------------
	// ------------- sequential timing	: 3.5 sec (two dimensional matrix pointing to doubles)
	// ------------- parallel timing	: 2.9 sec (two dimensional matrix pointing to doubles)
	// ------------------------------------------------------------------------------------------
	// intializing all four matrices
	#pragma omp parallel num_threads(NUM_THREADS)
	{
		#pragma omp sections
		{
			#pragma omp section
			{
				for(int i=0; i<MATRIX_SIZE; i++)
				{
					mat.push_back(new vector<double>(MATRIX_SIZE, (0)));
					mat_dup.push_back(new vector<double>(MATRIX_SIZE, (0)));
				}
			}
			#pragma omp section
			{
				for(int i=0; i<MATRIX_SIZE; i++)
					upper.push_back(new vector<double>(MATRIX_SIZE, (0)));
			}
			#pragma omp section
			{
				for(int i=0; i<MATRIX_SIZE; i++)
				{
					lower.push_back(new vector<double>(MATRIX_SIZE, (0)));
					(*lower[i])[i] = (1);
				}
			}
			#pragma omp section
			{
				for(int i=0; i<MATRIX_SIZE; i++)
					permutation[i] = i;
			}
		}
		// Initializing mat: 64000000 elements with random doubleing point values
		#pragma omp for schedule(static, 4) collapse(2)
		for(int i=0; i<MATRIX_SIZE; i++)
		{
			for(int j=0; j<MATRIX_SIZE; j++)
			{
				double random_num = (drand48());
				(*mat[i])[j] = random_num;
				(*mat_dup[i])[j] = random_num;
			}
		}
		// --------------------- only for debugging ----------------------------
	//	(*mat[0])[0] = 1; (*mat_dup[0])[0] = 1;
	//	(*mat[0])[1] = 2; (*mat_dup[0])[1] = 2;
	//	(*mat[0])[2] = 3; (*mat_dup[0])[2] = 3;
	//	(*mat[1])[0] = 3; (*mat_dup[1])[0] = 3;
	//	(*mat[1])[1] = 1; (*mat_dup[1])[1] = 1;
	//	(*mat[1])[2] = 2; (*mat_dup[1])[2] = 2;
	//	(*mat[2])[0] = 2; (*mat_dup[2])[0] = 2;
	//	(*mat[2])[1] = 3; (*mat_dup[2])[1] = 3;
	//	(*mat[2])[2] = 1; (*mat_dup[2])[2] = 1;
		// ---------------------------------------------------------------------
	}

//	double start_omp = omp_get_wtime();

	// outermost loop: non parallelizable due to data dependecies across iterations
	for(int k=0; k< MATRIX_SIZE; k++)
	{
		pair<int, double> max_of_column;
		max_of_column = get_maximum_element_index(mat, k, k, MATRIX_SIZE-1);
		int k_dash = max_of_column.first;
		double max_element = max_of_column.second;
		// if singular matrix then stop
		if(k_dash == -1)	{ cout << "Singular matrix.\n"; return 0;}
		
		// swapping k and k_dash row of matrix mat
		// This operation can not be parallelized because second last for loop assumes
		// two rows are swapped. If we do it parallely then second last for loop may
		// read non-swapped elements.
		vector<double>* temp_row_pointer = mat[k];
		mat[k] = mat[k_dash];
		mat[k_dash] = temp_row_pointer;
		// coping first (k-1) elements from kth row of lower triangular matrix because
		// while swapping elements in the loop every time going to that position through
		// pointers and copy the content into some temporary variable is time consuming
		// instead copy the whole subvector content to a temp vector and use that in 
		// swapping elements. This will reduce order(n^2) instruction.
		auto start = (*lower[k]).begin();
		auto end = (*lower[k]).begin()+k;
		vector<double> temp_sub_vector_lower_matrix(start, end);

		start = (*lower[k_dash]).begin();
		end = (*lower[k_dash]).begin() + k;
		vector<double> temp_sub_vector_lower_matrix_kdash(start, end);

		vector<double> temp_sub_vector_upper_kth_row(MATRIX_SIZE - 1 - k);
		start = (*mat[k]).begin() + (k+1);
		end = (*mat[k]).end();
		vector<double> temp_sub_vector_mat_kth_row(start, end);
		// starting a parallel section
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			// only master threads completes the following work
			#pragma omp master
			{
				// swapping k and k_dash elements of permutation vector
				int temp = permutation[k];
				permutation[k] = permutation[k_dash];
				permutation[k_dash] = temp;
				// updating upper triangular matrix's (k, k) entry
				(*upper[k])[k] = max_element;
			}
			// swapping lower triangular matrix elements
			#pragma omp for schedule(static, 4) nowait
			for(int j=0; j<=k-1; j++)
			{
				(*lower[k])[j] = temp_sub_vector_lower_matrix_kdash[j];
				(*lower[k_dash])[j] = temp_sub_vector_lower_matrix[j];
			}
			// updating fractions in lower triangular matrix
			// and elements of upper triangular matrix
			#pragma omp for schedule(static, 4)
			for(int j=k+1; j<MATRIX_SIZE; j++)
			{
				(*lower[j])[k] = (*mat[j])[k]/max_element;
				(*upper[k])[j] = temp_sub_vector_mat_kth_row[j - (k+1)];
			}
			// --------------------------------------------------------- barrier ------------------------------------------------------
			#pragma omp single
			{
				// coping kth row of upper matrix: constant in next for loop
				auto start = (*upper[k]).begin() + k+1;
				auto end = (*upper[k]).end();
				copy(start, end, temp_sub_vector_upper_kth_row.begin());
			}
			// --------------------------------------------------------- barrier -------------------------------------------------------
			// finally updating the input matrix
			#pragma omp for schedule(static, 4)
			for(int i=k+1; i<MATRIX_SIZE; i++)
			{
				double operand1 = (*lower[i])[k];
				for(int j=k+1; j<MATRIX_SIZE; j++)
					(*mat[i])[j] -= ((operand1 * temp_sub_vector_upper_kth_row[j-(k+1)]));
			}
		}
	}

//	double end_omp = omp_get_wtime();
//	double duration_omp = (double)(end_omp -start_omp);
//	cout << duration_omp << endl;

	// getting the end time
	auto end = high_resolution_clock::now();
	// get the duration
	auto duration = duration_cast<microseconds>(end - start);
	cout << "Time taken: " << duration.count() << " micro seconds" << endl;

//	print_matrix(lower, MATRIX_SIZE, MATRIX_SIZE);
//	print_matrix(upper, MATRIX_SIZE, MATRIX_SIZE);
//	print_vector(permutation, MATRIX_SIZE);
//	print_matrix(mat, MATRIX_SIZE, MATRIX_SIZE);

	return 0;
}
// This function takes original input, permutation matrix, lower and upper triangular matrix as input and calculates the error
// First rearranges the original input matrix a/c to permutation matrix. Then multiplies lower and upper triangular matrix and 
// subtracts the result with rearranged matrix. and prints the error.
void print_error(vector<vector<double>*> &mat_dup, vector<int> &permutation, vector<vector<double>*> &lower, vector<vector<double>*> &upper)
{

}

// This function finds the maximum element from a fix column and between raw index r1 nad r2
// If such non-zero element exists then it returns its (index, max_value)
// Else it returns (-1, 0).
// While comparing absolute values is compared against each other
// ------------------------------------------------------------------------------------------------------
// ------------- sequential time:	4.9 seconds (two dimensional matrix pointing to doubles)
// ------------- parallel time:		4.9 seconds (two dinensional matrix pointing to doubles)
// ------------- parallel time: 	2.2 seconds (one dimensional vector pointing to vector of doubles)
// ------------------------------------------------------------------------------------------------------
pair<int, double> get_maximum_element_index(vector<vector<double>* > &mat, int c, int r1, int r2)
{
	double max_element=0;
	int index=-1;
	for(int i= r1; i<= r2; i++)
	{
		if(abs((*mat[i])[c]) > max_element)
		{
			max_element = abs((*mat[i])[c]);
			index = i;
		}
	}
	return make_pair(index, max_element);
}

void print_matrix(vector<vector<double>* > mat, int r, int c)
{
	for(int i=0; i<r; i++)
	{
		for(int j=0; j<c; j++)	cout << (*mat[i])[j] << " ";
		cout << endl;
	}
	cout << endl;
}
void print_vector(vector<int> permutation, int size)
{
	for(int i=0; i<size; i++)	cout << permutation[i] << " ";
	cout << endl << endl;
}
