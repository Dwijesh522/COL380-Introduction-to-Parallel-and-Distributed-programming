// omp header
#include<omp.h>
#include<iostream>
#include<vector>
#include<chrono>

using namespace std;
using namespace std::chrono;

// CONSTANTS
// Default number of threads are 4
int NUM_THREADS = 4;

// create a vector of pointers with ith element pointing to value i
vector<int *> vector_of_pointers(int n)
{
	vector<int *> pi(n);
	for(int i=0; i<n; i++)
		pi[i] = new int(i);
	return pi;
}

// Atomic
void atomic()
{
	int count = 0;
	#pragma omp parallel num_threads(NUM_THREADS)
	{
		#pragma omp atomic
		count++;
	}
	cout << "count: " << count << endl;
}

// heirarchical pragma blocks
void heirarchical_pragmas()
{
	#pragma omp parallel num_threads(NUM_THREADS)
	{
		cout << "hi\n";
		#pragma omp barrier
		cout << "hey\n";
	}
}

// Critical
// parallelism matters for size > 10^8
void critical(int id, int size)
{
	// initializing the vector
	vector<int> vect;
	for(int i=0; i<size; i++) vect.push_back((i)%(size/2));
	// finding the maximum element sequentially
	if(id == 0)
	{
		int max=-1;
		for(int i=0; i<size; i++) if(vect[i] > max) max = vect[i];
		cout << "max element: " << max << endl;
	}
	// parallel but no scheduling
	else if(id == 1)
	{
		int max=-1;
		int thread_max=-1;
		#pragma omp parallel num_threads(NUM_THREADS) firstprivate(thread_max)
		{
			// four threads find their local maximas
			#pragma omp for
			for(int i=0; i<size; i++)
				if(vect[i] > thread_max) thread_max = vect[i];
			// updating the global maxima
			#pragma omp critical
			{
				if(thread_max > max) max = thread_max;
			}
		}
		cout << "max element: " << max << endl;
	}
	// parallel with static scheduling without chunk size
	else if(id == 2)
	{
		int max=-1;
		int thread_max=-1;
		#pragma omp parallel num_threads(NUM_THREADS) firstprivate(thread_max)
		{
			// four threads find their local maximas
			#pragma omp for schedule(static)
			for(int i=0; i<size; i++)
				if(vect[i] > thread_max) thread_max = vect[i];
			// updating the global maxima
			#pragma omp critical
			{
				if(thread_max > max) max = thread_max;
			}
		}
		cout << "max element: " << max << endl;
	}
	// parallel with static scheduling with chunk size
	else if(id == 3)
	{
		int max=-1;
		int thread_max=-1;
		#pragma omp parallel num_threads(NUM_THREADS) firstprivate(thread_max)
		{
			// four threads find their local maximas
			#pragma omp for schedule(static, 500)
			for(int i=0; i<size; i++)
				if(vect[i] > thread_max) thread_max = vect[i];
			// updating the global maxima
			#pragma omp critical
			{
				if(thread_max > max) max = thread_max;
			}
		}
		cout << "max element: " << max << endl;
	}
	// parallel with dynamic scheduling without chunk size:1
	else if(id == 4)
	{
		int max=-1;
		int thread_max=-1;
		#pragma omp parallel num_threads(NUM_THREADS) firstprivate(thread_max)
		{
			// four threads find their local maximas
			#pragma omp for schedule(dynamic)
			for(int i=0; i<size; i++)
				if(vect[i] > thread_max) thread_max = vect[i];
			// updating the global maxima
			#pragma omp critical
			{
				if(thread_max > max) max = thread_max;
			}
		}
		cout << "max element: " << max << endl;
	}
	// parallel with dynamic scheduling with chunk size:NUM_THREADS
	else if(id == 5)
	{
		int max=-1;
		int thread_max=-1;
		#pragma omp parallel num_threads(NUM_THREADS) firstprivate(thread_max)
		{
			// four threads find their local maximas
			#pragma omp for schedule(dynamic, 500)
			for(int i=0; i<size; i++)
				if(vect[i] > thread_max) thread_max = vect[i];
			// updating the global maxima
			#pragma omp critical
			{
				if(thread_max > max) max = thread_max;
			}
		}
		cout << "max element: " << max << endl;
	}
}

// reductions
void reduction(int vec_size, int id)
{
	// initialization
	int vec[vec_size];
	for(int i=0; i<vec_size; i++) vec[i]=i+1;
	// seq code
	if(id == 0)
	{
		int max_element=-1;
		for(int i=0; i<vec_size; i++)
			max_element = max_element > vec[i] ? max_element : vec[i];
		cout << "max_element: " << max_element << endl;
	}
	// par code with reduction applied
	else
	{
		int max_element=-1;
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			#pragma omp for schedule(static, NUM_THREADS) reduction(max: max_element)
			for(int i=0; i<vec_size; i++)	
				max_element = max_element > vec[i] ? max_element : vec[i];
		}
		cout << "max_element: " << max_element << endl;
	}
}

// printing a two dimensional vector
void print_matrix(vector<vector<int>> mat, int r, int c)
{
	for(int i=0; i<r; i++)
	{
		for(int j=0; j<c; j++)	cout << mat[i][j] << " ";
		cout << endl;
	}
}

// collapse
void collapse(int vec_size, int id)
{	
	// sequential execution
	if(id == 0)
	{
		// initialization
		vector<vector<int>> a, b, c;
		for(int i=0; i<vec_size; i++)
		{
			vector<int> temp_vec;
			for(int j=0; j<vec_size; j++) temp_vec.push_back(i);
			a.push_back(temp_vec);
			b.push_back(temp_vec);
			c.push_back(temp_vec);
		}
		// modifing matrix a
		for(int i=0; i<vec_size; i++)
			for(int j=0; j<vec_size; j++)
				a[i][j] = a[i][j] - b[i][j] - c[i][j];
	}
	else
	{
		// initialization
		vector<vector<int>> a(vec_size, vector<int>(vec_size));
		vector<vector<int>> b(vec_size, vector<int>(vec_size));
		vector<vector<int>> c(vec_size, vector<int>(vec_size));
		#pragma omp parallel for collapse(2) schedule(static, NUM_THREADS)
		for(int i=0; i<vec_size; i++)
			for(int j=0; j< vec_size; j++)
			{
				a[i][j] = i;
				b[i][j] = i;
				c[i][j] = i;
			}
		// modifing matrix a
		#pragma omp parallel for collapse(2) schedule(static, NUM_THREADS)
		for(int i=0; i<vec_size; i++)
			for(int j=0; j<vec_size; j++)
				a[i][j] = a[i][j] - b[i][j] - c[i][j];
	}
}

// ordered
void ordered_clause(int times, int id)
{
	// parallel code with and without ordering
	if(id)
	{
		#pragma omp parallel for ordered num_threads(NUM_THREADS)
		for(int i=0; i<times; i++)
		{
			cout << "unordered: " << i << endl;
			#pragma omp ordered
			cout << "ordered: " << i << endl;
		}
	}
}

// no wait
void no_wait(int vec_size, int id)
{
	// seq code
	if(id == 0)
	{
		int a;
		for(int i=0; i<vec_size; i++)	a = i%1000;
		for(int i=0; i<vec_size; i++)	a = i%1000;
		for(int i=0; i<vec_size; i++)	a = i%1000;
		for(int i=0; i<vec_size; i++)	a = i%1000;
	}
	// par code, default barriers
	else if(id == 1)
	{
		int a;
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			#pragma omp for
			for(int i=0; i<vec_size; i++)	a = (i)%1000;	
			#pragma omp for
			for(int i=0; i<vec_size; i++)	a = (i)%1000;
			#pragma omp for
			for(int i=0; i<vec_size; i++)	a = (i)%1000;
			#pragma omp for
			for(int i=0; i<vec_size; i++)	a = i%1000;
		}
	}
	// par code, with no wait clause
	else
	{
		int a;
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			#pragma omp for schedule(static, NUM_THREADS) nowait
			for(int i=0; i<vec_size; i++)	a = i%1000;
			#pragma omp for schedule(static, NUM_THREADS) nowait
			for(int i=0; i<vec_size; i++)	a = i%1000;
			#pragma omp for schedule(static, NUM_THREADS) nowait
			for(int i=0; i<vec_size; i++)	a = i%1000;
			#pragma omp for schedule(static, NUM_THREADS)
			for(int i=0; i<vec_size; i++)	a = i%1000;
		}
	}
}

// sections
void sections_clause(int vec_size, int id)
{
	vector<int> vec1(vec_size), vec2(vec_size), vec3(vec_size), vec4(vec_size);
	// par code, sharing across for loops
	if(id == 0)
	{
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			#pragma omp sections
			{
				#pragma omp section
				{
					for(int i=0; i<vec_size; i++)	vec1[i] = i%1000 + (i/10)%1000;
					for(int i=0; i<vec_size; i++)	vec1[i] = i%1000 + (i/10)%1000;
					for(int i=0; i<vec_size; i++)	vec1[i] = i%1000 + (i/10)%1000;
				}
				#pragma omp section
				{
					for(int i=0; i<vec_size; i++)	vec2[i] = i%1000 + (i/10)%1000;
					for(int i=0; i<vec_size; i++)	vec2[i] = i%1000 + (i/10)%1000;
					for(int i=0; i<vec_size; i++)	vec2[i] = i%1000 + (i/10)%1000;
				}
				#pragma omp section
				{
					for(int i=0; i<vec_size; i++)	vec3[i] = i%1000 + (i/10)%1000;
					for(int i=0; i<vec_size; i++)	vec3[i] = i%1000 + (i/10)%1000;
					for(int i=0; i<vec_size; i++)	vec3[i] = i%1000 + (i/10)%1000;
				}
				#pragma omp section
				{
					for(int i=0; i<vec_size; i++)	vec4[i] = i%1000 + (i/10)%1000;
					for(int i=0; i<vec_size; i++)	vec4[i] = i%1000 + (i/10)%1000;
					for(int i=0; i<vec_size; i++)	vec4[i] = i%1000 + (i/10)%1000;
				}
			}
		}
	}
	else
	{
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			#pragma omp for schedule(static, NUM_THREADS)
				for(int i=0; i<vec_size; i++)
				{
					vec1[i] = i%1000 + (i/10)%1000;
					vec1[i] = i%1000 + (i/10)%1000;
					vec1[i] = i%1000 + (i/10)%1000;
				}
			#pragma omp for schedule(static, NUM_THREADS)
				for(int i=0; i<vec_size; i++)
				{
					vec2[i] = i%1000 + (i/10)%1000;
					vec2[i] = i%1000 + (i/10)%1000;
					vec2[i] = i%1000 + (i/10)%1000;
				}
			#pragma omp for schedule(static, NUM_THREADS)
				for(int i=0; i<vec_size; i++)
				{
					vec3[i] = i%1000 + (i/10)%1000;
					vec3[i] = i%1000 + (i/10)%1000;
					vec3[i] = i%1000 + (i/10)%1000;
				}
			#pragma omp for schedule(static, NUM_THREADS)
				for(int i=0; i<vec_size; i++)
				{
					vec4[i] = i%1000 + (i/10)%1000;
					vec4[i] = i%1000 + (i/10)%1000;
					vec4[i] = i%1000 + (i/10)%1000;
				}
		}
	}
}

int main()
{
	int vec_size;
	cin >> vec_size >> NUM_THREADS;
	// getting the current time
	auto start = high_resolution_clock::now();
	
	sections_clause(vec_size, 0);
	
	// getting the end time
	auto end = high_resolution_clock::now();
	// get the duration
	auto duration = duration_cast<microseconds>(end - start);
	cout << "case 1: " << duration.count() << endl;
	cout << endl;
//	//-------------------------------------------------------------------------------------------------------------------------
	// getting the current time
	start = high_resolution_clock::now();
	
	sections_clause(vec_size, 1);
	
	// getting the end time
	end = high_resolution_clock::now();
	// get the duration
	duration = duration_cast<microseconds>(end - start);
	cout << "case 2: " << duration.count() << endl;
	cout << endl;
//	//-------------------------------------------------------------------------------------------------------------------------
	// getting the current time
//	start = high_resolution_clock::now();
//	
//	no_wait(vec_size, 2);
//	
//	// getting the end time
//	end = high_resolution_clock::now();
//	// get the duration
//	duration = duration_cast<microseconds>(end - start);
//	cout << "case 3: " << duration.count() << endl;
//	cout << endl;
//	//-------------------------------------------------------------------------------------------------------------------------
//	// getting the current time
//	auto start = high_resolution_clock::now();
//	
//	critical(0, vec_size);
//	
//	// getting the end time
//	auto end = high_resolution_clock::now();
//	// get the duration
//	auto duration = duration_cast<microseconds>(end - start);
//	cout << "Seq time: " << duration.count() << endl;
//	cout << endl;
//	//--------------------------------------------------------------------------------------------------------------------------
//	// getting the current time
//	start = high_resolution_clock::now();
//	
//	critical(1, vec_size);
//	
//	// getting the end time
//	end = high_resolution_clock::now();
//	// get the duration
//	duration = duration_cast<microseconds>(end - start);
//	cout << "par without scheduling:  " << duration.count() << endl;
//	cout << endl;
//	//-------------------------------------------------------------------------------------------------------------------------
//	// getting the current time
//	start = high_resolution_clock::now();
//	
//	critical(2, vec_size);
//	
//	// getting the end time
//	end = high_resolution_clock::now();
//	// get the duration
//	duration = duration_cast<microseconds>(end - start);
//	cout << "par static scheduling with default chunksize:  " << duration.count() << endl;
//	cout << endl;
//	//-------------------------------------------------------------------------------------------------------------------------
//	// getting the current time
//	start = high_resolution_clock::now();
//	
//	critical(3, vec_size);
//	
//	// getting the end time
//	end = high_resolution_clock::now();
//	// get the duration
//	duration = duration_cast<microseconds>(end - start);
//	cout << "par static scheduling with NUM_THREADS chunksize:  " << duration.count() << endl;
//	cout << endl;
//	//-------------------------------------------------------------------------------------------------------------------------
//	// getting the current time
//	start = high_resolution_clock::now();
//	
//	critical(4, vec_size);
//	
//	// getting the end time
//	end = high_resolution_clock::now();
//	// get the duration
//	duration = duration_cast<microseconds>(end - start);
//	cout << "par dynamic scheduling with 1 chunksize:  " << duration.count() << endl;
//	cout << endl;
//	//-------------------------------------------------------------------------------------------------------------------------
//	// getting the current time
//	start = high_resolution_clock::now();
//	
//	critical(5, vec_size);
//	
//	// getting the end time
//	end = high_resolution_clock::now();
//	// get the duration
//	duration = duration_cast<microseconds>(end - start);
//	cout << "par dynamic scheduling with NUM_THREADS chunksize:  " << duration.count() << endl;
//	cout << endl;
//	//-------------------------------------------------------------------------------------------------------------------------
	
	return 0;
}
