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
			#pragma omp for reduction(max: max_element)
			for(int i=0; i<vec_size; i++)	
				max_element = max_element > vec[i] ? max_element : vec[i];
		}
		cout << "max_element: " << max_element << endl;
	}
}

int main()
{
	int vec_size;
	cin >> vec_size >> NUM_THREADS;
	// getting the current time
	auto start = high_resolution_clock::now();
	
	reduction(vec_size, 0);
	
	// getting the end time
	auto end = high_resolution_clock::now();
	// get the duration
	auto duration = duration_cast<microseconds>(end - start);
	cout << "seq time: " << duration.count() << endl;
	cout << endl;
//	//-------------------------------------------------------------------------------------------------------------------------
	// getting the current time
	start = high_resolution_clock::now();
	
	reduction(vec_size, 1);
	
	// getting the end time
	end = high_resolution_clock::now();
	// get the duration
	duration = duration_cast<microseconds>(end - start);
	cout << "reduction time: " << duration.count() << endl;
	cout << endl;
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
