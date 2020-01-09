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
	else
	{
		int max=-1;
		int thread_max=-1;
		#pragma omp parallel num_threads(NUM_THREADS) private(thread_max)
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
			cout << "tid " << omp_get_thread_num() << " done.\n";
		}
		cout << "max element: " << max << endl;
	}
}

int main()
{
	int id, vec_size;
	cin >> id >> vec_size >> NUM_THREADS;
	
	// getting the current time
	auto start = high_resolution_clock::now();
	
	critical(id, vec_size);
	
	// getting the end time
	auto end = high_resolution_clock::now();
	// get the duration
	auto duration = duration_cast<microseconds>(end - start);
	cout << "Time taken in microseconds: " << duration.count() << endl;
	
	return 0;
}
