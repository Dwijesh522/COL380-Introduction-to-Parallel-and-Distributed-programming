#include<iostream>
#include<omp.h>
#include<chrono>

using namespace std;
using namespace std::chrono;

float fun1(int size)
{
	float a=1, b=1;
	for(int i=0; i<size; i++)
	{
		b = (b*1000/978)/98;
		a = (a*1000/978)/98*b;
	}	
	return a;
}

float fun2(int size)
{
	float a=1, b=1;
	#pragma omp parallel for num_threads(1) //schedule(static, 4)
	for(int i=0; i<size; i++)
	{
		b = (b*1000/978)/98;
		a = (a*1000/978)/98*b;
	}
	return a;
}
int main()
{
	
	int size, id;
	cin >> size >> id;
	// getting the current time
	auto start = high_resolution_clock::now();
	if(id == 0) 	fun1(size);
	else		fun2(size);

	// getting the end time
	auto end = high_resolution_clock::now();
	// get the duration
	auto duration = duration_cast<microseconds>(end - start);
	cout << "Time taken: " << duration.count() << " micro seconds" << endl;
	return 0;
}
