#include<iostream>
#include<omp.h>

using namespace std;

int main()
{
	int size;
	cin >> size;
	double arr[size][size];
//	#pragma omp parallel for num_threads(4) schedule(static, 4) collapse(2)
	for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			arr[i][j] = i+j;
	return 0;
}
