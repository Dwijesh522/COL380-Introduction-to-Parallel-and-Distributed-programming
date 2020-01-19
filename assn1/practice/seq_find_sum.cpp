#include<iostream>
using namespace std;
int main()
{
	int size = 1000000000;
	int *arr = new int[size];
	for(int i=0; i<size; i++)	arr[i] = i;
	int sum=0;
	for(int i=0; i<size; i++)	sum += arr[i];
	cout << "sum of array elements: " << sum << endl;
	return 0;
}
