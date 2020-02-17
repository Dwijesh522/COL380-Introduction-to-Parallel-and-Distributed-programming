#include<iostream>
#include<vector>
using namespace std;

int main()
{
	long long int first = 800, second=8;
	vector<vector<int>> arr(first, vector<int>(second));
//	for(int i=0; i<second; i++)
//		for(int j=0; j<first; j++)
//			arr[i][j] += i+j;
	for(int i=0; i<first; i++)
		for(int j=0; j<second; j++)
			arr[i][j] += i+j;
	return 0;
}
