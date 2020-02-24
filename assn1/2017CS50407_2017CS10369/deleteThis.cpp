#include "read_write.h"
#include<iostream>
#include<vector>

using namespace std;

int main()
{
	vector<vector<double>*> mat = read_matrix("B_100.txt", 100);
	return 0;
}
