#include<stdio.h>

int main()
{
	#pragma omp parallel
	{
		printf("hello world\n");
	}
	return 0;
}
