#include<iostream>
#include<pthread.h>
#include <unistd.h>
#include <sys/syscall.h>
#define gettid() syscall(SYS_gettid)
using namespace std;

void *print_hello_world(void *message);
void *find_sum(void *array);

struct data
{
	int *arr, num_threads, offset, size;
	int *ans;
	data(int size, int num_threads, int offset)
	{
		arr = new int[size];
		this->size = size;
		this->offset = offset;
		this->num_threads = num_threads;
		ans = new int[num_threads];
		// assuming remainder = 0
		// initializing arr
		for(int i=0; i<size; i++)	arr[i] = i;
	}
};

int main()
{
	pthread_t thread_1, thread_2, thread_3, thread_4;
	int size, num_threads, offset;
//	cout << "Enter <size> <num_threads> <offset>\n";
//	cin >> size >> num_threads >> offset;
	size = 1000000000;
	num_threads = 4;
	offset = 4;
	data thread_data = data(size, num_threads, offset);
	data *data_pointer = &thread_data;

	pthread_create(&thread_1, NULL, find_sum, (void *)data_pointer);
	pthread_create(&thread_3, NULL, find_sum, (void *)data_pointer);
	pthread_create(&thread_4, NULL, find_sum, (void *)data_pointer);
	pthread_create(&thread_2, NULL, find_sum, (void *)data_pointer);

	pthread_join(thread_1, NULL);
	pthread_join(thread_2, NULL);
	pthread_join(thread_3, NULL);
	pthread_join(thread_4, NULL);

	cout << "Sum of all elements: " << thread_data.ans[0] + thread_data.ans[1] + thread_data.ans[2] + thread_data.ans[3]<< endl;
	exit(0);
}

void *find_sum(void *data_pointer)
{
	data* thread_data = (data *)data_pointer;
	pid_t tid = syscall(SYS_gettid);
	
	int *arr = thread_data->arr;
	int num_threads = thread_data->num_threads;
	int offset = thread_data->offset;
	int int_tid = tid%num_threads;
	int size = thread_data->size;

	// number of rows
	int x = (num_threads*offset);
	int c = size/x;
	int sum=0;
	for(int i=0; i<c; i++)
	{
		int start_index = (i*x) + (int_tid*offset);
		int end_index = start_index + offset-1;
		for(int j=start_index; j<= end_index; j++)	sum += arr[j];
	}
	thread_data->ans[int_tid] = sum;
}

void *print_hello_world(void *ptr_to_message)
{
	string message = *((string *)ptr_to_message);
	cout << message;
}
