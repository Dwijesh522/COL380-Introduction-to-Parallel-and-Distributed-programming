#include<vector>
#include<fstream>
#include<iostream>
#include<string.h>
#include <boost/algorithm/string.hpp>
using namespace std;

vector<vector<double>*> read_matrix(string filename, int n)
{
	vector<vector<double>*> mat;
	// initializing first verticle vector of matrix
	for(int i=0; i<n; i++)	mat.push_back(new vector<double>(n));
	int matrix_counter = 0;

	string line;
	fstream matrix_file(filename);
	if(matrix_file.is_open())
	{
		while(getline(matrix_file, line))
		{
			vector<string> string_elements;
			boost::split(string_elements, line, boost::is_any_of(" "));
			for(int i=0; i<n; i++)	(*mat[matrix_counter])[i] = stod(string_elements[i]);
			matrix_counter++;
		}
		matrix_file.close();
	}
	else	cout << "error while opening the file\n";

	return mat;
}


void output_files_omp(vector<vector<double>* > mat, int r, int c, string file_name)
{
  ofstream file;
  file.open(file_name);
  for(int i=0; i<r; i++)
{
for(int j=0; j<c; j++) file << (*mat[i])[j] << " ";
file << endl;
}
file.close();
}

void output_P_omp(vector<int> P, int size, string file_name)
{
  ofstream file;
  file.open(file_name);
  for(int i=0; i<size; i++) file << P[i] << " ";
  file.close();
}

void output_files_pthread(double** mat, int r, int c, string file_name)
{
  ofstream file;
  file.open(file_name);
  for(int i=0; i<r; i++)
{
for(int j=0; j<c; j++) file << mat[i][j] << " ";
file << endl;
}
file.close();
}

void output_P_pthread(int* P, int r, int c, string file_name)
{
  ofstream file;
  file.open(file_name);
  for(int i=0; i<size; i++) file << P[i] << " ";
  file.close();
}
