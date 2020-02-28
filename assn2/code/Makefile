one:
	mpicc blocking_point_to_point.c -lm
	mpirun -n 4 ./a.out
two:
	mpicc non_blocking_point_to_point.c -lm
	mpirun -n 4 ./a.out
three:
	mpicc collective_communication.c
	mpirun -n 4 ./a.out
