/*******************************************************************************
 *	EECS587 Problem Set 3                                                  *
 *	Introduction to MPI                                                    *
 *	Ze Jia Zhang (3184 7771)                                               *
 *                                                                             *
 *  Compiler instructions:                                                     *
 *  mpicc -O3 -o PS3 ps3.c                                                     *
 *                                                                             *
 *  Run instructions:                                                          *
 *  mpirun -np <number of processors> PS3 <size of matrix, n>                  *
 *                                                                             *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main (int argc, char * argv[]) {
	// declare variables
	double ** A; // local A matrix and previous iteration matrix
	double ** A0; 
	double * A_t; // top halo cells
	double * A_r; // right halo cells
	double * A_l; // left halo cells for sending (column->row)
	double * A_b; // bottom halo cells
	
	int pid, np, nn; // current id, total number of procs, matrix size
	int i, j, i_global, j_global; // local and global indices for loops
	int i_max, j_max; // max dimension of local matrix
	int np_root; // square root of np
	int p_row, p_col; //location of proc in mesh

	int t, k; // iteration counters
	int t_max = 10; // maximum iterations

	double t_start, t_end; // for timing the program
	double sum, min, sum_global, min_global; // verification values

	double A_ip, A_ipp, A_jp, x; // to store A(i',j), etc.

	MPI_Status status;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&np);
	MPI_Comm_rank(MPI_COMM_WORLD,&pid);

	// parse input
	if ( argc != 2 ) {
		if ( pid == 0 ){
        		printf("Error: Enter size of matrix n \n");
		}
		return 1;
        } else {
		sscanf(argv[1], "%d", &nn);
	}

	// check input
	if ( np > pow((double)nn/3,2) ) {
		if ( pid == 0 ){
        		printf("Error: p too large \n");
		}
		return 1;
        }

	np_root = (int)sqrt(np);
	if ( pid == 0 ){
		printf("Beginning computation for %d by %d matrix ",nn,nn);		
		printf("on a mesh of %d by %d processors...\n",
			np_root,np_root);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// compute matrix size
	p_row = (int)floor(pid/np_root);
	p_col = (int)pid - p_row*np_root; 
	if ( p_col == np_root-1 ) {
		j_max = nn/np_root + nn%np_root;
	} else {
		j_max = nn/np_root;
	}
	if ( p_row == np_root-1 ) {
		i_max = nn/np_root + nn%np_root;
	} else {
		i_max = nn/np_root;
	}

	// allocate memory
	A = malloc((i_max)*sizeof(double *));
	A0 = malloc((i_max)*sizeof(double *));
	A[0] = malloc(i_max*j_max*sizeof(double));
	A0[0] = malloc(i_max*j_max*sizeof(double));
	for (i = 0; i < i_max; i++) {
		A[i]=&(A[0][i*j_max]);
		A0[i]=&(A0[0][i*j_max]);
	}
	A_t = malloc((j_max)*sizeof(double));
	A_b = malloc((j_max)*sizeof(double));
	A_r = malloc((i_max)*sizeof(double));
	A_l = malloc((i_max)*sizeof(double));

	// initialise matrix
	for (i = 0; i < i_max; i++) {
		i_global = p_row*(nn/np_root) + i;
		for (j = 0; j < j_max; j++) {
			j_global = p_col*(nn/np_root) + j;
			if (i_global == j_global) {
				A[i][j] = i_global*sin(sqrt(i_global));
			} else {
				A[i][j] = pow(i_global+j_global,1.1);
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	// start timing
	if ( pid == 0 ){
		t_start = MPI_Wtime();
	}

	for (t = 0; t < t_max; t++) {
		// back up matrix
		for (i = 0; i < i_max; i++) {
			for (j = 0; j < j_max; j++) {
				A0[i][j] = A[i][j];
			}
		}

		// construct leftmost column into a row
		if (p_col > 0) {
			for (i = 0; i < i_max; i++) {
				A_l[i] = A[i][0];
			}
			// send left; messages are tagged "1"
			MPI_Send(A_l, i_max, MPI_DOUBLE, pid-1, 1, MPI_COMM_WORLD);
		}
		if (p_col < np_root-1) { // not cyclical
			// receive from right; messages are tagged "1"
			MPI_Recv(A_r, i_max, 
				MPI_DOUBLE, pid+1, 1, MPI_COMM_WORLD,&status);
		}
		// send top; messages are tagged "2"
		if (p_row > 0) {
			MPI_Send(A[0], j_max, MPI_DOUBLE, pid-np_root, 2, MPI_COMM_WORLD);
		} else {
			MPI_Send(A[0], j_max, MPI_DOUBLE, 
				np_root*(np_root-1) + pid, 2, MPI_COMM_WORLD);
		}
		// receive from bottom; messages are tagged "2"
		if (p_row < np_root-1) {
			MPI_Recv(A_b, j_max, 
				MPI_DOUBLE, pid+np_root, 2, MPI_COMM_WORLD,&status);
		} else {
			MPI_Recv(A_b, j_max, 
				MPI_DOUBLE, pid%np_root, 2, MPI_COMM_WORLD,&status);
		}
		// send bottom; messages are tagged "3"
		if (p_row < np_root-1) {
			MPI_Send(A[i_max-1], j_max, 
				MPI_DOUBLE, pid+np_root, 3, MPI_COMM_WORLD);
		} else {
			MPI_Send(A[i_max-1], j_max, 
				MPI_DOUBLE, pid%np_root, 3, MPI_COMM_WORLD);
		}
		// receive from top; messages are tagged "3"
		if (p_row > 0) {
			MPI_Recv(A_t, j_max, 
				MPI_DOUBLE, pid-np_root, 3, MPI_COMM_WORLD,&status);
		} else {
			MPI_Recv(A_t, j_max, MPI_DOUBLE, 
				np_root*(np_root-1) + pid, 3, MPI_COMM_WORLD,&status);
		}
		
		// all halos are ready, now update A
		for (i = 0; i < i_max; i++) {
			i_global = p_row*(nn/np_root) + i;
			for (j = 0; j < j_max; j++) {
				j_global = p_col*(nn/np_root) + j;
				// update appropriate cells only
				if (i_global != j_global && j_global != nn-1) {
					if (i < i_max-1) {
						A_ip = A0[i+1][j];
					} else {
						A_ip = A_b[j];					
					}
					if (i > 0) {
						A_ipp = A0[i-1][j];
					} else {
						A_ipp = A_t[j];
					}
					if (j < j_max-1) {		
						A_jp = A0[i][j+1];
					} else {
						A_jp = A_r[i];					
					}
					x = 0.0;
					for (k = 1; k <= 10; k++) {
						x = x + pow(fabs(0.5+A_ip),1.0/k) - 
							pow(fabs(A_ipp),1.0/(k+1))*
							pow(fabs(A_jp),1.0/(k+2));
					}
					A[i][j] = fmax(-10.0,fmin(10.0,x));
				}
			}
		}
	}

	// compute local sum and min
	min = A[0][0];
	sum = 0.0;
	for (i = 0; i < i_max; i++) {
		for (j = 0; j < j_max; j++) {
			sum += fabs(A[i][j]);
			if (A[i][j] < min) {
				min = A[i][j];
			}
		}
	}

	// reduce
	MPI_Reduce(&sum, &sum_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&min, &min_global, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

	// stop timing
	if ( pid == 0 ){
		t_end = MPI_Wtime();
		printf("Elapsed time is %.6lf seconds\n", t_end-t_start );
		printf("Sum of A and minimum in A are below\n");
		printf("%lf \n",sum_global);
		printf("%lf \n",min_global);
	}
    free(A[0]);
	free(A);
    free(A0[0]);
	free(A0);
	free(A_t);
	free(A_b);
	free(A_r);
	free(A_l);
 	MPI_Finalize();
	return 0;
}
