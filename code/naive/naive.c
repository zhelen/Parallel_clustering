/*******************************************************************************
 *	EECS587 Final                                                   *
 *	Ze Jia Zhang (3184 7771)                                               *
 *                                                                             *
 *  Compiler instructions:                                                     *
 *  mpicc naive.c -O3 -lm -o naive                                                     *
 *                                                                             *
 *  Run instructions:                                                          *
 *  mpirun -np <number of processors> naive <input file>                  *
 *                                                                             *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

#define SCATTER_X 0
#define SCATTER_N 1
#define SCATTER_D 2

double distance(double ** X, int j, int k, int d);
double ** allocate_double(int r, int c);

int main (int argc, char * argv[]) {
	// declare variables
	double ** X; 
	double * hs;
	double h1, beta, dist;
	int * netSize;
	int * numCentres;
	int * indCentres;
	int * numNb;
	int * indNb;
	int * isNb;
	int * isCentre;
	int i, j, k, noX, tmp;
	int N, N0, d, netDim, S;
	int currCentre,currProc;
	FILE *dataFile, *listC, *numC;
	
	int pid,np;
	int Nlocal,p,noXglobal;
	int * numGlobal;
	int * procDone;
	double * centreGlobal;
	int prow, pcol; //location of proc in mesh

	double timeStart, timeEnd; // for timing the program
	double timeTotal;
	
	MPI_Status status;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&np);
	MPI_Comm_rank(MPI_COMM_WORLD,&pid);

	timeStart = clock();
	// parse input
	if ( argc != 2 ) {
		if ( pid == 0 ){
        	printf("Error: Enter input file \n");
        }
		return 1;
    } else {
    	if ( pid == 0 ){
    		timeStart = MPI_Wtime();
			dataFile = fopen(argv[1], "r");
			if (dataFile == NULL) {
        		printf("Error: Invalid input file\n");
        		return 1;
    		}
    		tmp = fscanf(dataFile, "%d %d %d",&N, &d, &netDim);
    		

						netSize = malloc(netDim*sizeof(int));
						
			for (i = 0; i < netDim; i++){
				tmp = fscanf(dataFile, "%d", &netSize[i]);	
			}
			tmp = fscanf(dataFile, "%lf %lf %d",&h1, &beta, &S);
			
			Nlocal = N/np + N%np;
			X = allocate_double(Nlocal,d);
/*			for (i = 0; i < N; i++) {
				X[i]=&(X[0][i*d]);
			}*/
			
			Nlocal = N/np;
//			MPI_Bcast(&Nlocal,1,MPI_INT,0,MPI_COMM_WORLD);

			for (p = 1;p<np;p++) {		
//			printf("%d %d\n",pid,p);	
				MPI_Send(&Nlocal,1,MPI_INT,p,SCATTER_N,MPI_COMM_WORLD);
				MPI_Send(&d,1,MPI_INT,p,SCATTER_D,MPI_COMM_WORLD);
				for (i = 0; i < Nlocal; i++){
					for (j = 0; j < d; j++){
						tmp = fscanf(dataFile, "%lf", &X[i][j]);	
						if (tmp < 1){
							printf("Error: Incorrect input file\n");
							return 1;
						}
					}
				}
//				printf("%d %d %d\n",pid,Nlocal,p);
				MPI_Send(&(X[0][0]),Nlocal*d,MPI_DOUBLE,p,SCATTER_X,
					MPI_COMM_WORLD);
			}
			
			Nlocal = N/np + N%np;
//			printf("%d %d %d\n",pid,Nlocal,d);
			for (i = 0; i < Nlocal; i++){
				for (j = 0; j < d; j++){
					tmp = fscanf(dataFile, "%lf", &X[i][j]);
//					printf("%d %d %lf\n",pid,Nlocal,X[i][j]);	
					if (tmp < 1){
						printf("Error: Incorrect input file\n");
						return 1;
					}
				}
			}
			N0 = N/np*(np-1);
			numGlobal = malloc(S*sizeof(int));
			procDone = malloc(np*sizeof(int));
//			printf("%d %d %d\n",pid,Nlocal,123);
			
	//		fclose(dataFile);

			printf("Beginning computation for %d points of dimension %d ",N,d);		
			printf("on a %d dimensional mesh of %d processors...\n",netDim,np);
    	}
//    	printf("%d %d %d\n",pid,Nlocal,130);
	}

	if (pid > 0) {
	
		MPI_Recv(&Nlocal,1,MPI_INT,0,SCATTER_N,MPI_COMM_WORLD,&status);
		MPI_Recv(&d,1,MPI_INT,0,SCATTER_D,MPI_COMM_WORLD,&status);
//		printf("%d %d %d\n",pid,Nlocal,d);
//		MPI_Probe(0, SCATTER_X, MPI_COMM_WORLD, &status);
//		MPI_Get_count(&status, MPI_DOUBLE, &n);
//		d = n/Nlocal;
//printf("%d %d %d\n",pid,Nlocal,n);
		X = allocate_double(Nlocal,d);
/*		for (i = 0; i < N; i++) {
			X[i]=&(X[0][i*d]);
		}*/
//printf("%d %d %d %lf\n",pid,Nlocal,d,X[0][0]);
		MPI_Recv(&(X[0][0]),Nlocal*d,MPI_DOUBLE,0,SCATTER_X,MPI_COMM_WORLD,
			&status);
//					printf("%d %d %d\n",pid,Nlocal,d);
		N0 = (pid-1)*Nlocal;
	}

	MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&h1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&beta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&S,1,MPI_INT,0,MPI_COMM_WORLD);

	indCentres = malloc(Nlocal*sizeof(int));
	numNb = malloc(Nlocal*sizeof(int));
	indNb = malloc(Nlocal*sizeof(int));
	isNb = malloc(Nlocal*sizeof(int));
	isCentre = malloc(Nlocal*sizeof(int));
//			netSize = malloc(netDim*sizeof(int));
			
    hs = malloc(S*sizeof(double));
	numCentres = malloc(S*sizeof(int));
	
	centreGlobal = malloc(d*sizeof(double));
	
/*printf("%d %d %d\n",pid,Nlocal,167);
	
	for (i = 0; i < Nlocal; i++){
        for (j = 0; j < d; j++){
    		printf("%d %lf ",i,X[i][j]);	
        }
        printf("\n");	
    }*/


	memset(isCentre, 0, Nlocal*sizeof(int));
	noX = 0;
    for (i = 0; i < S; i++){
        hs[i] = h1*pow(beta,i);
        if (i > 0){
        	numCentres[i] = numCentres[i-1];
    	} else {
			numCentres[i] = 0;
		}
    	
//    	printf("%lf \n",hs[i]);	
		memset(isNb, 0, Nlocal*sizeof(int));
		noX = 0;
		noXglobal = 0;
		currProc = 0;
		while (noXglobal == 0){
		
			MPI_Gather(&noX, 1, MPI_INT, procDone, 1, MPI_INT, 0, 
			MPI_COMM_WORLD);
			
			if(pid==0) {
				noXglobal = 1;
				for (p = 0;p<np;p++) {		
					if(procDone[p] == 0){
						noXglobal = 0;
						currProc = p;
						break;
					}
				}
			}
			
			MPI_Bcast(&noXglobal,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(&currProc,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			
			if(pid==currProc){
				for (j=0;j<Nlocal;j++){
					if(isCentre[j]==0 && isNb[j]==0){
						isCentre[j] = 1;
						tmp = numCentres[i];
						indCentres[tmp] = j+N0;
					}
				}
			}
		
					for (k=0;k<N;k++){
						if(isCentre[k]==0 && isNb[k]==0){
							dist = distance(X,j,k,d);
//							printf("%d %d %lf\n",j,k,dist);
							if(dist < hs[i]){
								isNb[k] = 1;
//								printf("%lf %lf\n",dist,hs[i]);
//								indNb[numNb[numCentres[i]]] = k;
//								numNb[numCentres[i]] += 1;
//								printf("%d \n",numNb[numCentres[i]]);
							}
						}
					}
					numCentres[i] += 1;

//				if (noX != 1) {
//					numNb[numCentres[i]] = numNb[numCentres[i]-1];
//				}
			}		
		}	
//		MPI_Barrier(MPI_COMM_WORLD);		    					
    } 
//    MPI_Reduce(numCentres, numGlobal, S, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
/*    for (i=0;i<numCentres[S-1];i++){
    	printf("%d \n",indCentres[i]);	
    }

	printf("%d %d %d \n",numCentres[0],numCentres[1],numCentres[2]);
	printf("%lf \n",hs[0]);*/
	
/*	numC = fopen("numcentres.out", "w");
	listC = fopen("listcentres.out", "w");

	if (numC == NULL || listC == NULL) {
  		printf("Error: Cannot create output file\n");
  		return 1;
	}
	
	for (i=0;i<numCentres[S-1];i++){
    	fprintf(listC,"%d \n",indCentres[i]);	
    }
    
    for (i=0;i<S;i++){
    	fprintf(numC,"%d \n",numCentres[i]);	
    }
    
    fclose(listC);
    fclose(numC);
    fclose(dataFile);*/

	if(pid==0){
		timeEnd = MPI_Wtime();
		timeTotal = timeEnd - timeStart;
		printf("Finished.\n");
		printf("Elapsed time: %.6lf s.\n",timeTotal);
	}
	
	
    free(X[0]);
	free(X);
	free(hs);
//	free(netSize);
	free(indCentres);
	free(numCentres);
	free(indNb);
	free(numNb);
	free(isNb);
	free(isCentre);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	
	return 0;
}

double distance(double ** X, int j, int k, int d) {
    int D;
    double dists = 0.0;
    for (D = 0;D<d;D++){
    	dists = dists + (X[j][D]-X[k][D])*(X[j][D]-X[k][D]);
//    	printf("%lf %lf ",X[j][D],X[k][D]);
	}
	dists = sqrt(dists);
//	printf("%lf \n",dists);
	return dists;
}

double ** allocate_double(int r, int c) {
	int i;
    double * data = (double *)malloc(r*c*sizeof(double));
    double ** array= (double **)malloc(r*sizeof(double*));
    for (i=0; i<r; i++)
        array[i] = &(data[c*i]);
    return array;
}
