/*******************************************************************************
 *	EECS587 Final                                                   *
 *	Ze Jia Zhang (3184 7771)                                               *
 *                                                                             *
 *  Compiler instructions:                                                     *
 *  gcc serial.c -O3 -lm -o serial                                                     *
 *                                                                             *
 *  Run instructions:                                                          *
 *  ./serial <input file>                  *
 *                                                                             *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

double distance(double ** X, int j, int k, int d);

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
	int N, d, netDim, S;
	int currCentre;
	FILE *dataFile, *listC, *numC;

	clock_t timeStart, timeEnd; // for timing the program
	double timeTotal;

	timeStart = clock();
	// parse input
	if ( argc != 2 ) {
        printf("Error: Enter input file \n");
		return 1;
    } else {
		dataFile = fopen(argv[1], "r");
	}

    if (dataFile == NULL) {
        printf("Error: Invalid input file\n");
        return 1;
    }
    tmp = fscanf(dataFile, "%d %d %d",&N, &d, &netDim);

    netSize = malloc(netDim*sizeof(int));

    indCentres = malloc(N*sizeof(int));
    numNb = malloc(N*sizeof(int));
    indNb = malloc(N*sizeof(int));
    isNb = malloc(N*sizeof(int));
    isCentre = malloc(N*sizeof(int));
    

    X = malloc(N*sizeof(double *));
    X[0] = malloc(N*d*sizeof(double));
	for (i = 0; i < N; i++) {
		X[i]=&(X[0][i*d]);
	}

    for (i = 0; i < netDim; i++){
    	tmp = fscanf(dataFile, "%d", &netSize[i]);	
    }

    tmp = fscanf(dataFile, "%lf %lf %d",&h1, &beta, &S);
    
    hs = malloc(S*sizeof(double));
    numCentres = malloc(S*sizeof(int));
        
    for (i = 0; i < N; i++){
        for (j = 0; j < d; j++){
    		tmp = fscanf(dataFile, "%lf", &X[i][j]);	
    		if (tmp < 1){
    			printf("Error: Incorrect input file\n");
        		return 1;
    		}
        }
    }

    fclose(dataFile);

//	printf("Beginning computation for %d points of dimension %d ",N,d);		
//	printf("on a %d dimensional mesh of %d processors...\n",netDim,1);

/*	
	for (i = 0; i < N; i++){
        for (j = 0; j < d; j++){
    		printf("%lf ",X[i][j]);	
        }
        printf("\n");	
    }
*/

	memset(isCentre, 0, N*sizeof(int));
	noX = 0;
    for (i = 0; i < S; i++){
        hs[i] = h1*pow(beta,i);
        if (i > 0){
        	numCentres[i] = numCentres[i-1];
    	} else {
			numCentres[i] = 0;
		}
    	
//    	printf("%lf \n",hs[i]);	
		memset(isNb, 0, N*sizeof(int));
		noX = 0;
		while (noX == 0){
			for (j=0;j<N;j++){
				noX = 1;
				if(isCentre[j]==0 && isNb[j]==0){
					noX = 0;
					isCentre[j] = 1;
					tmp = numCentres[i];
					indCentres[tmp] = j;
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
				}
//				if (noX != 1) {
//					numNb[numCentres[i]] = numNb[numCentres[i]-1];
//				}
			}		
		}	
					    					
    }
    
/*    for (i=0;i<numCentres[S-1];i++){
    	printf("%d \n",indCentres[i]);	
    }

	printf("%d %d %d \n",numCentres[0],numCentres[1],numCentres[2]);
	printf("%lf \n",hs[0]);*/
	
	numC = fopen("numcentres.out", "w");
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

	timeEnd = clock();
	timeTotal = ((double) (timeEnd - timeStart)) / CLOCKS_PER_SEC;
	printf("Finished.\n");
	printf("Elapsed time: %lf s.\n",timeTotal);
	
    free(X[0]);
	free(X);
	free(hs);
	free(netSize);
	free(indCentres);
	free(numCentres);
	free(indNb);
	free(numNb);
	free(isNb);
	free(isCentre);
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
