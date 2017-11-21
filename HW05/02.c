#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	int i, width, partialto, partial_remain;
    double N, totalsum = 0, partialsum = 0;
    time_t Start, End;
    double *vec = NULL;
    double *totalvec = NULL;
    double *partialvec = NULL;
    if(rank == 0)
    {
        int power;
        fflush(stdin);
        printf("It Will Compute N = 10^(power), Input power\n");
        scanf("%d",&power);
        N = pow(10.0,(double)power);
        printf("N = 10^%d = %.1lf\n",power,N);
		int compute_size = size*((int)(N/size));
        vec = (double *)malloc(compute_size * sizeof(double));
        for(i = 0 ; i < compute_size ; i++)
            vec[i] = i + 1;
    }
    MPI_Bcast(&N,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    width = (int)(N/size);
    partialto = (rank+1)*width;
    partialvec = (double*)malloc(width * sizeof(double));
    if(rank == size-1)
    {
        totalvec = (double *)malloc(size * sizeof(double));
        partial_remain = (int)N - partialto;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == size-1)
        Start = clock();
    MPI_Scatter(vec, width, MPI_DOUBLE, partialvec, width, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for(i = 0 ; i < width ; i++)
        partialsum += sqrt(partialvec[i])/N;
    MPI_Gather(&partialsum, 1, MPI_DOUBLE, totalvec, 1, MPI_DOUBLE, size-1, MPI_COMM_WORLD);
    if(rank == size-1)
    {
        for(i=0;i<size;i++)
            totalsum+=totalvec[i];
        if(partial_remain!=0)
            for(i = partial_remain;i>0;i--)
                totalsum += sqrt((double)(i+partialto))/N;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == size-1)
    {
        End = clock();
        double compute_time = ((double)(End-Start)/CLOCKS_PER_SEC)*1000;
        printf("Total sum = %lf and Time = %lf ms\n", totalsum, compute_time);
        free(totalvec);
    }
    if(rank == 0)
        free(vec);
    free(partialvec);
    MPI_Finalize();
}
