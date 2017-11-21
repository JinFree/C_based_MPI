#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int i, width, partialto;
    double N, partialsum = 0, totalsum = 0;
    time_t Start, End;
    if(rank == 0)
    {
        int power;
        fflush(stdin);
        printf("It Will Compute N = 10^(power), Input power\n");
        scanf("%d",&power);
        N = pow(10.0,(double)power);
        printf("N = 10^%d = %.1lf\n",power,N);
    }
    MPI_Bcast(&N,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    width = (int)(N/size);
    partialto = (rank+1)*width;
    if(rank == size-1)
        partialto = (int)N;
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
        Start = clock();
    for(i = rank*width+1 ; i <= partialto ; i++)
        partialsum += sqrt((double)i)/N;
    MPI_Reduce(&partialsum,&totalsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank ==0)
    {
        End = clock();
        double compute_time = ((double)(End-Start)/CLOCKS_PER_SEC)*1000;
        printf("Total sum = %lf and Time = %lf ms\n", totalsum, compute_time);
    }
    MPI_Finalize();
}