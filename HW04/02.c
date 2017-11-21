#include <stdio.h>
#include <mpi.h>
#define FINAL 10000
int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int width = (int)(FINAL/size);
    int partialto = (rank+1)*width;
    int partialfrom = rank*width;
    int partialsum=0;
    int totalsum;
    int i;
    if(rank == size-1)
        partialto = FINAL+1;
    for(i = partialfrom ; i < partialto ; i++)
    {
        partialsum += i;
    }
    MPI_Reduce(&partialsum, &totalsum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(rank == 0)
    {
        printf("total sum = %d\n", totalsum);
    }
    MPI_Finalize();
}
