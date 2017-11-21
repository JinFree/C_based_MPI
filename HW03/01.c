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
    int temp;
    if(rank == size-1)
        partialto = FINAL+1;
    for(int i = partialfrom ; i < partialto ; i++)
    {
        partialsum += i;
    }
    for(int i=1;i<size;i++)
    {
        if(rank == i)
            MPI_Send(&partialsum,1, MPI_INT, 0, i, MPI_COMM_WORLD);
        else if(rank == 0)
        {
            MPI_Recv(&temp, 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            partialsum += temp;
        }
    }
    if(rank == 0)
    {
        printf("total sum = %d\n", partialsum);
    }
    MPI_Finalize();
}
