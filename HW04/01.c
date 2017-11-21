#include <stdio.h>
#include <stdlib.h>
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
    int partialfrom = rank*width+1;
    int partial_remain = 0;
    int partialsum=0;
    int totalsum=0;
    int i;
    int *vec = NULL;
    int *partialvec = (int*)malloc(width*sizeof(int));
    int *totalvec = NULL;
    if(rank == size -1)
    {
        vec = (int *)malloc(partialto*sizeof(int));
        totalvec = (int *)malloc(size*sizeof(int));
        partial_remain = FINAL-partialto;
        for(i=0;i<partialto;i++)
            vec[i]=i+1;
    }
    MPI_Scatter(vec, width, MPI_INT, partialvec, width, MPI_INT, size-1, MPI_COMM_WORLD);
    for(i=0;i<width;i++)
        partialsum+=partialvec[i];
    MPI_Gather(&partialsum, 1, MPI_INT, totalvec, 1, MPI_INT, size-1, MPI_COMM_WORLD);
    if(rank==size-1)
    {
        for(i = 0 ; i < size ; i++)
            totalsum += totalvec[i];
        if(partial_remain!=0)
            for(i = partial_remain ; i>0;i--)
                totalsum += i+partialto;
        printf("total sum = %d\n", totalsum);
        free(vec);
        free(totalvec);
    }
    free(partialvec);
    MPI_Finalize();
}
