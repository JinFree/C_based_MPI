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
    int i, j, N, N_fix, cols_p, rows_p, sum = 0;
    double tic, toc;
    if(rank == 0)
    {
        fflush(stdin);
        printf("It will compute A * x = y\n");
        printf("A will be N by N Matrix, x and y will be N by 1 Matrix.\n");
        printf("Input N\n");
        scanf("%d", &N);
        N_fix = N;
        while(N_fix%size!=0)
            N_fix++;
    }
    MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&N_fix,1,MPI_INT,0,MPI_COMM_WORLD);
    cols_p = N;
    rows_p = N_fix/size;
    if(rank == size-1)
        rows_p = N-rows_p*(size-1);
    int *A_p = (int *)calloc(rows_p * cols_p, sizeof(int));
    int *x_p = (int *)calloc(rows_p * 1, sizeof(int));
    int *temp_p = (int *)calloc(1 * cols_p, sizeof(int));
    int *y_p = (int *)calloc(rows_p * 1, sizeof(int));
    int *recvcounts = malloc(sizeof(int)*size);
    int *displs = malloc(sizeof(int)*size);
    // Define Matrix
    srand(time(NULL));
    for( j = 0 ; j < rows_p ; j++)
    {
        for( i = 0 ; i < cols_p ; i++)
        {
            int coord = j * cols_p + i;
            A_p[coord] = (rand()/(rank+1))%10;
        }
        x_p[j] = (rand()/(rank+1))%10;
    }
    for( i = 0 ; i < size-1 ; i++ )
    {
        recvcounts[i] = N_fix/size;
        displs[i] = sum;
        sum += recvcounts[i];
    }
    recvcounts[size-1] = cols_p-sum;
    displs[size-1] = sum;
    printf("rank = %d, N = %d, N_fix = %d, rows_p = %d\n", rank, N, N_fix, rows_p);
    if( rank == 0 || rank == size-1)
    {
        for( i = 0 ; i < size ; i++ )
        {
            printf("rank = %d, recvcounts[%d] = %d, displs[%d] = %d\n", rank, i, recvcounts[i], i, displs[i]);
        }
    }
    int rank_check;
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
        tic = MPI_Wtime();
    for(rank_check = 0 ; rank_check < size ; rank_check++)
    { 
        MPI_Gatherv(x_p, recvcounts[rank], MPI_INT, temp_p, recvcounts, displs, MPI_INT, rank_check, MPI_COMM_WORLD);
    }
    for( j = 0 ; j < rows_p ; j++ )
    {
        for( i = 0 ; i < cols_p ; i++)
        {
            int coord = j * cols_p + i;
            y_p[j] += A_p[coord] * temp_p[i];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
    {
        toc = MPI_Wtime();
        printf("Total Excute Time = %lf seconds\n", toc-tic);
    }
    // Print A, x, y
    printf("rank = %d, ", rank);
    for( j = 0 ; j < rows_p ; j++)
    {
        for( i = 0 ; i < cols_p ; i++)
        {
            int coord = j * cols_p + i;
            printf("A[%d,%d] = %d, ", j,i, A_p[coord]);
        }
        printf("x[%d] = %d, ", j, x_p[j]);
        printf("y[%d] = %d, ", j, y_p[j]);
    }
    printf("\n");
    free(A_p);
    free(x_p);
    free(y_p);
    MPI_Finalize();
}
