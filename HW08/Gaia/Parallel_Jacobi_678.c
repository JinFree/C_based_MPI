#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#define PI 3.14159265359
#define Tol 0.0001
void Elliptic(int N, int rank, int size);
void FileWriter(double *U, int N, double delta, double time, int size);
double L2Error(double *U, int N, double delta);
void ExactWriter(int N, double delta);
int main(int argc, char **argv)
{
    int rank, size, N;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	/*if(rank == 0)
    {
        fflush(stdin);
        printf("It will Elliptic, Jacobi Method \n");
        printf("Input N\n");
        scanf("%d", &N);
    }
	MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);*/
    Elliptic(700, rank, size);
    Elliptic(650, rank, size);
    Elliptic(550, rank, size);
	if(rank < 7)
	{
        Elliptic(700, rank, 7);
		Elliptic(650, rank, 7);
		Elliptic(600, rank, 7);
		Elliptic(550, rank, 7);
    }
	if(rank < 6)
	{
        Elliptic(700, rank, 6);
		Elliptic(650, rank, 6);
		Elliptic(600, rank, 6);
		Elliptic(550, rank, 6);
	}
    MPI_Finalize();
}
void Elliptic(int N, int rank, int size)
{
    int NN = N*N, N_fix = N*N;
    int rows_p, index, i, j, pos, iter=0;
    double tic, toc, x, y, f, delta = 2.0/((double)(N-1));
    double Error=0.0, Error_p=0.0;
    while(N_fix%size!=0)
        N_fix++;
    rows_p = (int)(N_fix/size);
    if(rank == size-1)
        rows_p = NN-rows_p*(size-1);
    double *U = (double *)malloc(NN*sizeof(double));
    double *U_p = (double *)calloc(rows_p, sizeof(double));
    int *counts = malloc(sizeof(int)*size);
    int *displs = malloc(sizeof(int)*size);
    int sum = 0;
    for( i = 0 ; i < size-1 ; i++ )
    {
        counts[i] = (int)(N_fix/size);
        displs[i] = sum;
        sum += counts[i];
    }
    counts[size-1] = NN-sum;
    displs[size-1] = sum;
    for( index = 0 ; index < counts[rank] ; index++ )
    {
        pos = index + displs[rank];
        i = (int)(pos % N);
        j = (int)(pos / N);
        if( pos != j * N + i )
        {
            printf("Wrong index\n");
        }
        if( i != 0 && j != 0 && i != N-1 && j != N-1 )
            U_p[index] = 0.0;
        else if( i == 0 )   // 0,y
            U_p[index] = sin(PI+delta*j);
        else if( i == N-1 ) // 2,y
            U_p[index] = sin(PI+delta*j);  
        else if( j == 0 )   // x,0
            U_p[index] = 0.0;
        else if( j == N-1 ) // x,2
            U_p[index] = cos(PI*delta*i)*sin(PI+2.0);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if( rank == 0 )
        tic = MPI_Wtime();
    Error = 1.0;
    while(Error > Tol)
    {
        iter++;
        MPI_Allgatherv(U_p, counts[rank], MPI_DOUBLE, U, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        for( index = 0 ; index < counts[rank] ; index++ )
        {
            pos = index + displs[rank];
            i = (int)(pos % N);
            j = (int)(pos / N);
            x = delta * i;
            y = delta * j;
            f = cos(PI*x)*sin(y)+PI*PI*cos(PI*x)*sin(y);
            if( i != 0 && j != 0 && i != N-1 && j != N-1 )
                U_p[index] = (U[pos+1]+U[pos-1]+U[pos+N]+U[pos-N]-delta*delta*f)*0.25;
            else if( i == 0 )   // 0,y
                U_p[index] = sin(PI+delta*j);
            else if( i == N-1 ) // 2,y
                U_p[index] = sin(PI+delta*j);  
            else if( j == 0 )   // x,0
                U_p[index] = 0.0;
            else if( j == N-1 ) // x,2
                U_p[index] = cos(PI*delta*i)*sin(PI+2.0);
        }
        Error_p = 0.0;
        for( index = 0 ; index < counts[rank] ; index++ )
        {
            pos = index + displs[rank];
            Error_p += fabs(U_p[index] - U[pos]);
        }
        MPI_Allreduce(&Error_p,&Error,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if( rank == 0 )
			printf("\rL1=%.5lf",Error);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if( rank == 0 )
    {
        toc = MPI_Wtime();
        printf("\n");
    }
    MPI_Gatherv(U_p, counts[rank], MPI_DOUBLE, U, counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if( rank == 0 )
    {
        FileWriter(U, N, delta, toc-tic, size);
        ExactWriter(N, delta);
    }
    free(U);
    free(U_p);
    free(counts);
    free(displs);
}
void FileWriter(double *U, int N, double delta, double time, int size)
{
    double L2 = L2Error(U, N, delta);
    FILE *Solve;
    char Solname[100];
    int i,j,pos;
    double x,y;
    sprintf(Solname, "Elliptic, N=%d, L2=%.13lf, np=%d, t=%lf, Jacobi.csv",N,L2, size, time);
    Solve = fopen(Solname, "w");
    fprintf(Solve, "X,Y,U\n");
    for( j = 0 ; j < N ; j++ )
    {
        for( i = 0 ; i < N ; i++ )
        {
            pos = j * N + i;
            x = delta*i;
            y = delta*j;
            fprintf(Solve, "%lf,%lf,%lf\n", x, y, U[pos]);
        }
    }
    fclose(Solve);
    return;
}
double L2Error(double *U, int N, double delta)
{
    double error = 0.0, x, y;
    int i,j;
    for( j = 0 ; j < N ; j++ )
    {
        for( i = 0 ; i < N ; i++ )
        {
            int pos = j*N+i;
            x = delta*i;
            y = delta*j;
            double abs_error = fabs(U[pos]-cos(PI*x)*sin(PI+y));
            error += abs_error * abs_error;
        }
    }
    error = sqrt(error)/(double)(N*N);
    return error;
}
void ExactWriter(int N, double delta)
{
    FILE *fp;
    char name[50];
    int i, j;
    double x, y, Exact;
    sprintf(name, "Elliptic, N=%d, Exact.csv", N);
    fp = fopen(name, "w");
    fprintf(fp, "X,Y,U\n");
    for( j = 0 ; j < N ; j++ )
    {
        for( i = 0 ; i < N ; i++ )
        {
            x = delta*i;
            y = delta*j;
            Exact = cos(PI*x)*sin(PI+y);
            fprintf(fp, "%lf,%lf,%lf\n", x, y, Exact);
        }
    }
    fclose(fp);
}