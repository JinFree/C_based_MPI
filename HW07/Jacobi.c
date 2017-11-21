#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265359
#define Tol 0.0001
void Elliptic(int N);
void Bound_Cond(double *U, double *Unew, int N, double delta);
double CallJacobi(double *U, double *Unew, int N, double delta);
void Jacobi(double *U, double *Unew, int N, double delta);
double LinfError(double *U, double *Unew, int N);
void FileWriter(double *U, int N, double delta);
double L2Error(double *U, int N, double delta);
void ExactWriter(int N, double delta);
int main(int argc, char **argv)
{
    Elliptic(50);
    Elliptic(100);
    Elliptic(200);
    return;
}
void Elliptic(int N)
{
    double Linf = 0, delta = 2.0/((double)(N-1));
    double *U = (double *)calloc(sizeof(double),N*N);
    double *Unew = (double *)calloc(sizeof(double),N*N);
    Bound_Cond(U, Unew, N, delta);
    Linf = CallJacobi(U, Unew, N, delta);
    FileWriter(U, N, delta);
    ExactWriter(N, delta);
    free(U);
    free(Unew);
}
void Bound_Cond(double *U, double *Unew, int N, double delta)
{
    int i,j;
    for( i = 0 ; i < N ; i++ )
    {
       U[i] = Unew[i] = 0.0;                                         // x,0
       U[(N-1)*N+i] = Unew[(N-1)*N+i] = cos(PI*delta*i)*sin(PI+2.0); // x,2
       U[i*N] = Unew[i*N] = sin(PI+delta*i);                         // 0,y
       U[i*N+N-1] = Unew[i*N+N-1] = sin(PI+delta*i);                 // 2,y
    }
}
double CallJacobi(double *U, double *Unew, int N, double delta)
{
    double Error = 1.0;
    int i,j, iter=0;
    while(Error > Tol)
    {
        iter++;
        Jacobi(U, Unew, N, delta);
        Error = LinfError(U, Unew, N);
        for( j = 1 ; j < N-1 ; j++ )
            for( i = 1 ; i < N-1 ; i++ )
                U[j*N+i] = Unew[j*N+i];
        printf("\rN=%d, iter=%d, Linf=%.13lf",N,iter,Error);
    };
    printf("\n");
    return Error;
}
void Jacobi(double *U, double *Unew, int N, double delta)
{
    int i,j, pos;
    double x, y, f;
    for( j = 1 ; j < N-1 ; j++ )
    {
        for( i = 1 ; i < N-1 ; i++ )
        {
            pos = j * N + i;
            x = delta*i;
            y = delta*j;
            f = cos(PI*x)*sin(y)+PI*PI*cos(PI*x)*sin(y);
            Unew[pos] = (U[pos+1]+U[pos-1]+U[pos+N]+U[pos-N]-delta*delta*f)*0.25;
        }
    }
}
double LinfError(double *U, double *Unew, int N)
{
    double error = 0.0;
    int i,j;
    for( j = 1 ; j < N-1 ; j++ )
    {
        for( i = 1 ; i < N-1 ; i++ )
        {
            int pos = j*N+i;
            double abs_error = fabs(Unew[pos]-U[pos]);
            error = (error < abs_error) ? abs_error : error;
        }
    }
    return error;
}
void FileWriter(double *U, int N, double delta)
{
    double L2 = L2Error(U, N, delta);
    FILE *fp;
    char name[50];
    int i, j, pos;
    double x, y;
    sprintf(name, "Elliptic, N=%d, L2=%lf, Jacobi.csv", N, L2);
    fp = fopen(name, "w");
    fprintf(fp, "X,Y,U\n");
    for( j = 0 ; j < N ; j++ )
    {
        for( i = 0 ; i < N ; i++ )
        {
            pos = j * N + i;
            x = delta*i;
            y = delta*j;
            fprintf(fp, "%lf,%lf,%lf\n", x, y, U[pos]);
        }
    }
    fclose(fp);
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
    error = sqrt(error)/(double)(N-1);
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
