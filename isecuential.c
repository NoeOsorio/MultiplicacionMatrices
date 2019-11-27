#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

// Modifica esto para el tamaño de matriz
#define N 1024
#define TRUE 1
#define FALSE 0

void imprimematriz(int *matriz, int n);
void MatrixMultiply(int n, int *a, int *b, int *c);
int *create_array_as_matrix(int r, int c);
void populate_array_as_matrix(int *arr, int r, int c);
int array_as_matrix_equals(int *a, int *b, int r, int c);

int *a;
int *b;
int *c;

void MatrixMultiply(int n, int *a, int *b, int *c)
{
    int i, j, k;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++)
                c[i * n + j] += a[i * n + k] * b[k * n + j];
}

void imprimematriz(int *matriz, int n)
{
    int i;
    int j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%d\t", matriz[i * n + j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    // Inicializa las matrices
    struct timeval start2, stop;
    a = create_array_as_matrix(N, N);
    b = create_array_as_matrix(N, N);
    c = create_array_as_matrix(N, N);

    populate_array_as_matrix(&a[0], N, N);
    populate_array_as_matrix(&b[0], N, N);

    /* Si queremos imprimir las matrices A y B */
    // printf("A \n");
    // imprimematriz(&a[0], N);
    // printf("B \n");
    // imprimematriz(&b[0], N);

    MPI_Init(&argc, &argv);
    int nro_procesos;
    MPI_Comm_size(MPI_COMM_WORLD, &nro_procesos);

    int mi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mi_rank);

    if (mi_rank == 0)
    {

        /* -----SECUENCIAL------ */

        /* Mide el tiempo del algoritmo secuencial */
        gettimeofday(&start2, 0);
        MatrixMultiply(N, a, b, c);
        /* Si queremos podemos imprimirla */
        // imprimematriz(&c[0], N);
        gettimeofday(&stop, 0);
        // Se detiene el tiempo y se muestran los resultados
        fprintf(stdout, "%0.10f\n",
                (stop.tv_sec + stop.tv_usec * 1e-6) - (start2.tv_sec + start2.tv_usec * 1e-6));
    }

    MPI_Finalize();
    return 0;
}

int *create_array_as_matrix(int r, int c)
{
    int *mat = calloc(r * c, sizeof(int));
    return mat;
}

void populate_array_as_matrix(int *arr, int r, int c)
{
    int j;
    for (j = 0; j < r * c; j++)
    {
        arr[j] = rand() % 2000 + 1000;
    }
}

int array_as_matrix_equals(int *a, int *b, int r, int c)
{
    int i = 0;
    for (i = 0; i < r * c; i++)
    {
        if (a[i] != b[i])
        {
            return FALSE;
        }
    }
    return TRUE;
}