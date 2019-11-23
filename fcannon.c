#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "util/matrixFloat.h"

#define N 4
void imprimematriz(double *matriz, int n);
void MatrixMultiply(int n, double *a, double *b, double *c);
double *a;
double *b;
double *c;

// double bloque_a[25][N * N];
// double bloque_b[25][N * N];
// double bloque_c[25][N * N];

void MatrixMatrixMultiply(int n, double *a, double *b, double *c, double *c_grande, MPI_Comm comm)
{
    int i;
    int nlocal;
    int npes, dims[2], periods[2];
    int myrank, my2drank, mycoords[2];
    int uprank, downrank, leftrank, rightrank, coords[2];
    int shiftsource, shiftdest;
    MPI_Status status;
    MPI_Comm comm_2d;

    /* Get the communicator related information */
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &myrank);

    /* Set up the Cartesian topology */
    dims[0] = dims[1] = sqrt(npes);

    /* Set the periods for wraparound connections */
    periods[0] = periods[1] = 1;

    /* Create the Cartesian topology, with rank reordering */
    MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);

    /* Get the rank and coordinates with respect to the new topology */
    MPI_Comm_rank(comm_2d, &my2drank);
    MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

    /* Compute ranks of the up and left shifts */
    MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
    MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);

    /* Determine the dimension of the local matrix block */
    nlocal = n / dims[0];

    /* Perform the initial matrix alignment. First for A and then for B */
    MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(a, nlocal * nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);
    MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(b, nlocal * nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);

    /* Get into the main computation loop */
    for (i = 0; i < dims[0]; i++)
    {
        MatrixMultiply(nlocal, a, b, c); /*c=c+a*b*/
        /* Shift matrix a left by one */
        MPI_Sendrecv_replace(a, nlocal * nlocal, MPI_INT, leftrank, 1, rightrank, 1, comm_2d, &status);
        /* Shift matrix b up by one */
        MPI_Sendrecv_replace(b, nlocal * nlocal, MPI_INT, uprank, 1, downrank, 1, comm_2d, &status);
    }

    /* Restore the original distribution of a and b */
    MPI_Cart_shift(comm_2d, 1, +mycoords[0], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(a, nlocal * nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);
    MPI_Cart_shift(comm_2d, 0, +mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(b, nlocal * nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);
    MPI_Comm_free(&comm_2d); /* Free up communicator */

    int ind;
    int n_grande = n;
    int ind_col;
    for (ind = 0; ind < nlocal; ind++)
    {
        for (ind_col = 0; ind_col < nlocal; ind_col++)
        {
            int fila_grande = mycoords[0] * nlocal + ind;
            int columna_grande = mycoords[1] * nlocal + ind_col;
            c_grande[fila_grande * n_grande + columna_grande] = c[ind * nlocal + ind_col];
        }
    }

    if (myrank != 0)
    {
        MPI_Reduce(c_grande, c_grande, n_grande * n_grande, MPI_INT, MPI_SUM, 0, comm);
    }
    else
    {
        MPI_Reduce(MPI_IN_PLACE, c_grande, n_grande * n_grande, MPI_INT, MPI_SUM, 0, comm);
    }
}

/* This function performs a serial matrix-matrix multiplication c = a*b */
void MatrixMultiply(int n, double *a, double *b, double *c)
{
    int i, j, k;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++)
                c[i * n + j] += a[i * n + k] * b[k * n + j];
}

void imprimematriz(double *matriz, int n)
{
    int i;
    int j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%f\t", matriz[i * n + j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    // int *a;
    // int *b;
    // int *c;

    // // a = &A[0][0];
    // printf("Matrix B %d \n", &b);
    // // b = &B[0][0];
    // printf("Matrix C %d \n", &c);
    // c = &C[0][0];

    // printf("A \n");
    // imprimematriz(&a[0], N);
    // printf("B \n");
    // imprimematriz(&b[0], N);

    a = create_array_as_matrix(N, N);
    b = create_array_as_matrix(N, N);
    c = create_array_as_matrix(N, N);
    // printf("Matrix A, %p \n", &a[0]);
    // printf("Matrix B, %p \n", &b[0]);
    // printf("Matrix C, %p \n", &c[0]);
    populate_array_as_matrix(&a[0], N, N);
    populate_array_as_matrix(&b[0], N, N);

    printf("A \n");
    imprimematriz(&a[0], N);
    printf("B \n");
    imprimematriz(&b[0], N);

    // printf("Matriz creadas \n");

    MPI_Init(&argc, &argv);
    int nro_procesos;
    MPI_Comm_size(MPI_COMM_WORLD, &nro_procesos);

    printf("MPI inicializado %d\n", nro_procesos);
    int n_chica = N / sqrt(nro_procesos);
    printf("N chica %d\n", n_chica);
    printf("mini matrices de %d x %d", nro_procesos, n_chica * n_chica);

    int max_fila_bloque = n_chica;
    int max_columna_bloque = n_chica;

    int fila_bloque;
    int columna_bloque;
    int fila;
    int columna;

    int n_bloques = N / n_chica;

    printf("Creando subbloques \n");
    double bloque_a[nro_procesos][n_chica * n_chica];
    double bloque_b[nro_procesos][n_chica * n_chica];
    double bloque_c[nro_procesos][n_chica * n_chica];

    int mi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mi_rank);
    printf("Ranks listos \n");
    printf("Rank %d", mi_rank);
    columna_bloque = mi_rank % n_bloques;
    fila_bloque = (mi_rank - columna_bloque) / n_bloques;

    int indice_bloque = 0;

    for (fila = fila_bloque * n_chica; fila < fila_bloque * n_chica + n_chica; fila++)
    {
        for (columna = columna_bloque * n_chica; columna < columna_bloque * n_chica + n_chica; columna++)
        {
            bloque_a[mi_rank][indice_bloque] = a[fila * N + columna];
            bloque_b[mi_rank][indice_bloque] = b[fila * N + columna];
            bloque_c[mi_rank][indice_bloque] = 0;
            indice_bloque++;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    MatrixMatrixMultiply(N, bloque_a[mi_rank], bloque_b[mi_rank], bloque_c[mi_rank], &c[0], MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();

    if (mi_rank == 0)
    {
        // Serial multiplication. In order to check te results.
        double *d = create_array_as_matrix(N, N);
        MatrixMultiply(N, a, b, d);
           imprimematriz(&d[0], N);
        printf("Cannon \n\n");
           imprimematriz(&c[0], N);

        int equal = array_as_matrix_equals(&d[0], &c[0], N, N);
        if (equal)
        {
            printf("\nSon iguales\n");
        }
        else{
            printf("\n No son iguales\n");
        }

        printf("\nTiempo: %.4f segundos\n", (end - start));
    }

    MPI_Finalize();
    return 0;
}

void create_matrix(Matrix *m, int nrow, int ncol)
{
    int i;

    m->nrow = nrow;
    m->ncol = ncol;
    m->data = malloc(nrow * sizeof(int *));
    for (i = 0; i < ncol; i++)
    {
        m->data[i] = calloc(ncol, sizeof(int));
    }
}

void populate_matrix(Matrix *m)
{
    int i, j;
    for (i = 0; i < m->nrow; i++)
    {
        for (j = 0; j < m->ncol; j++)
        {
            m->data[i][j] = rand() % 10 + 1;
        }
    }
}

void print_matrix(Matrix *m, char iden)
{
    int i, j;

    for (i = 0; i < m->nrow; i++)
    {
        for (j = 0; j < m->ncol; j++)
        {
            printf("%c[%d][%d] = %d  ", iden, i, j, m->data[i][j]);
        }
        printf("\n");
    }
}

/**
 * shift the given matrix left.
 *
 * @param m: the matrix to shift.
 * @param initial: a value > 0 indicates that it is a first shift, otherwise is
 *                 normal shift.
 */
void shift_matrix_left(Matrix *m, int block_sz, int initial)
{
    int i, j, k, s, step = block_sz;
    Matrix aux;

    create_matrix(&aux, 1, m->ncol);
    for (k = 0, s = 0; k < m->ncol; k += block_sz, s++)
    {
        for (i = k; i < (k + block_sz); i++)
        {
            if (initial > 0)
            {
                step = s * block_sz;
            }
            for (j = 0; j < m->ncol; j++)
            {
                aux.data[0][j] = m->data[i][(j + step) % m->ncol];
            }
            for (j = 0; j < m->ncol; j++)
            {
                m->data[i][j] = aux.data[0][j];
            }
        }
    }
}

/**
 * shift the given matrix up.
 *
 * @param m: the matrix to shift.
 * @param initial: a value > 0 indicates that it is a first shift, otherwise is
 *                 normal shift.
 */
void shift_matrix_up(Matrix *m, int block_sz, int initial)
{
    int i, j, k, s, step = block_sz;
    Matrix aux;

    create_matrix(&aux, 1, m->nrow);
    for (k = 0, s = 0; k < m->nrow; k += block_sz, s++)
    {
        for (i = k; i < (k + block_sz); i++)
        {
            if (initial > 0)
            {
                step = s * block_sz;
            }
            for (j = 0; j < m->nrow; j++)
            {
                aux.data[0][j] = m->data[(j + step) % m->nrow][i];
            }
            for (j = 0; j < m->nrow; j++)
            {
                m->data[j][i] = aux.data[0][j];
            }
        }
    }
}

/**
 * Matrix multiplication
 */
void matrix_product(Matrix *c, Matrix *a, Matrix *b)
{
    int r, s, k;

    for (r = 0; r < a->nrow; r++)
    {
        for (s = 0; s < b->ncol; s++)
        {
            for (k = 0; k < a->ncol; k++)
            {
                c->data[r][s] += a->data[r][k] * b->data[k][s];
            }
        }
    }
}

double *create_array_as_matrix(int r, int c)
{
    // printf("Creando matriz\n");
    //  double *mat = (int *)malloc(r * c * sizeof(*mat));
    double *mat = calloc(r * c, sizeof(*mat));
    // printf("YA se creo: %f \n", *mat);
    return mat;
}

void populate_array_as_matrix(double *arr, int r, int c)
{
    // printf("Cargando matriz de %p\n", &arr[0]);
    int j;
    for (j = 0; j < r * c; j++)
    {
        // arr[j] = rand() % 10 + 1;
        arr[j] = j * 10;
        // printf("%d ", arr[j]);
        // printf("%d ", j);
    }
    // printf("Matriz %p con %f \n ", &arr[0], arr[r * c - 1]);
}

int array_as_matrix_equals(double *a, double *b, int r, int c)
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