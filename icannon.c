#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

// Modifica esto para el tamaño de matriz
#define N 2000
#define TRUE 1
#define FALSE 0

void imprimematriz(int *matriz, int n);
void MatrixMultiply(int n, int *a, int *b, int *c);
int* create_array_as_matrix(int r, int c);
void populate_array_as_matrix(int *arr, int r, int c);
int array_as_matrix_equals(int *a, int *b, int r, int c);

int *a;
int *b;
int *c;

int bloque_a[4][N * N];
int bloque_b[4][N * N];
int bloque_c[4][N * N];

void MatrixMatrixMultiply(int n, int *a, int *b, int *c, int *c_grande, MPI_Comm comm)
{   
    int i;
    int nlocal;
    int npes, dims[2], periods[2];
    int myrank, my2drank, mycoords[2];
    int uprank, downrank, leftrank, rightrank, coords[2];
    int shiftsource, shiftdest;
    MPI_Status status;
    MPI_Comm comm_2d;

    /* Obtiene toda la información relacionada al canal de comunicación */
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &myrank);

    /* Crea la topología cartesiana */
    dims[0] = dims[1] = sqrt(npes);

    /* Crea los periodos para las conecciones*/
    periods[0] = periods[1] = 1;

    /* Crea la topología cartesiana, reordenando las tareas */
    MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);

    /* Obtiene las tareas y coordenadas segun la topología*/
    MPI_Comm_rank(comm_2d, &my2drank);
    MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

    /* Aqui se mueven intercambia el contenido de las tareas a la izquierda o arriba
        segun sea la matriz A o B.

        Esto es parte fundamental del algoritmo
     */
    MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
    MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);

    /* Determina la dimension de la matriz de bloque local */
    nlocal = n / dims[0];

    /* Inicializa el orden inicial de las matrices A y B*/
    MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(a, nlocal * nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);
    MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(b, nlocal * nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);

    /* Entra en el loop */
    for (i = 0; i < dims[0]; i++)
    {
        /*Aqui multiplica como lo hace cannon*/
        MatrixMultiply(nlocal, a, b, c); /*c=c+a*b*/
        /* Mueve las columnas de A a la izquierda una casilla */
        MPI_Sendrecv_replace(a, nlocal * nlocal, MPI_INT, leftrank, 1, rightrank, 1, comm_2d, &status);
        /* Mueve las filas de B arriba una casilla */
        MPI_Sendrecv_replace(b, nlocal * nlocal, MPI_INT, uprank, 1, downrank, 1, comm_2d, &status);
    }

    /* Regresa la distribución original de A y B */
    MPI_Cart_shift(comm_2d, 1, +mycoords[0], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(a, nlocal * nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);
    MPI_Cart_shift(comm_2d, 0, +mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(b, nlocal * nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);
    /* Limpia el comunicador */
    MPI_Comm_free(&comm_2d);

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

    /* Esta parte es la que va relacionando las submatrices con la matriz más grande */

    if (myrank != 0)
    {
        MPI_Reduce(c_grande, c_grande, n_grande * n_grande, MPI_INT, MPI_SUM, 0, comm);
    }
    else
    {
        MPI_Reduce(MPI_IN_PLACE, c_grande, n_grande * n_grande, MPI_INT, MPI_SUM, 0, comm);
    }
}

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

    int n_chica = N / sqrt(nro_procesos);

    int max_fila_bloque = n_chica;
    int max_columna_bloque = n_chica;

    int fila_bloque;
    int columna_bloque;
    int fila;
    int columna;

    int n_bloques = N / n_chica;
    int mi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mi_rank);
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
    /* Comienza a correr el tiempo */
    gettimeofday(&start2, 0);
    double start = MPI_Wtime();
    /* Aqui se encuentra el algoritmo de cannon com tal */
    MatrixMatrixMultiply(N, bloque_a[mi_rank], bloque_b[mi_rank], bloque_c[mi_rank], &c[0], MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (mi_rank == 0)
    {
        // Se detiene el tiempo y se muestran los resultados

        /* -----CANNON----- */


        double end = MPI_Wtime();
        // printf("\nTiempo Cannon: %.4f segundos\n", (end - start));
        gettimeofday(&stop, 0);
        fprintf(stdout, "%.8f\n",
                (stop.tv_sec + stop.tv_usec * 1e-6) - (start2.tv_sec + start2.tv_usec * 1e-6));





        // /* -----SECUENCIAL------ */
        // /* SOLO PARA COMPROBAR */

        // int *d = create_array_as_matrix(N, N);
        // /* Mide el tiempo del algoritmo secuencial */
        // gettimeofday(&start2, 0);
        // start = MPI_Wtime();
        // MatrixMultiply(N, a, b, d);
        // /* Si queremos podemos imprimirla */
        // //    imprimematriz(&d[0], N);
        // printf("Cannon \n\n");
        // //    imprimematriz(&c[0], N);

        // int equal = array_as_matrix_equals(&d[0], &c[0], N, N);
        // /* Este punto es importantisimo

        //     Se evalua si efectivamente las dos matrices (Cannon y Secuencial)
        //     estan dando el mismo resultado.


        //  */
        // if (equal)
        // {
        //     /* Si son iguales, quiere decir que cannon esta bien ejecutado 
            
        //         y Muestra el tiempo del Algoritmo secuencial
            
        //      */
        //     printf("\nSon iguales\n");
        //     end = MPI_Wtime();
        //     printf("\nTiempo Secuencial: %.4f segundos\n", (end - start));
        //     gettimeofday(&stop, 0);
        //     fprintf(stdout, "Time = %.6f\n\n",
        //             (stop.tv_sec + stop.tv_usec * 1e-6) - (start2.tv_sec + start2.tv_usec * 1e-6));
        // }
        // else
        // {
        //     printf("\n No son iguales\n");
        // }
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