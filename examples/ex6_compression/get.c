#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "dataspaces.h"
#include "mpi.h"

int main(int argc, char** argv)
{
    int rank, nprocs; 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm gcomm = MPI_COMM_WORLD;

    int N = 8;
    int M = 8;


    double *data = (double*) malloc(N*M*sizeof(double));
    
  
    dspaces_init(nprocs, 2, &gcomm, NULL);

    for(int timestep=0; timestep<10, timestep++)
    {
        dspaces_lock_on_read("my_test_lock", NULL);

        char var_name[128];
        sprintf(var_name, "ex6_sample_data");

        int ndim = 2;

        uint64_t lb[2] = {0}, ub[2] = {0};

        ub[0] = N-1; 
        ub[1] = M-1;


        dspaces_get(var_name, timestep, sizeof(double), ndim, lb, ub, data);

        dspaces_unlock_on_read("my_test_lock", NULL);
    }

    for(int i=0; i<N;i++)
    {
        for (int j = 0; j < M; j++)
        {
            printf("%lf ",*(data+i*N+j));
        }
        printf("\n");       
    }
}