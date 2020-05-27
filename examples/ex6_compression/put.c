#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "zfp_conf.h"
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

    for(int i=0; i<N;i++)
    {
        for (int j = 0; j < M; j++)
        {
            *(data+i*N+j) = 1.0*i*N+j;
            printf("%lf ",*(data+i*N+j));
        }
        printf("\n");       
    }

    dspaces_init(nprocs, 1, &gcomm, NULL);

    for(int timestep=0; timestep<10, timestep++)
    {
        dspaces_lock_on_write("my_test_lock", NULL);

        char var_name[128];
        sprintf(var_name, "ex6_sample_data");

        int ndim = 2;

        uint64_t lb[2] = {0}, ub[2] = {0};

        ub[0] = N-1; 
        ub[1] = M-1;

        zfp_conf conf = {
            .type = zfp_type_double,
            .rate = 0,
            .precision = 0,
            .tolerance = 1e-1,
            .dims = ndim
        };

        dspaces_put(var_name, timestep, sizeof(double), ndim, lb, ub, data, 1, &conf);

        dspaces_unlock_on_write("my_test_lock", NULL);
    }

    printf("finished!\n");

    return 0;


}