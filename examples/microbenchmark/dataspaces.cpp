#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <random>
#include "dataspaces.h"
#include "Timer.hpp"
#include "mpi.h"

const int timestep = 10;

int main(int argc, char** argv) {
    int rank, nprocs; 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm gcomm = MPI_COMM_WORLD;

    int dim0 = 128;
    int dim1 = 128;
    int dim2 = 256;
    int bytesize = dim0*dim1*dim2*sizeof(double);

    double *array = (double*) malloc(bytesize);

    double lower_bound = 0;
    double upper_bound = 10000;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    for(int i=0; i<dim0*dim1*dim2; i++) {
        array[i] = unif(re);
    }

    std::ofstream log;

    if(rank == 0) {
        log.open("dataspaces.log", std::ofstream::out | std::ofstream::trunc);
        log << "step\twrite_ms" << std::endl;
    }

    Timer timer_write;

    dspaces_init(nprocs, 1, &gcomm, NULL);

    for(int ts=0; ts<timestep; ts++) {

        MPI_Barrier(gcomm);
        timer_write.start();

        uint64_t lb[3];
        uint64_t ub[3];
        lb[0] = (rank/8)*128;
        lb[1] = ((rank%8)/2)*128;
        lb[2] = (rank%2)*256;
        ub[0] = lb[0]+127;
        ub[1] = lb[1]+127;
        ub[2] = lb[2]+255;
        void* src = reinterpret_cast<void*> (array);
        dspaces_lock_on_write("my_test_lock", &gcomm);
        int err = dspaces_put("random_data", ts, sizeof(double), 3, lb, ub, src);
        dspaces_unlock_on_write("my_test_lock", &gcomm);

        if(err) {
            printf("Dataspaces: write failed \n");
        }

        double time_write = timer_write.stop();
        MPI_Barrier(gcomm);
        double *avg_time = NULL;

        if(rank == 0) {
            avg_time = (double*) malloc(sizeof(double)*nprocs);
        }

        MPI_Gather(&time_write, 1, MPI_DOUBLE, avg_time, 1, MPI_DOUBLE, 0, gcomm);

        if(rank == 0) {
            double avg = 0;
            for(int i=0; i<nprocs; i++) {
                avg += avg_time[i];
            }
            avg = avg/nprocs;
            log << ts << "\t" << avg << std::endl;
            free(avg_time);
        }

    }

    double totaltime = timer_write.elapsed();

    double *avg_totaltime = NULL;

    if(rank == 0) {
        avg_totaltime = (double*) malloc(sizeof(double)*nprocs);
    }

    MPI_Gather(&totaltime, 1, MPI_DOUBLE, avg_totaltime, 1, MPI_DOUBLE, 0, gcomm);

    if(rank == 0) {
        double avg = 0;
        for(int i=0; i<nprocs; i++) {
            avg += avg_totaltime[i];
        }
        avg = avg/nprocs;
        log << "total\t" << avg << std::endl;
        free(avg_totaltime);
        log.close();
    }

    MPI_Finalize();

    return 0;

}