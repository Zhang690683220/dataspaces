#include <fstream>
#include <random>
#include <cstdlib>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdio>
#include "Timer.hpp"
#include "mpi.h"

const int timestep = 10;


int main(int argc, char**argv) {

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
        log.open("stdio.log", std::ofstream::out | std::ofstream::trunc);
        log << "step\twrite_ms" << std::endl;
    }

    Timer timer_write;

    for(int ts=0; ts<timestep; ts++) {

        MPI_Barrier(gcomm);
        timer_write.start();

        std::fstream file_strm;
        std::string path = "timestep_" + std::to_string(ts) + "rank_" + std::to_string(rank) + ".bin";

        if (file_strm.is_open()) {
            printf("file was left open...closing\n");
            file_strm.close();
        }

        file_strm.open(path.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);

        char* src = reinterpret_cast<char*> (array);

        file_strm.write(&src[0], bytesize);
        
        if (file_strm.fail()) {
            printf("StdFile: write failed \n");
        }

        file_strm.close();

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