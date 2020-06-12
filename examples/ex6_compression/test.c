#include "stdio.h"
#include "stdlib.h"
#include "unistd.h"
#include "stdint.h"
#include "zfp_conf.h"

int main(int argc, char** argv)
{
    int N = 8;
    int M = 8;
    int ndim = 2;


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

    zfp_type type = conf.type;      /* array scalar type */
    zfp_field *field;               /* array meta data */
    zfp_stream *zfp;                 /* compressed stream */
    void *buffer;                   /* storage for compressed stream */
    size_t bufsize;                 /* byte size of compressed buffer */
    bitstream *stream;              /* bit stream to write to or read from */
    size_t zfpsize;                 /* byte size of compressed stream */



    /* allocate meta data for the n-D array a */
    switch (conf.dims)
    {
    case 1:
        field = zfp_field_1d(data, type, ub[0]-lb[0]);
        break;
        
    case 2:
        field = zfp_field_2d(data, type, ub[0]-lb[0], ub[1]-lb[1]);
        break;
        
    case 3:
        field = zfp_field_3d(data, type, ub[0]-lb[0], ub[1]-lb[1], ub[2]-lb[2]);
        break;

    case 4:
        field = zfp_field_4d(data, type, ub[0]-lb[0], ub[1]-lb[1], ub[2]-lb[2], ub[3]-lb[3]);
        break;
        
    default:
        fprintf(stderr, "zfp only support up to 4 dimension compression\n");
        exit(1);
        break;
    }

    /* allocate meta data for a compressed stream */
    zfp = zfp_stream_open(NULL);

        
        
    /* set compression mode and parameters via one of three functions */
    if (conf.rate !=0)
    {
        zfp_stream_set_rate(zfp, conf.rate, type, conf.dims, 0);
    }
    else if(conf.precision !=0)
    {
        zfp_stream_set_precision(zfp, conf.precision);
    }
    else if(conf.tolerance !=0)
    {
        /* find the max and min of the data */
            double data_max, data_min;
            switch (type)
            {
            case zfp_type_int32:
            {
                int32_t max = *((int32_t*) data);
                int32_t min = *((int32_t*) data);
                int i; // traverse all data elements
                for(i=0;i<N*M;i++)
                {
                    if(max < *((int32_t*) (data+i*sizeof(int32_t))))
                    {
                        max = *((int32_t*) (data+i*sizeof(int32_t)));
                    }

                    if(min > *((int32_t*) (data+i*sizeof(int32_t))))
                    {
                        min = *((int32_t*) (data+i*sizeof(int32_t)));
                    }
                    data_max = (double) max;
                    data_min = (double) min;
                }
                break;
            }
            case zfp_type_int64:
            {
                int64_t max = *((int64_t*) data);
                int64_t min = *((int64_t*) data);
                int i; // traverse all data elements
                for(i=0;i<N*M;i++)
                {
                    if(max < *((int64_t*) (data+i*sizeof(int64_t))))
                    {
                        max = *((int64_t*) (data+i*sizeof(int64_t)));
                    }

                    if(min > *((int64_t*) (data+i*sizeof(int64_t))))
                    {
                        min = *((int64_t*) (data+i*sizeof(int64_t)));
                    }
                    data_max = (double) max;
                    data_min = (double) min;
                }
                break;
            }
            case zfp_type_float:
            {
                float max = *((float*) data);
                float min = *((float*) data);
                int i; // traverse all data elements
                for(i=0;i<N*M;i++)
                {
                    if(max < *((float*) (data+i*sizeof(float))))
                    {
                        max = *((float*) (data+i*sizeof(float)));
                    }

                    if(min > *((float*) (data+i*sizeof(float))))
                    {
                        min = *((float*) (data+i*sizeof(float)));
                    }
                    data_max = (double) max;
                    data_min = (double) min;
                }
                break;
            }
            case zfp_type_double:
            {
                double max = *((double*) data);
                double min = *((double*) data);
                int i; // traverse all data elements
                for(i=0;i<N*M;i++)
                {
                    if(max < *((double*) (data+i*sizeof(double))))
                    {
                        max = *((double*) (data+i*sizeof(double)));
                    }

                    if(min > *((double*) (data+i*sizeof(double))))
                    {
                        min = *((double*) (data+i*sizeof(double)));
                    }
                    data_max = (double) max;
                    data_min = (double) min;
                }
                break;
            }
            default:
                fprintf(stderr, "zfp_type error!\n");
                exit(1);
                break;
            }
            conf.max = data_max;
            conf.min = data_min;
            zfp_stream_set_accuracy(zfp, (data_max-data_min)*conf.tolerance);
        }

        /* allocate buffer for compressed data */
        bufsize = zfp_stream_maximum_size(zfp, field);
        buffer = malloc(bufsize);

        /* associate bit stream with allocated buffer */
        stream = stream_open(buffer, bufsize);
        zfp_stream_set_bit_stream(zfp, stream);
        zfp_write_header(zfp, field, ZFP_HEADER_FULL);
        zfp_stream_rewind(zfp);
        

        /* compress array and output compressed stream */
        zfpsize = zfp_compress(zfp, field);
        if (!zfpsize) {
            fprintf(stderr, "compression failed!\n");
            exit(1);
        }
        else {
            
            double *rdata = (double*) malloc(N*M*sizeof(double));
            zfp_field *rfield;

            switch (conf.dims)
                {
                case 1:
                        rfield = zfp_field_1d(rdata, type, ub[0]-lb[0]);
                        break;
        
                case 2:
                        rfield = zfp_field_2d(rdata, type, ub[0]-lb[0], ub[1]-lb[1]);
                        break;
        
                case 3:
                        rfield = zfp_field_3d(rdata, type, ub[0]-lb[0], 
                                            ub[1]-lb[1], 
                                            ub[2]-lb[2]);
                        break;

                case 4:
                        rfield = zfp_field_4d(rdata, type, ub[0]-lb[0], 
                                            ub[1]-lb[1], 
                                            ub[2]-lb[2], 
                                            ub[3]-lb[3]);
                        break;
        
                default:
                        fprintf(stderr, "zfp only support up to 4 dimension compression\n");
                        exit(1);
                        break;
                }
                
            //stream = stream_open(buffer, bufsize);
            //zfp_stream_set_bit_stream(zfp, stream);
            int header_bits = zfp_write_header(zfp, rfield, ZFP_HEADER_FULL);
            zfp_stream_rewind(zfp);

            if (!zfp_decompress(zfp, rfield)) {
                        fprintf(stderr, "decompression failed\n");
                        exit(1);
                }

                for(int i=0; i<8;i++)
                {
                        for (int j = 0; j < 8; j++)
                        {
                            printf("%lf ", *(rdata+i*8+j));
                        }
                        printf("\n");       
                }
        }
            




}