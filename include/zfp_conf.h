
#ifndef __ZFP_CONF_H
#define __ZFP_CONF_H
#include "zfp.h"

typedef struct 
{
    zfp_type        type;
    double          rate;
    uint            precision;
    double          tolerance;
    int             dims;
    double          min, max;
} zfp_conf;

#endif
