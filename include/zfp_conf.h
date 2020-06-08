
#ifndef __ZFP_CONF_H
#define __ZFP_CONF_H
#include "zfp.h"
#include "bbox.h"

typedef struct 
{
    zfp_type        type;
    double          rate;
    uint            precision;
    double          tolerance;
    int             dims;
    double          min, max;
    struct bbox     parent_bbox;       
} zfp_conf;

#endif
