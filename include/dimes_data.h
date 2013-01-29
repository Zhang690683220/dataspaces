/*
 * Copyright (c) 2009, NSF Cloud and Autonomic Computing Center, Rutgers University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided
 * that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this list of conditions and
 * the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided with the distribution.
 * - Neither the name of the NSF Cloud and Autonomic Computing Center, Rutgers University, nor the names of its
 * contributors may be used to endorse or promote products derived from this software without specific prior
 * written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*
 * Fan Zhang (2012) TASSL Rutgers University
 * zhangfan@cac.rutgers.edu
 */
#ifndef __DIMES_DATA_H__
#define __DIMES_DATA_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "dart.h"
#include "ss_data.h"

struct ptlid_map;

struct box_2_pointers{
        void *ptr1;
        void *ptr2;
};

/****************************************************************************/
// Header structure for dimes_put request.
struct hdr_dimes_put {
        int has_rdma_data;
        struct obj_descriptor odsc;
} __attribute__((__packed__));

// Header structure for dimes_obj_get request.
struct hdr_dimes_obj_get {
        int qid;
        int rank;
        int rc;
        union {
                struct {
                        /* Number of directory entries. */
                        int                     num_de;
                        struct obj_descriptor   odsc;
                } o;
                struct {
                        /* Number of versions available. */
                        int                     num_vers;
                        int                     versions[1];
                } v;
        } u;
} __attribute__((__packed__));

/****************************************************************************/
struct cmd_data {
        struct  list_head       cmd_entry;
        struct  rpc_cmd         cmd;
};

struct cmd_storage {
        int     num_cmd;
        int     size_hash;
        /* List of rpc_cmd objects */
        struct list_head        cmd_hash[1];
};


struct cmd_storage *cmd_s_alloc(int max_versions);
int cmd_s_free(struct cmd_storage *s);
int cmd_s_find_all(struct cmd_storage *s, struct obj_descriptor *odsc,
                struct list_head *out_list, int *out_num_cmd);
struct cmd_data *
cmd_s_find_no_version(struct cmd_storage *s, struct obj_descriptor *odsc);
void cmd_s_add(struct cmd_storage *s, struct rpc_cmd *cmd);

/****************************************************************************/

#ifdef __cplusplus
}
#endif

#endif
