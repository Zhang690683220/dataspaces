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
*  Fan Zhang (2013) TASSL Rutgers University
*  zhangfan@cac.rutgers.edu
*/
#ifndef __DART_RDMA_UGNI_H__
#define __DART_RDMA_UGNI_H__

#include "config.h"

#ifdef DS_HAVE_DIMES
#include "dart_rpc_gni.h"

enum dart_memory_type {
    dart_memory_non_rdma = 0,
    dart_memory_rdma,
#ifdef DS_HAVE_DIMES_SHMEM
    dart_memory_shmem_non_rdma,
    dart_memory_shmem_rdma,
#endif
};

enum dart_rdma_tran_type {
    dart_rdma_tran_type_read = 0,
    dart_rdma_tran_type_write
};

struct dart_rdma_mem_handle {
	gni_mem_handle_t mdh;
	void *base_addr;
	size_t size;
    enum dart_memory_type mem_type;
};

struct dart_rdma_op {
	struct list_head entry;
	int tran_id;
	gni_post_descriptor_t post_desc;
	size_t src_offset;
	size_t dst_offset;
	size_t bytes;
	int ret; // TODO(fan): Need this?
};

struct dart_rdma_tran {
	struct list_head entry;
	struct list_head rdma_op_list;
	int tran_id;
	int num_rdma_op;
	int f_all_rdma_op_posted;
    enum dart_rdma_tran_type type;
	struct node_id *remote_peer;
	struct dart_rdma_mem_handle src;
	struct dart_rdma_mem_handle dst;
};

struct dart_rdma_handle {
	struct rpc_server *rpc_s;
	gni_cq_handle_t post_cq_hndl;
	struct list_head read_tran_list;
    struct list_head write_tran_list;
	int num_rdma_op_posted;
};

int dart_rdma_init(struct rpc_server *rpc_s);
int dart_rdma_finalize();

int dart_rdma_register_mem(struct dart_rdma_mem_handle *mem_hndl,
                           void *data, size_t bytes);
int dart_rdma_deregister_mem(struct dart_rdma_mem_handle *mem_hndl);

int dart_rdma_set_memregion_to_cmd(struct dart_rdma_mem_handle *mem_hndl,
                                   struct rpc_cmd *cmd);
int dart_rdma_get_memregion_from_cmd(struct dart_rdma_mem_handle *mem_hndl,
                                   struct rpc_cmd *cmd);

int dart_rdma_create_read_tran(struct node_id *remote_peer,
                struct dart_rdma_tran **pp);
int dart_rdma_delete_read_tran(int tran_id);

int dart_rdma_schedule_read(int tran_id, size_t src_offset, size_t dst_offset,
                size_t bytes);
int dart_rdma_perform_reads(int tran_id);
int dart_rdma_process_reads();
int dart_rdma_check_reads(int tran_id);

#ifdef DS_HAVE_DIMES_SHMEM
int dart_rdma_perform_reads_local(int tran_id);
#endif

/*
 DART gni as-a-service implementation does not create the GNI end point for peers
 of remote application. Because GNI end point is required to post one-sided GET 
 operation on RDMA memory buffer created by remote peer, we added the following 
 APIs to dynamically create/bind and unbind/delete GNI end point for peer of 
 remote applications.
*/
struct node_id* dart_rdma_create_remote_peer(struct ptlid_map *ptlmap);
int dart_rdma_delete_remote_peer(struct node_id* peer);
#endif

#endif
