//##############################################################################
//
// process.h
//

#ifndef PROCESS_H
#define PROCESS_H

#if PARALLEL
#include <mpi.h>

#if DP_ON
#define MPI_DATATYPE MPI_DOUBLE
#else
#define MPI_DATATYPE MPI_FLOAT
#endif

#endif

#include "lattice.h"

void process_init( lattice_ptr lattice, int argc, char **argv);
void process_compute_local_params( lattice_ptr lattice);
void process_send_recv_begin( lattice_ptr lattice, const int subs);
void process_send_recv_end( lattice_ptr lattice, const int subs);
void rho_send_recv_begin( lattice_ptr lattice, const int subs);
void rho_send_recv_end( lattice_ptr lattice, const int subs);
void solid_send_recv_begin( lattice_ptr lattice, const int subs);
void solid_send_recv_end( lattice_ptr lattice, const int subs);
void process_dump_pdfs_to_recv( lattice_ptr lattice, char *comment_str);
void process_dump_pdfs_to_send( lattice_ptr lattice, char *comment_str);
void process_dump_pdfs( lattice_ptr lattice, char *comment_str, real *pdfs);
void gather_north_pointing_pdfs(
       lattice_ptr lattice,
       real *north,
       const int subs,
       const int k,
       const int which_pdf);
void gather_south_pointing_pdfs(
       lattice_ptr lattice,
       real *south,
       const int subs,
       const int k,
       const int which_pdf);
void process_reduce_real_sum( lattice_ptr lattice, real *arg_x);
void process_reduce_int_sum( lattice_ptr lattice, int *arg_n);
void process_barrier();
void process_finalize();
void process_exit( int exit_val);

//##############################################################################
//
// Accessor methods for process struct.
//

int get_proc_id( lattice_ptr lattice);
int get_num_procs( lattice_ptr lattice);
int* get_proc_id_ptr( lattice_ptr lattice);
int* get_num_procs_ptr( lattice_ptr lattice);

int is_on_root_proc( lattice_ptr lattice);

int get_g_LX( lattice_ptr lattice);
int get_g_LY( lattice_ptr lattice);
int get_g_LZ( lattice_ptr lattice);
int get_g_SX( lattice_ptr lattice);
int get_g_SY( lattice_ptr lattice);
int get_g_SZ( lattice_ptr lattice);
int get_g_EX( lattice_ptr lattice);
int get_g_EY( lattice_ptr lattice);
int get_g_EZ( lattice_ptr lattice);
int get_g_NumNodes( lattice_ptr lattice);
int get_g_StartNode( lattice_ptr lattice);

#if PARALLEL
int get_g_SX( lattice_ptr lattice);
int get_g_SY( lattice_ptr lattice);
int get_g_SZ( lattice_ptr lattice);
int get_g_EX( lattice_ptr lattice);
int get_g_EY( lattice_ptr lattice);
int get_g_EZ( lattice_ptr lattice);

void set_g_LX( lattice_ptr lattice, const int arg_LX);
void set_g_LY( lattice_ptr lattice, const int arg_LY);
void set_g_LZ( lattice_ptr lattice, const int arg_LZ);
void set_g_SX( lattice_ptr lattice, const int arg_SX);
void set_g_SY( lattice_ptr lattice, const int arg_SY);
void set_g_SZ( lattice_ptr lattice, const int arg_SZ);
void set_g_EX( lattice_ptr lattice, const int arg_EX);
void set_g_EY( lattice_ptr lattice, const int arg_EY);
void set_g_EZ( lattice_ptr lattice, const int arg_EZ);

void set_g_NumNodes( lattice_ptr lattice, const int arg_NumNodes);
void set_g_StartNode( lattice_ptr lattice, const int arg_n);
#endif

#endif
