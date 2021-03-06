#ifndef LB3D_PRIME_H
#define LB3D_PRIME_H
//##############################################################################
//
// lbgpu_prime.h
//
//  - Header file for lbgpu_prime.
//

#include <stdio.h>
#include <stdlib.h>
//#include <malloc.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifndef CLK_TCK
// Old versions of gcc use CLK_TCK.
#define CLK_TCK CLOCKS_PER_SEC
#endif

//#include "flags.h"   //included in lbgpu_prime.cu
//#include "bc_flags.h"
#include "process.h"
#include "lattice.h"   //included in process.h
#include "forward_declarations.h"

// D3Q19
////
////             0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
////             |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
//int vx[19] = { 0,-1, 1, 0, 0, 0, 0,-1, 1,-1, 1,-1, 1,-1, 1, 0, 0, 0, 0};
//int vy[19] = { 0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1};
//int vz[19] = { 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1};
////             C  W  E  N  S  T  B  N  N  S  S  T  T  B  B  T  T  B  B
////                                  W  E  W  E  W  E  W  E  N  S  N  S

#if 0
// A better configuration - minimizes the number of memcpys required
// to fill the MPI buffers. D2Q9 still contiguous at the beginning,
// and opposites separated by a stride (4 starting from E and 5 starting
// from T. Before you implement this, do a test to see how much of a time 
// reduction is involved... Yes there is a significant reduction, but this
// is too simplistic - we're only transferring boundary cells from a few
// directions, so even if those directions are adjecent in the f array,
// we won't be able to do the transfer with a single memcpy.
//             0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
//             |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
int vx[19] = { 0, 1, 0, 1,-1,-1, 0,-1, 1, 0, 1,-1, 0, 0, 0,-1, 1, 0, 0};
int vy[19] = { 0, 0, 1, 1, 1, 0,-1,-1,-1, 0, 0, 0, 1,-1, 0, 0, 0,-1, 1};
int vz[19] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1};
//             C  E  N  N  N  W  S  S  S  T  T  T  T  T  B  B  B  B  B
//                      E  W        W  E     E  W  N  S     W  E  S  N

#define C   0
#define E   1
#define N   2
#define NE  3
#define NW  4
#define W   5
#define S   6
#define SW  7
#define SE  8
#define T   9
#define TE 10
#define TW 11
#define TN 12
#define TS 13
#define B  14
#define BW 15
#define BE 16
#define BS 17
#define BN 18

#endif

// Rearranged so that opposing directions are adjacent and the D2Q9 directions
// are contiguous at the beginning.
//             0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
//             |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
int vx[19] = { 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0};
int vy[19] = { 0, 0, 0, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1};
int vz[19] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1};
//             C  E  W  N  S  N  S  N  S  T  B  T  B  T  B  T  B  T  B
//                            E  W  W  E        E  W  W  E  N  S  S  N

real wt[19];
int cumul_stride[20];

#define C  0
#define E  1
#define W  2
#define N  3
#define S  4
#define NE 5
#define SW 6
#define NW 7
#define SE 8
#define T  9
#define B  10
#define TE 11
#define BW 12
#define TW 13
#define BE 14
#define TN 15
#define BS 16
#define TS 17
#define BN 18

#define EPS .0000000000001
#define PI  3.1415926535

#define HRULE0 "- - - - - - - - - - - - - - - - - - - - " \
               "- - - - - - - - - - - - - - - - - - - - "

#define HRULE1 "----------------------------------------" \
               "----------------------------------------"

#define HRULE2 "========================================" \
               "========================================"

#if DO_NOT_STORE_SOLIDS
// TODO: Incorporate version that omits storage of interior solids (which
// are not involved in flow).
//#include "min_nodes/compute.c"
//#include "min_nodes/stream.c"
//#include "min_nodes/bcs.c"
//#include "min_nodes/collide.c"
//#include "min_nodes/lbio.c"
//#include "min_nodes/latman.c"
#else /* !( DO_NOT_STORE_SOLIDS) */
#include "lattice.c"
#include "params.h"
#include "process.c"
#include "compute.c"
//#include "stream.c"
#ifdef __CUDACC__
#include "boundary_kernels.c"
#include "system_boundary_kernels.c"
#include "k_stream_collide_stream.c"
#include "k_collide.c"
#else
#include "stream_collide_stream.c"
#include "collide.c"
#endif
#include "bcs.c"
#include "lbio.c"
#include "latman.c"
#endif /* DO_NOT_STORE_SOLIDS */

#endif /* LB3D_PRIME_H */
