//##############################################################################
//
// lbgpu_prime.c
//
//  - Lattice Boltzmann
//
//  - D3Q19
//

typedef double real; // To be relocated (in flags.in?).

#ifdef __CUDACC__
  real* f_mem_d;
  int f_mem_size;
  real* mv_mem_d;
  int mv_mem_size;
  unsigned char* solids_mem_d;
  int solids_mem_size;
  int* is_end_of_frame_mem_d;
#endif

#include "lbgpu_prime.h"

int main( int argc, char **argv)
{
  int subs;
  int time, frame;

  struct lattice_struct *lattice;

  setbuf( stdout, (char*)NULL); // Don't buffer screen output.

#if VERBOSITY_LEVEL > 0
  printf("%s %d >> lbgpu_prime.c: main() -- Hello.\n", __FILE__, __LINE__);
#endif /* VERBOSITY_LEVEL > 0 */

  construct_lattice( &lattice, argc, argv);

  process_tic( lattice);

  init_problem( lattice);

  set_time( lattice, 0);
  set_frame( lattice, 0);

  output_frame( lattice);

  process_toc( lattice);
  display_etime( lattice);
  process_tic( lattice);

#ifdef __CUDACC__
  int blocksize = 3; // TODO: What is the best value for this?

  if( get_NumNodes( lattice) % blocksize)
  {
    printf("\n%s %d ERROR: NumNodes %d must be divisible by blocksize %d. "
           "NumNodes %% blocksize = %d. (Exiting.)\n\n"
          , __FILE__
          , __LINE__
          , get_NumNodes( lattice)
          , blocksize
          , get_NumNodes( lattice) % blocksize);
    process_exit(1);
  }

  dim3 blockDim( blocksize, 1, 1);
  dim3 gridDim( get_NumNodes( lattice) / blocksize, 1, 1);
#endif

#ifdef __CUDACC__
#if 0
  cudaMemcpy( f_mem_d
            , get_fptr(lattice, 0, 0, 0, 0, 0)
            , get_NumVelDirs(lattice)
             *get_NumNodes(lattice)
             *get_NumSubs(lattice)
             *sizeof(real)
            , cudaMemcpyHostToDevice);
#else
  cudaMemcpy( solids_mem_d
            , get_solids_ptr(lattice, 0)
            , get_NumNodes(lattice)*sizeof(unsigned char)
            , cudaMemcpyHostToDevice);

  for( subs = 0; subs < get_NumSubs(lattice); subs++)
  {
    cudaMemcpy( f_mem_d + subs*get_NumNodes(lattice)*get_NumVelDirs(lattice)
              , get_fptr(lattice, subs, 0,0,0, 0,0,0, 0)
              , get_NumVelDirs(lattice)
               *get_NumNodes(lattice)
               *sizeof(real)
              , cudaMemcpyHostToDevice);
  }

  int temp = 0;
  cudaMemcpy( is_end_of_frame_mem_d
            , &temp
            , sizeof(int)
            , cudaMemcpyHostToDevice);
#endif
#endif

  for( frame = 0, time=1; time<=get_NumTimeSteps(lattice); time++)
  {
    set_time( lattice, time);

    // Time steps are combined in a somewhat awkward way in order to minimize
    // the amount of temporary storage required and to allow the nodes to
    // be updated independently of neighbor nodes. This facilitates an
    // efficient GPU implementation.
    //
    // c.f., Bailey, Myre, Walsh et. al., Accelerating Lattice Boltzmann Fluid
    // Flow Simulations using Graphics Processors, ICPP 2009.
    // http://www.tc.umn.edu/~bail0253/

#ifdef __CUDACC__
    k_stream_collide_stream
    <<<
        gridDim
      , blockDim
      , blocksize*sizeof(real)*( get_NumVelDirs(lattice)
                               + get_NumDims(lattice)+1 )
    >>>( f_mem_d, mv_mem_d, solids_mem_d);
#else
    stream_collide_stream(lattice);
#endif

    set_time( lattice, ++time);

#ifdef __CUDACC__
    if( is_end_of_frame(lattice,time))
    {
      int temp = 1;
      cudaMemcpy( is_end_of_frame_mem_d
                , &temp
                , sizeof(int)
                , cudaMemcpyHostToDevice);
    }
    k_collide
    <<<
        gridDim
      , blockDim
      , blocksize*sizeof(real)*( get_NumVelDirs(lattice)
                               + get_NumDims(lattice)+1 )
    >>>( f_mem_d, mv_mem_d, solids_mem_d, is_end_of_frame_mem_d);
    if( is_end_of_frame(lattice,time))
    {
      int temp = 0;
      cudaMemcpy( is_end_of_frame_mem_d
                , &temp
                , sizeof(int)
                , cudaMemcpyHostToDevice);
    }
#else
    collide( lattice);
#endif

    if( is_end_of_frame(lattice,time))
    {
      set_frame( lattice, ++frame);

#ifdef __CUDACC__
      for( subs = 0; subs < get_NumSubs(lattice); subs++)
      {
        printf("Transferring subs %d "
               "macrovars from device to host for output. \n", subs);
        cudaMemcpy( get_rho_ptr(lattice, subs, 0)
                  , mv_mem_d
                    + subs*get_NumNodes(lattice)*( 1 + get_NumDims(lattice))
                    + 0*get_NumNodes(lattice)
                  , get_NumNodes(lattice)*sizeof(real)
                  , cudaMemcpyDeviceToHost);
        cudaMemcpy( get_ux_ptr(lattice, subs, 0)
                  , mv_mem_d
                    + subs*get_NumNodes(lattice)*( 1 + get_NumDims(lattice))
                    + 1*get_NumNodes(lattice)
                  , get_NumNodes(lattice)*sizeof(real)
                  , cudaMemcpyDeviceToHost);
        cudaMemcpy( get_uy_ptr(lattice, subs, 0)
                  , mv_mem_d
                    + subs*get_NumNodes(lattice)*( 1 + get_NumDims(lattice))
                    + 2*get_NumNodes(lattice)
                  , get_NumNodes(lattice)*sizeof(real)
                  , cudaMemcpyDeviceToHost);
       if(get_NumDims(lattice) == 3)
       {
        cudaMemcpy( get_uz_ptr(lattice, subs, 0)
                  , mv_mem_d
                    + subs*get_NumNodes(lattice)*( 1 + get_NumDims(lattice))
                    + 3*get_NumNodes(lattice)
                  , get_NumNodes(lattice)*sizeof(real)
                  , cudaMemcpyDeviceToHost);
       }
       if( do_write_f_txt_file( lattice))
       {
        cudaMemcpy( get_fptr(lattice, subs, 0,0,0, 0,0,0, 0)
                  , f_mem_d + subs*get_NumNodes(lattice)*get_NumVelDirs(lattice)
                  , get_NumVelDirs(lattice)
                   *get_NumNodes(lattice)
                   *sizeof(real)
                  , cudaMemcpyDeviceToHost);
       }
      }
#endif

      output_frame( lattice);
      process_toc( lattice);
      display_etime( lattice);

#ifdef __CUDACC__
      printf("**************\n");
      printf("  __CUDACC__\n");
      printf("**************\n");
#endif
    }

  } /* for( time=1; time<=lattice->NumTimeSteps; time++) */

  process_barrier();
  process_toc( lattice);
  display_etime( lattice);

  destruct_lattice( lattice);

#if VERBOSITY_LEVEL > 0
  printf("\n");
  printf("%s %d >> lbgpu_prime.c: main() -- Terminating normally.\n",
      __FILE__, __LINE__);
  printf("\n");
#endif /* VERBOSITY_LEVEL > 0 */

  return 0;

} /* int main( int argc, char **argv) */
