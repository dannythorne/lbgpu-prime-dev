//##############################################################################
//
// lbgpu_prime.cu
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

  if( get_LX( lattice) % get_BX(lattice) != 0)
  {
    printf("\n%s %d ERROR: LX = %d must be divisible by BX = %d. (Exiting.)\n\n"
          , __FILE__
          , __LINE__
          , get_LX( lattice)
          , get_BX( lattice));
    process_exit(1);
  }
  if( get_LY( lattice) % get_BY(lattice) != 0)
  {
    printf("\n%s %d ERROR: LY = %d must be divisible by BY = %d. (Exiting.)\n\n"
          , __FILE__
          , __LINE__
          , get_LY( lattice)
          , get_BY( lattice));
    process_exit(1);
  }
  if( get_LZ( lattice) % get_BZ(lattice) != 0)
  {
    printf("\n%s %d ERROR: LX = %d must be divisible by BX = %d. (Exiting.)\n\n"
          , __FILE__
          , __LINE__
          , get_LZ( lattice)
          , get_BZ( lattice));
    process_exit(1);
  }


  dim3 blockDim( get_BX( lattice), get_BY( lattice), get_BZ( lattice));
  dim3 gridDim( get_LX( lattice) / get_BX( lattice),
                get_LY( lattice) / get_BY( lattice),
                get_LZ( lattice) / get_BZ( lattice) );

#endif

#ifdef __CUDACC__
#if 0
  cudaError_t err; 

  cudaMemcpy( f_mem_d
            , get_fptr(lattice, 0,0,0,0, 0, 0, 0, 0)
            , get_NumVelDirs(lattice)
             *get_NumNodes(lattice)
             *get_NumSubs(lattice)
             *sizeof(real)
            , cudaMemcpyHostToDevice);
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");
#else
  cudaMemcpy( solids_mem_d
            , get_solids_ptr(lattice, 0)
            , get_NumNodes(lattice)*sizeof(unsigned char)
            , cudaMemcpyHostToDevice);

  checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

  for( subs = 0; subs < get_NumSubs(lattice); subs++)
  {
    cudaMemcpy( f_mem_d + subs*get_NumNodes(lattice)*get_NumVelDirs(lattice)
              , get_fptr(lattice, subs, 0,0,0, 0,0,0, 0)
              , get_NumVelDirs(lattice)
               *get_NumNodes(lattice)
               *sizeof(real)
              , cudaMemcpyHostToDevice);

    checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

    cudaMemcpy( mv_mem_d + subs*get_NumNodes(lattice)*(1 + get_NumDims(lattice))
              , get_rho_ptr(lattice, subs, 0)
              , (1. + get_NumDims(lattice))
               *get_NumNodes(lattice)
               *sizeof(real)
              , cudaMemcpyHostToDevice);

    checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

  }

  int temp = 0;
  cudaMemcpy( is_end_of_frame_mem_d
            , &temp
            , sizeof(int)
            , cudaMemcpyHostToDevice);

  checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

#endif
  
#if TEXTURE_FETCH
  cudaBindTexture(0, tex_solid, solids_mem_d);
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
      , sizeof(real)*( get_NumVelDirs(lattice)
                    + get_NumDims(lattice)+1 )
                    * get_BX( lattice)
                    * get_BY( lattice)
                    * get_BZ( lattice)
    >>>( f_mem_d, mv_mem_d, solids_mem_d);
    cudaThreadSynchronize();
    checkCUDAError( __FILE__, __LINE__, "k_stream_collide_stream");
    
#else
    stream_collide_stream(lattice);
#endif

    set_time( lattice, ++time);

#ifdef __CUDACC__
#if 1
    if( is_end_of_frame(lattice,time))
    {
      int temp = 1;
      cudaMemcpy( is_end_of_frame_mem_d
                , &temp
                , sizeof(int)
                , cudaMemcpyHostToDevice);

      checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

    }
#endif
    k_collide
    <<<
        gridDim
      , blockDim
      , sizeof(real)*( get_NumVelDirs(lattice)
                    + get_NumDims(lattice)+1 )
                    * get_BX( lattice)
                    * get_BY( lattice)
                    * get_BZ( lattice)
    >>>( f_mem_d, mv_mem_d, solids_mem_d, is_end_of_frame_mem_d);
    cudaThreadSynchronize();
    checkCUDAError( __FILE__, __LINE__, "k_collide");

#if 1
    if( is_end_of_frame(lattice,time))
    {
      int temp = 0;
      cudaMemcpy( is_end_of_frame_mem_d
                , &temp
                , sizeof(int)
                , cudaMemcpyHostToDevice);
      checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

    }
#endif
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
        checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        cudaMemcpy( get_ux_ptr(lattice, subs, 0)
                  , mv_mem_d
                    + subs*get_NumNodes(lattice)*( 1 + get_NumDims(lattice))
                    + 1*get_NumNodes(lattice)
                  , get_NumNodes(lattice)*sizeof(real)
                  , cudaMemcpyDeviceToHost);
        checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        cudaMemcpy( get_uy_ptr(lattice, subs, 0)
                  , mv_mem_d
                    + subs*get_NumNodes(lattice)*( 1 + get_NumDims(lattice))
                    + 2*get_NumNodes(lattice)
                  , get_NumNodes(lattice)*sizeof(real)
                  , cudaMemcpyDeviceToHost);
        checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        if(get_NumDims(lattice) == 3)
        {
          cudaMemcpy( get_uz_ptr(lattice, subs, 0)
                    , mv_mem_d
                    + subs*get_NumNodes(lattice)*( 1 + get_NumDims(lattice))
                    + 3*get_NumNodes(lattice)
                    , get_NumNodes(lattice)*sizeof(real)
                    , cudaMemcpyDeviceToHost);
          checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        }
        if( do_write_f_txt_file( lattice))
        {
          cudaMemcpy( get_fptr(lattice, subs, 0,0,0, 0,0,0, 0)
                    , f_mem_d + subs*get_NumNodes(lattice)*get_NumVelDirs(lattice)
                    , get_NumVelDirs(lattice)
                    *get_NumNodes(lattice)
                    *sizeof(real)
                    , cudaMemcpyDeviceToHost);
          checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

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
