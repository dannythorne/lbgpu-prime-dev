//##############################################################################
//
// lbgpu_prime.cu
//
//  - Lattice Boltzmann
//
//  - D2Q9 / D3Q19
//
#include "flags.h"

#if DP_ON
typedef double real; // To be relocated (in flags.in?).
#else
typedef float real; // To be relocated (in flags.in?).
#endif

#ifdef __CUDACC__
real* f_mem_d;
int f_mem_size;
real* mv_mem_d;
int mv_mem_size;
unsigned char* solids_mem_d;
real* ns_mem_d;
int solids_mem_size;

#if BOUNDARY_KERNEL
real* pos_dir_send_ptr_d;
real* pos_dir_recv_ptr_d;
real* neg_dir_send_ptr_d;
real* neg_dir_recv_ptr_d;
#endif
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

#if 1
  if( get_LX( lattice) % get_BX( lattice) != 0)
  {
    printf("\n%s %d ERROR: LX = %d must be divisible by BX = %d. (Exiting.)\n\n"
        , __FILE__
        , __LINE__
        , get_LX( lattice)
        , get_BX( lattice));
    process_exit(1);
  }
  if( get_LY( lattice) % get_BY( lattice) != 0)
  {
    printf("\n%s %d ERROR: LY = %d must be divisible by BY = %d. (Exiting.)\n\n"
        , __FILE__
        , __LINE__
        , get_LY( lattice)
        , get_BY( lattice));
    process_exit(1);
  }
  if( get_LZ( lattice) % get_BZ( lattice) != 0)
  {
    printf("\n%s %d ERROR: LX = %d must be divisible by BX = %d. (Exiting.)\n\n"
        , __FILE__
        , __LINE__
        , get_LZ( lattice)
        , get_BZ( lattice));
    process_exit(1);
  }
#endif

  dim3 blockDim( get_BX( lattice), get_BY( lattice), get_BZ( lattice));
  dim3 gridDim( get_LX( lattice) / get_BX( lattice),
      get_LY( lattice) / get_BY( lattice),
      get_LZ( lattice) / get_BZ( lattice) );

#if BOUNDARY_KERNEL  

  dim3 blockboundDim(1, 1, 1);
  dim3 gridboundDim(1, 1, 1);

  if( get_NumDims( lattice) == 2)
  { 
    blockboundDim.x = get_BX( lattice); 
    gridboundDim.x = get_LX( lattice) / get_BX( lattice); 
  }
  if( get_NumDims( lattice) == 3)
  {  
    blockboundDim.x = get_BX( lattice); 
    gridboundDim.x = get_LX( lattice) / get_BX( lattice); 
    blockboundDim.y = get_BY( lattice); 
    gridboundDim.y = get_LY( lattice) / get_BY( lattice); 
  }


#endif

  cudaMemcpy( solids_mem_d + get_EndBoundSize( lattice)
      , get_solids_ptr(lattice, 0)
      , get_NumNodes( lattice)*sizeof(unsigned char)
      , cudaMemcpyHostToDevice);

  checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

  // Do boundary swap for solids.  This only needs to be done
  // once, of course.

#if PARALLEL
  // We'll want something involving the MPI buffers here.

  solid_send_recv_begin( lattice, 0);
  solid_send_recv_end( lattice, 0);

  // North / Top end
  cudaMemcpy( solids_mem_d + get_EndBoundSize( lattice) + get_NumNodes( lattice)
      , lattice->process.neg_dir_solid_to_recv
      , get_EndBoundSize( lattice)*sizeof(unsigned char)
      , cudaMemcpyHostToDevice);

  checkCUDAError( __FILE__, __LINE__, "solid swap");

  // South / Bottom end
  cudaMemcpy( solids_mem_d
      , lattice->process.pos_dir_solid_to_recv
      , get_EndBoundSize( lattice)*sizeof(unsigned char)
      , cudaMemcpyHostToDevice);

  checkCUDAError( __FILE__, __LINE__, "solid swap");

#else   // if !(PARALLEL)
  // North / Top end
  cudaMemcpy( solids_mem_d + get_EndBoundSize( lattice) + get_NumNodes( lattice)
      , solids_mem_d + get_EndBoundSize( lattice)
      , get_EndBoundSize( lattice)*sizeof(unsigned char)
      , cudaMemcpyDeviceToDevice);

  checkCUDAError( __FILE__, __LINE__, "solid swap");

  // South / Bottom end
  cudaMemcpy( solids_mem_d
      , solids_mem_d + get_NumNodes( lattice)
      , get_EndBoundSize( lattice)*sizeof(unsigned char)
      , cudaMemcpyDeviceToDevice);

  checkCUDAError( __FILE__, __LINE__, "solid swap");
#endif  // if PARALLEL

  // Same as above, but for the array of ns values
#if WALSH_NS_ON
  cudaMemcpy( ns_mem_d + get_EndBoundSize( lattice)
      , get_ns_ptr(lattice, 0)
      , get_NumNodes( lattice)*sizeof(real)
      , cudaMemcpyHostToDevice);

  checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

  // Do boundary swap for ns.  This only needs to be done
  // once, of course.

#if PARALLEL
  // We'll want something involving the MPI buffers here.

  ns_send_recv_begin( lattice, 0);
  ns_send_recv_end( lattice, 0);

  // North / Top end
  cudaMemcpy( ns_mem_d + get_EndBoundSize( lattice) + get_NumNodes( lattice)
      , lattice->process.neg_dir_ns_to_recv
      , get_EndBoundSize( lattice)*sizeof(real)
      , cudaMemcpyHostToDevice);

  checkCUDAError( __FILE__, __LINE__, "ns swap");

  // South / Bottom end
  cudaMemcpy( ns_mem_d
      , lattice->process.pos_dir_ns_to_recv
      , get_EndBoundSize( lattice)*sizeof(real)
      , cudaMemcpyHostToDevice);

  checkCUDAError( __FILE__, __LINE__, "ns swap");

#else   // if !(PARALLEL)
  // North / Top end
  cudaMemcpy( ns_mem_d + get_EndBoundSize( lattice) + get_NumNodes( lattice)
      , ns_mem_d + get_EndBoundSize( lattice)
      , get_EndBoundSize( lattice)*sizeof(real)
      , cudaMemcpyDeviceToDevice);

  checkCUDAError( __FILE__, __LINE__, "ns swap");

  // South / Bottom end
  cudaMemcpy( ns_mem_d
      , ns_mem_d + get_NumNodes( lattice)
      , get_EndBoundSize( lattice)*sizeof(real)
      , cudaMemcpyDeviceToDevice);

  checkCUDAError( __FILE__, __LINE__, "ns swap");
#endif  // if PARALLEL

#endif


  int dir;
#if 1
  for( subs = 0; subs < get_NumSubs( lattice); subs++)
  {
    // In the final version, the CPU arrays may well have
    // the same 'additional boundary' structure as the GPU
    // arrays, in which case a single memcpy can be performed.
    cudaMemcpy( f_mem_d + subs * cumul_stride[get_NumVelDirs( lattice)]
        , get_fptr(lattice, subs, 0,0,0, 0,0,0, 0)
        , get_NumNodes( lattice)
        *(get_NumUnboundDirs( lattice))
        *sizeof(real)
        , cudaMemcpyHostToDevice);

    checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

    for( dir=get_NumUnboundDirs( lattice); dir < get_NumVelDirs( lattice); dir+=2)
    {

      cudaMemcpy( f_mem_d + cumul_stride[dir]
          + subs * cumul_stride[get_NumVelDirs( lattice)]
          , get_fptr(lattice, subs, 0,0,0, 0,0,0, dir)
          , 2*get_NumNodes( lattice)
          *sizeof(real)
          , cudaMemcpyHostToDevice);

      checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");
    }


    cudaMemcpy( mv_mem_d + subs*get_NumNodes( lattice)*(1 + get_NumDims( lattice))
        , get_rho_ptr(lattice, subs, 0)
        , (1. + get_NumDims( lattice))
        *get_NumNodes( lattice)
        *sizeof(real)
        , cudaMemcpyHostToDevice);

    checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

  }
#endif
  int temp_frame_switch = 0;
  cudaMemcpyToSymbol( is_end_of_frame_mem_c, &temp_frame_switch, sizeof(int));

#if TEXTURE_FETCH
  cudaBindTexture(0, tex_solid, solids_mem_d);
#endif

#endif  // ifdef __CUDACC__


  //cudaEvent_t start, stop;
  //real timertime;
  //real totaltime = 0.;
  //cudaEventCreate(&start);
  //cudaEventCreate(&stop);

  read_PEST_in_files( &lattice, argc, argv);

  for( frame = 0, time=1; time<=get_NumTimeSteps( lattice); time++)
  {
    set_time( lattice, time);


#if INAMURO_SIGMA_COMPONENT
    cudaMemcpyToSymbol( time_c, &(lattice->time), sizeof(int));
    checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
#endif

    // TODO: All of the boundary swap stuff should be part of the
    // mpi sendrecv functions.  For the case of CUDA without MPI, 
    // the periodic boundary conditions implemented within get_f1
    // should be restored using #if !(PARALLEL)
    // Do boundary swaps
#if PARALLEL
#ifdef __CUDACC__
    //cudaEventRecord(start,0);
#if BOUNDARY_KERNEL
    k_bound_DtH_1
      <<<
      gridboundDim
      , blockboundDim
      >>>( f_mem_d, pos_dir_send_ptr_d, neg_dir_send_ptr_d);

    cudaThreadSynchronize();
    checkCUDAError( __FILE__, __LINE__, "boundary kernel");

#if !(POINTER_MAPPING)
    cudaMemcpy( lattice->process.pos_dir_pdf_to_send
        , pos_dir_send_ptr_d
        , get_NumBoundDirs( lattice) * get_EndBoundSize( lattice)
        * get_NumSubs( lattice) * sizeof(real) / 2
        , cudaMemcpyDeviceToHost);
    checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
    cudaMemcpy( lattice->process.neg_dir_pdf_to_send
        , neg_dir_send_ptr_d
        , get_NumBoundDirs( lattice) * get_EndBoundSize( lattice)
        * get_NumSubs( lattice) * sizeof(real) / 2
        , cudaMemcpyDeviceToHost);
    checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
#endif  // !(POINTER_MAPPING)
#else   // !(BOUNDARY_KERNEL)
    for( subs = 0; subs < get_NumSubs( lattice); subs++)
    {
#if 0
      for( dir=get_NumUnboundDirs( lattice); dir < get_NumVelDirs( lattice); dir++)
      {
        pdf_boundary_parallel( lattice, f_mem_d, cumul_stride
            , subs, dir, time, DEVICE_TO_HOST);
      }
#else
      // For D2Q9 and D3Q19, NumUnboundDirs is odd
      for( dir=get_NumUnboundDirs( lattice); dir < get_NumVelDirs( lattice); dir+=2)
      {
        pdf_bound_DtH_odd_1( lattice, f_mem_d, cumul_stride
            , subs, dir);
      }
      for( dir=get_NumUnboundDirs( lattice)+1; dir < get_NumVelDirs( lattice); dir+=2)
      {
        pdf_bound_DtH_even_1( lattice, f_mem_d, cumul_stride
            , subs, dir);
      }
#endif
    }

#endif  //BOUNDARY_KERNEL
    //cudaEventRecord(stop, 0);
    //cudaEventSynchronize(stop);
    //cudaEventElapsedTime(&timertime,start,stop);
    //totaltime += timertime;
#endif  // #ifdef __CUDACC__
    // Buffers now large enough to deal with all substances in one transfer.
    // Currently, send a receives are activated if __CUDACC__ not defined,
    // but that's bad because for that case they're totally broken.
    process_send_recv_begin( lattice, 0);
    process_send_recv_end(lattice, 0);

#ifdef __CUDACC__
    //cudaEventRecord(start,0);
#if BOUNDARY_KERNEL
#if !(POINTER_MAPPING)
    cudaMemcpy( pos_dir_recv_ptr_d
        , lattice->process.pos_dir_pdf_to_recv
        , get_NumBoundDirs( lattice) * get_EndBoundSize( lattice)
        * get_NumSubs( lattice) * sizeof(real) / 2
        , cudaMemcpyHostToDevice);
    checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
    cudaMemcpy( neg_dir_recv_ptr_d
        , lattice->process.neg_dir_pdf_to_recv
        , get_NumBoundDirs( lattice) * get_EndBoundSize( lattice)
        * get_NumSubs( lattice) * sizeof(real) / 2
        , cudaMemcpyHostToDevice);
    checkCUDAError( __FILE__, __LINE__, "boundary swap"); 

#endif  // !(POINTER_MAPPING)

    k_bound_HtD_1
      <<<
      gridboundDim
      , blockboundDim
      >>>( f_mem_d, pos_dir_recv_ptr_d, neg_dir_recv_ptr_d);

    cudaThreadSynchronize();
    checkCUDAError( __FILE__, __LINE__, "boundary kernel");
#else   // !(BOUNDARY_KERNEL)
    for( subs = 0; subs < get_NumSubs( lattice); subs++)
    {
#if 0
      for( dir=get_NumUnboundDirs( lattice); dir < get_NumVelDirs( lattice); dir++)
      {
        pdf_boundary_parallel( lattice, f_mem_d, cumul_stride
            , subs, dir, time, HOST_TO_DEVICE);
      }
#else
      // For D2Q9 and D3Q19, NumUnboundDirs is odd
      for( dir=get_NumUnboundDirs( lattice); dir < get_NumVelDirs( lattice); dir+=2)
      {
        pdf_bound_HtD_odd_1( lattice, f_mem_d, cumul_stride
            , subs, dir);
      }
      for( dir=get_NumUnboundDirs( lattice)+1; dir < get_NumVelDirs( lattice); dir+=2)
      {
        pdf_bound_HtD_even_1( lattice, f_mem_d, cumul_stride
            , subs, dir);
      }
#endif
    }

#endif  // BOUNDARY_KERNEL
    //cudaEventRecord(stop, 0);
    //cudaEventSynchronize(stop);
    //cudaEventElapsedTime(&timertime,start,stop);
    //totaltime += timertime;

#endif  // #ifdef __CUDACC__


#else   // if !(PARALLEL)
    // Implement this inside get_f1_d functions,
    // as per the old method
#if 0   
#ifdef __CUDACC__
    for( subs = 0; subs < get_NumSubs( lattice); subs++)
    {
      for( dir=get_NumUnboundDirs( lattice); dir < get_NumVelDirs( lattice); dir++)
      {
        pdf_boundary_swap( lattice, f_mem_d, cumul_stride
            , subs, dir, time);
      }
    }
#endif
#endif
#endif  // if (PARALLEL)

#if 1
    // Only implemented in 2D right now
    if( get_NumDims( lattice) == 2)
    {
      // Implement boundary conditions
      bcs_1( lattice, f_mem_d, solids_mem_d, ns_mem_d); 
    }
#endif

    // Time steps are combined in a somewhat awkward way in order to minimize
    // the amount of temporary storage required and to allow the nodes to
    // be updated independently of neighbor nodes. This facilitates an
    // efficient GPU implementation.
    //
    // c.f., Bailey, Myre, Walsh et. al., Accelerating Lattice Boltzmann Fluid
    // Flow Simulations using Graphics Processors, ICPP 2009.
    // http://www.tc.umn.edu/~bail0253/

#ifdef __CUDACC__
    if( is_end_of_frame(lattice,time))
    {
      temp_frame_switch = 1;
      cudaMemcpyToSymbol( is_end_of_frame_mem_c, &temp_frame_switch, sizeof(int));
    }

    PEST_gpu_switch( &lattice, mv_mem_d, argc, argv);

    k_stream_collide_stream
      <<<
      gridDim
      , blockDim
      , sizeof(real)*( get_NumVelDirs( lattice)
          + get_NumDims( lattice)+1 )
      * get_BX( lattice)
      * get_BY( lattice)
      * get_BZ( lattice)
      >>>( f_mem_d, mv_mem_d, solids_mem_d, ns_mem_d
         );

    cudaThreadSynchronize();
    checkCUDAError( __FILE__, __LINE__, "k_stream_collide_stream");

    write_PEST_out_data( &lattice, mv_mem_d, argc, argv);

    if( temp_frame_switch)
    {

      // Although it is necessary to copy device arrays to host at this point,
      // it is inefficient to write them to file here. Rather, we write to file
      // immediately after the first kernel, to allow concurrent host and
      // device execution.  On the last frame of course, it is necessary to
      // write to file here.
      set_frame( lattice, ++frame);

      for( subs = 0; subs < get_NumSubs( lattice); subs++)
      {
        printf("Transferring subs %d "
            "macrovars from device to host for output. \n", subs);
        cudaMemcpy( get_rho_ptr(lattice, subs, 0)
            , mv_mem_d
            + subs*get_NumNodes( lattice)*( 1 + get_NumDims( lattice))
            + 0*get_NumNodes( lattice)
            , get_NumNodes( lattice)*sizeof(real)
            , cudaMemcpyDeviceToHost);
        checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        cudaMemcpy( get_ux_ptr(lattice, subs, 0)
            , mv_mem_d
            + subs*get_NumNodes( lattice)*( 1 + get_NumDims( lattice))
            + 1*get_NumNodes( lattice)
            , get_NumNodes( lattice)*sizeof(real)
            , cudaMemcpyDeviceToHost);
        checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        cudaMemcpy( get_uy_ptr(lattice, subs, 0)
            , mv_mem_d
            + subs*get_NumNodes( lattice)*( 1 + get_NumDims( lattice))
            + 2*get_NumNodes( lattice)
            , get_NumNodes( lattice)*sizeof(real)
            , cudaMemcpyDeviceToHost);
        checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        if(get_NumDims( lattice) == 3)
        {
          cudaMemcpy( get_uz_ptr(lattice, subs, 0)
              , mv_mem_d
              + subs*get_NumNodes( lattice)*( 1 + get_NumDims( lattice))
              + 3*get_NumNodes( lattice)
              , get_NumNodes( lattice)*sizeof(real)
              , cudaMemcpyDeviceToHost);
          checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        }
      }
      output_frame( lattice);
      process_toc( lattice);
      display_etime( lattice);
      temp_frame_switch = 0;
      cudaMemcpyToSymbol( is_end_of_frame_mem_c, &temp_frame_switch, sizeof(int));
    }
#else   // ifndef __CUDACC__
    stream_collide_stream( lattice);
#endif  // ifdef __CUDACC__

    set_time( lattice, ++time);

#if INAMURO_SIGMA_COMPONENT
    cudaMemcpyToSymbol( time_c, &(lattice->time), sizeof(int));
    checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
#endif


    // Do boundary swaps.
#if PARALLEL
#ifdef __CUDACC__
    //cudaEventRecord(start, 0);
#if BOUNDARY_KERNEL
    k_bound_DtH_2
      <<<
      gridboundDim
      , blockboundDim
      >>>( f_mem_d, pos_dir_send_ptr_d, neg_dir_send_ptr_d);

    cudaThreadSynchronize();
    checkCUDAError( __FILE__, __LINE__, "boundary kernel");

#if !(POINTER_MAPPING)
    cudaMemcpy( lattice->process.pos_dir_pdf_to_send
        , pos_dir_send_ptr_d
        , get_NumBoundDirs( lattice) * get_EndBoundSize( lattice)
        * get_NumSubs( lattice) * sizeof(real) / 2
        , cudaMemcpyDeviceToHost);
    checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
    cudaMemcpy( lattice->process.neg_dir_pdf_to_send
        , neg_dir_send_ptr_d
        , get_NumBoundDirs( lattice) * get_EndBoundSize( lattice)
        * get_NumSubs( lattice) * sizeof(real) / 2
        , cudaMemcpyDeviceToHost);
    checkCUDAError( __FILE__, __LINE__, "boundary swap"); 

#endif  // !(POINTER_MAPPING)
#else   // !(BOUNDARY_KERNEL)
    for( subs = 0; subs < get_NumSubs( lattice); subs++)
    {
#if 0
      for( dir=get_NumUnboundDirs( lattice); dir < get_NumVelDirs( lattice); dir++)
      {
        pdf_boundary_parallel( lattice, f_mem_d, cumul_stride
            , subs, dir, time, DEVICE_TO_HOST);
      }
#else
      // For D2Q9 and D3Q19, NumUnboundDirs is odd
      for( dir=get_NumUnboundDirs( lattice); dir < get_NumVelDirs( lattice); dir+=2)
      {
        pdf_bound_DtH_odd_2( lattice, f_mem_d, cumul_stride
            , subs, dir);
      }
      for( dir=get_NumUnboundDirs( lattice)+1; dir < get_NumVelDirs( lattice); dir+=2)
      {
        pdf_bound_DtH_even_2( lattice, f_mem_d, cumul_stride
            , subs, dir);
      }
#endif
    }
#endif  // BOUNDARY KERNEL
    //cudaEventRecord(stop, 0);
    //cudaEventSynchronize(stop);
    //cudaEventElapsedTime(&timertime,start,stop);
    //totaltime += timertime;
#endif  // #ifdef __CUDACC__

    // Buffers now large enough to deal with all substances in one transfer.
    // Currently, send a receives are activated if __CUDACC__ not defined,
    // but that's bad because because for that case they're totally broken.
    process_send_recv_begin( lattice, 0);
    process_send_recv_end(lattice, 0);

#ifdef __CUDACC__
    //cudaEventRecord(start, 0);
#if BOUNDARY_KERNEL
#if !(POINTER_MAPPING)
    cudaMemcpy( pos_dir_recv_ptr_d
        , lattice->process.pos_dir_pdf_to_recv
        , get_NumBoundDirs( lattice) * get_EndBoundSize( lattice)
        * get_NumSubs( lattice) * sizeof(real) / 2
        , cudaMemcpyHostToDevice);
    checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
    cudaMemcpy( neg_dir_recv_ptr_d
        , lattice->process.neg_dir_pdf_to_recv
        , get_NumBoundDirs( lattice) * get_EndBoundSize( lattice)
        * get_NumSubs( lattice) * sizeof(real) / 2
        , cudaMemcpyHostToDevice);
    checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
#endif  // !(POINTER_MAPPING)

    k_bound_HtD_2
      <<<
      gridboundDim
      , blockboundDim
      >>>( f_mem_d, pos_dir_recv_ptr_d, neg_dir_recv_ptr_d);

    cudaThreadSynchronize();
    checkCUDAError( __FILE__, __LINE__, "boundary kernel");

#else   // !(BOUNDARY_KERNEL)
    for( subs = 0; subs < get_NumSubs( lattice); subs++)
    {
#if 0
      for( dir=get_NumUnboundDirs( lattice); dir < get_NumVelDirs( lattice); dir++)
      {
        pdf_boundary_parallel( lattice, f_mem_d, cumul_stride
            , subs, dir, time, HOST_TO_DEVICE);
      }
#else
      // For D2Q9 and D3Q19, NumUnboundDirs is odd
      for( dir=get_NumUnboundDirs( lattice); dir < get_NumVelDirs( lattice); dir+=2)
      {
        pdf_bound_HtD_odd_2( lattice, f_mem_d, cumul_stride
            , subs, dir);
      }
      for( dir=get_NumUnboundDirs( lattice)+1; dir < get_NumVelDirs( lattice); dir+=2)
      {
        pdf_bound_HtD_even_2( lattice, f_mem_d, cumul_stride
            , subs, dir);
      }
#endif
    }
#endif  // BOUNDARY_KERNEL
    //cudaEventRecord(stop, 0);
    //cudaEventSynchronize(stop);
    //cudaEventElapsedTime(&timertime,start,stop);
    //totaltime += timertime;

#if 0
    // Only implemented in 2D right now
    if( get_NumDims( lattice) == 2)
    {
      // Implement boundary conditions
      bcs_2( lattice, f_mem_d, solids_mem_d); 
    }
#endif

#endif  // #ifdef __CUDACC__
#else   // if !(PARALLEL)
    // Implement this inside get_f1_d functions,
    // as per the old method
#if 0   

#ifdef __CUDACC__
    for( subs = 0; subs < get_NumSubs( lattice); subs++)
    {
      for( dir=get_NumUnboundDirs( lattice); dir < get_NumVelDirs( lattice); dir++)
      {
        pdf_boundary_swap( lattice, f_mem_d, cumul_stride
            , subs, dir, time);
      }
    }
#endif
#endif
#endif  // if PARALLEL

#if 1
    // Only implemented in 2D right now
    if( get_NumDims( lattice) == 2)
    {
      // Implement boundary conditions
      bcs_2( lattice, f_mem_d, solids_mem_d, ns_mem_d); 
    }
#endif

#ifdef __CUDACC__
    if( is_end_of_frame(lattice,time))
    {
      temp_frame_switch = 1;
      cudaMemcpyToSymbol( is_end_of_frame_mem_c, &temp_frame_switch, sizeof(int));
    }


    PEST_gpu_switch( &lattice, mv_mem_d, argc, argv);

    k_collide
      <<<
      gridDim
      , blockDim
      , sizeof(real)*( get_NumVelDirs( lattice)
          + get_NumDims( lattice)+1 )
      * get_BX( lattice)
      * get_BY( lattice)
      * get_BZ( lattice)
      >>>( f_mem_d, mv_mem_d, solids_mem_d, ns_mem_d);
    cudaThreadSynchronize();
    checkCUDAError( __FILE__, __LINE__, "k_collide");

    write_PEST_out_data( &lattice, mv_mem_d, argc, argv);

    if( temp_frame_switch)
    {
      // Although it is necessary to copy device arrays to host at this point,
      // it is inefficient to write them to file here. Rather, we write to file
      // immediately after the first kernel, to allow concurrent host and
      // device execution.  On the last frame of course, it is necessary to
      // write to file here.
      set_frame( lattice, ++frame);

      for( subs = 0; subs < get_NumSubs( lattice); subs++)
      {
        printf("Transferring subs %d "
            "macrovars from device to host for output. \n", subs);
        cudaMemcpy( get_rho_ptr(lattice, subs, 0)
            , mv_mem_d
            + subs*get_NumNodes( lattice)*( 1 + get_NumDims( lattice))
            + 0*get_NumNodes( lattice)
            , get_NumNodes( lattice)*sizeof(real)
            , cudaMemcpyDeviceToHost);
        checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        cudaMemcpy( get_ux_ptr(lattice, subs, 0)
            , mv_mem_d
            + subs*get_NumNodes( lattice)*( 1 + get_NumDims( lattice))
            + 1*get_NumNodes( lattice)
            , get_NumNodes( lattice)*sizeof(real)
            , cudaMemcpyDeviceToHost);
        checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        cudaMemcpy( get_uy_ptr(lattice, subs, 0)
            , mv_mem_d
            + subs*get_NumNodes( lattice)*( 1 + get_NumDims( lattice))
            + 2*get_NumNodes( lattice)
            , get_NumNodes( lattice)*sizeof(real)
            , cudaMemcpyDeviceToHost);
        checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        if(get_NumDims( lattice) == 3)
        {
          cudaMemcpy( get_uz_ptr(lattice, subs, 0)
              , mv_mem_d
              + subs*get_NumNodes( lattice)*( 1 + get_NumDims( lattice))
              + 3*get_NumNodes( lattice)
              , get_NumNodes( lattice)*sizeof(real)
              , cudaMemcpyDeviceToHost);
          checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        }
      }



      output_frame( lattice);
      process_toc( lattice);
      display_etime( lattice);

      temp_frame_switch = 0;
      cudaMemcpyToSymbol( is_end_of_frame_mem_c, &temp_frame_switch, sizeof(int));
    }


#else   // ifndef __CUDACC__
    collide( lattice);
#endif  // ifdef __CUDACC__

#if 0
    if( is_end_of_frame(lattice,time))
    { 
      // Although it is necessary to copy device arrays to host at this point,
      // it is inefficient to write them to file here. Rather, we write to file
      // immediately after the first kernel, to allow concurrent host and
      // device execution.  On the last frame of course, it is necessary to
      // write to file here.
      set_frame( lattice, ++frame);

#ifdef __CUDACC__
      for( subs = 0; subs < get_NumSubs( lattice); subs++)
      {
        printf("Transferring subs %d "
            "macrovars from device to host for output. \n", subs);
        cudaMemcpy( get_rho_ptr(lattice, subs, 0)
            , mv_mem_d
            + subs*get_NumNodes( lattice)*( 1 + get_NumDims( lattice))
            + 0*get_NumNodes( lattice)
            , get_NumNodes( lattice)*sizeof(real)
            , cudaMemcpyDeviceToHost);
        checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        cudaMemcpy( get_ux_ptr(lattice, subs, 0)
            , mv_mem_d
            + subs*get_NumNodes( lattice)*( 1 + get_NumDims( lattice))
            + 1*get_NumNodes( lattice)
            , get_NumNodes( lattice)*sizeof(real)
            , cudaMemcpyDeviceToHost);
        checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        cudaMemcpy( get_uy_ptr(lattice, subs, 0)
            , mv_mem_d
            + subs*get_NumNodes( lattice)*( 1 + get_NumDims( lattice))
            + 2*get_NumNodes( lattice)
            , get_NumNodes( lattice)*sizeof(real)
            , cudaMemcpyDeviceToHost);
        checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        if(get_NumDims( lattice) == 3)
        {
          cudaMemcpy( get_uz_ptr(lattice, subs, 0)
              , mv_mem_d
              + subs*get_NumNodes( lattice)*( 1 + get_NumDims( lattice))
              + 3*get_NumNodes( lattice)
              , get_NumNodes( lattice)*sizeof(real)
              , cudaMemcpyDeviceToHost);
          checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

        }
        if( do_write_f_txt_file( lattice))
        {  
          // In the final version, the CPU arrays may well have
          // the same 'additional boundary' structure as the GPU
          // arrays, in which case a single memcpy can be performed.
          cudaMemcpy( get_fptr(lattice, subs, 0,0,0, 0,0,0, 0)
              , f_mem_d + subs * cumul_stride[get_NumVelDirs( lattice)]
              , get_NumNodes( lattice)
              *(get_NumUnboundDirs( lattice))
              *sizeof(real)
              , cudaMemcpyDeviceToHost);

          checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");

          for( dir=get_NumUnboundDirs( lattice); dir < get_NumVelDirs( lattice); dir+=2)
          {

            cudaMemcpy( get_fptr(lattice, subs, 0,0,0, 0,0,0, dir)
                , f_mem_d + cumul_stride[dir]
                + subs * cumul_stride[get_NumVelDirs( lattice)]
                , 2*get_NumNodes( lattice)
                *sizeof(real)
                , cudaMemcpyDeviceToHost);

            checkCUDAError( __FILE__, __LINE__, "cudaMemcpy");
          }
        }
      }
#endif  // ifdef __CUDACC__
      if( frame == get_NumFrames( lattice))
      {
        output_frame( lattice);
      }
      process_toc( lattice);
      display_etime( lattice);

#ifdef __CUDACC__
      printf("**************\n");
      printf("  __CUDACC__\n");
      printf("**************\n");
#endif
    }
#endif

  } /* for( time=1; time<=lattice->NumTimeSteps; time++) */

  write_PEST_out_file( &lattice, argc, argv);

  // Explicit GPU thread exit call
  cudaThreadExit();

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

  //printf(" \n\n\n Time taken for loop is %f \n\n", totaltime);
  //return 0;

} /* int main( int argc, char **argv) */
