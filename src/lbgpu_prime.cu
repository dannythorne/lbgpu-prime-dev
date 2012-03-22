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

#if 0
__constant__ int vx_c[19];     //
__constant__ int vy_c[19];     // Enough space for D3Q19; first 9
__constant__ int vz_c[19];     // components constitute D2Q9 model
__constant__ real wt_c[19];    // 
__constant__ int cumul_stride_c[20]; // For variable stride created by boundary regions

__constant__ int is_end_of_frame_mem_c;
__constant__ int pest_output_flag_c;

__constant__ real tau_c[2];    // Hardcoded for up to 2 fluid components
__constant__ real gaccel_c[6]; // Hardcoded for up to 2 3D fluid components

__constant__ int numsubs_c;
__constant__ int numdims_c;
__constant__ int numdirs_c;
__constant__ int numbounddirs_c;
__constant__ int numunbounddirs_c;
__constant__ int end_bound_c;
__constant__ int proc_id_c;
__constant__ int num_procs_c;
__constant__ int ni_c;
__constant__ int nj_c;
__constant__ int nk_c;
__constant__ int kloop_c;
__constant__ int bixbj_c;
__constant__ int bjxbk_c;
__constant__ int bixbk_c;
__constant__ int nixnj_c;
__constant__ int blocksize_c;
__constant__ int numnodes_c;
__constant__ int subs_c;

__constant__ real fixed_bound_var_c;

__constant__ real sink_c;
__constant__ real pfmul_c;
__constant__ real pfadd_c;

#if INAMURO_SIGMA_COMPONENT
__constant__ int time_c;
__constant__ int sigma_t_on_c;
__constant__ int sigma_t_off_c;
#endif
#endif

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



  read_PEST_in_files( &lattice, argc, argv);

  for( frame = 0, time=1; time<=get_NumTimeSteps( lattice); time++)
  {
    set_time( lattice, time);


#if INAMURO_SIGMA_COMPONENT
    cudaMemcpyToSymbol( time_c, &(lattice->time), sizeof(int));
    checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
#endif

//------------------------------------------------------//    
// Write boundary cells to buffers and send to host     //
//------------------------------------------------------//
#if PARALLEL
#ifdef __CUDACC__
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
    }

#endif  //BOUNDARY_KERNEL
#endif  // #ifdef __CUDACC__

//------------------------------------------------------//    
// Send boundary buffers to neighboring processes       //
//------------------------------------------------------//

    // Buffers now large enough to deal with all substances in one transfer.
    // Currently, send a receives are activated if __CUDACC__ not defined,
    // but that's bad because for that case they're totally broken.
    process_send_recv_begin( lattice, 0);
    process_send_recv_end(lattice, 0);

//------------------------------------------------------//    
// Write boundary cells to buffers and send to device   //
//------------------------------------------------------//
#ifdef __CUDACC__
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
    }

#endif  // BOUNDARY_KERNEL
#endif  // #ifdef __CUDACC__
#endif  // if (PARALLEL)


//------------------------------------------------------//    
// Apply system boundary conditions                     //
//------------------------------------------------------//
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

    if( is_end_of_frame(lattice,time))
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

      if( is_final_frame( lattice, frame))
      {
        printf("\n breaking out of time loop on frame %d\n", frame);
        break;
      }
    }
#else   // ifndef __CUDACC__
    stream_collide_stream( lattice);
#endif  // ifdef __CUDACC__
 
        printf("\n breaking out of time loop failed on frame %d\n", frame);
    set_time( lattice, ++time);

#if INAMURO_SIGMA_COMPONENT
    cudaMemcpyToSymbol( time_c, &(lattice->time), sizeof(int));
    checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
#endif


//------------------------------------------------------//    
// Write boundary cells to buffers and send to host     //
//------------------------------------------------------//
#if PARALLEL
#ifdef __CUDACC__
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
    }
#endif  // BOUNDARY KERNEL
#endif  // #ifdef __CUDACC__


//------------------------------------------------------//    
// Send buffers to neighboring processes                //
//------------------------------------------------------//

    // Buffers now large enough to deal with all substances in one transfer.
    // Currently, send a receives are activated if __CUDACC__ not defined,
    // but that's bad because because for that case they're totally broken.
    process_send_recv_begin( lattice, 0);
    process_send_recv_end(lattice, 0);


//------------------------------------------------------//    
// Write boundary cells to buffers and send to device   //
//------------------------------------------------------//

#ifdef __CUDACC__
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
    }
#endif  // BOUNDARY_KERNEL
#endif  // #ifdef __CUDACC__
#endif  // if PARALLEL

//------------------------------------------------------//    
// Apply system boundary conditions                     //
//------------------------------------------------------//

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

    if( is_end_of_frame(lattice,time))
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



//      output_frame( lattice);
//      process_toc( lattice);
//      display_etime( lattice);

      temp_frame_switch = 0;
      cudaMemcpyToSymbol( is_end_of_frame_mem_c, &temp_frame_switch, sizeof(int));
    }


#else   // ifndef __CUDACC__
    collide( lattice);
#endif  // ifdef __CUDACC__

    if( is_end_of_frame(lattice,time))
    { 


      // Although it is necessary to copy device arrays to host at this point,
      // it is inefficient to write them to file here. Rather, we write to file
      // immediately after the first kernel, to allow concurrent host and
      // device execution.  On the last frame of course, it is necessary to
      // write to file here.

#ifdef __CUDACC__
      for( subs = 0; subs < get_NumSubs( lattice); subs++)
      {

#if 0
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
#endif
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

        printf("\n broken out of time loop\n");
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
