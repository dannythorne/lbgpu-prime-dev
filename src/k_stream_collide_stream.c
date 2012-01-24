extern __shared__ real fptr[];

__global__
void k_stream_collide_stream(
  real* f_mem_d
, real* mv_mem_d
, unsigned char* solids_mem_d )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = threadIdx.y + blockIdx.y*blockDim.y;
  
  int a, subs, k, n, b;
  for( subs=0; subs<numsubs_c; subs++)
  {

    #if __CUDA_ARCH__ < 200
    int klc;
    for( klc=0; klc < kloop_c; klc++)
    {
    k = threadIdx.z + klc*blockDim.z;
    #else
    k = threadIdx.z + blockIdx.z*blockDim.z;
    #endif

    n = i + j * ni_c + k * nixnj_c;
    
//    #if !(COMPUTE_ON_SOLIDS)
//    if( d_is_not_solid( solids_mem_d, n)
//    {
//    #endif
    b = threadIdx.x + threadIdx.y*blockDim.x
        + threadIdx.z*bixbj_c;

    // Populate shared memory for a given node with global memory values from
    // surrounding nodes.  This is a streaming operation. Splitting to odd and
    // even parts is necessary for the boundary condition implementation.

    fptr[b]
      = get_f1d_d( f_mem_d, solids_mem_d
                 , subs
                 , i,j,k,n
                 , -vx_c[0],-vy_c[0],-vz_c[0]
                 , 0, 0);

    for( a=1; a<numdirs_c; a+=2)
    { 
      fptr[b + a*blocksize_c]
        = get_f1d_d( f_mem_d, solids_mem_d
                   , subs
                   , i,j,k,n
                   , -vx_c[a],-vy_c[a],-vz_c[a]
                   , a, 1);
    }
    for( a=2; a<numdirs_c; a+=2)
    { 
      fptr[b + a*blocksize_c]
        = get_f1d_d( f_mem_d, solids_mem_d
                   , subs
                   , i,j,k,n
                   , -vx_c[a],-vy_c[a],-vz_c[a]
                   , a, -1);
    }



#if 1
    // Initialize shared memory values for calculating macro vars.
    fptr[b + (numdirs_c+0)*blocksize_c] = 0.;
    fptr[b + (numdirs_c+1)*blocksize_c] = 0.;
    fptr[b + (numdirs_c+2)*blocksize_c] = 0.;
   if( numdims_c==3)
   {
    fptr[b + (numdirs_c+3)*blocksize_c] = 0.;
   }

    // Calculate macroscopic variables.
    for( a=0; a<numdirs_c; a++)
    {
      fptr[b + (numdirs_c+0)*blocksize_c]
        += fptr[b + a*blocksize_c];

      if( /*debug*/0)
      {
        fptr[b + (numdirs_c+0)*blocksize_c] = 8.;
      }

      fptr[b + (numdirs_c+1)*blocksize_c]
        += vx_c[a]*fptr[b + a*blocksize_c];

      fptr[b + (numdirs_c+2)*blocksize_c]
        += vy_c[a]*fptr[b + a*blocksize_c];

     if( numdims_c==3)
     {
      fptr[b + (numdirs_c+3)*blocksize_c]
        += vz_c[a]*fptr[b + a*blocksize_c];
     }
    }

    fptr[b + (numdirs_c+1)*blocksize_c] /=
      fptr[b + (numdirs_c+0)*blocksize_c];

    fptr[b + (numdirs_c+2)*blocksize_c] /=
      fptr[b + (numdirs_c+0)*blocksize_c];

    if( numdims_c==3)
    {
      fptr[b + (numdirs_c+3)*blocksize_c] /=
        fptr[b + (numdirs_c+0)*blocksize_c];
    }

    if( !d_skip_updating_macrovars())
    {
      // Store macroscopic variables in global memory.
      for( a=0; a<=numdims_c; a++)
      {
        set_mv_d( mv_mem_d
                , subs, n, a
                , fptr[ b + (numdirs_c + a)*blocksize_c] );
        if( /*debug*/0)
        {
          set_mv_d( mv_mem_d, subs, n, a, 6.);
        }
      }
    }

    if( !d_skip_collision_step())
    {
      if( !d_skip_body_force_term())
      {
        // Modify macroscopic variables with a body force
        for( a=1; a<=numdims_c; a++)
        {
          apply_accel_mv( subs, a, b, blocksize_c, fptr);
        }
      }

      // Calculate u-squared since it is used many times
      real usq = fptr[b + (numdirs_c+1)*blocksize_c]
               * fptr[b + (numdirs_c+1)*blocksize_c]

               + fptr[b + (numdirs_c+2)*blocksize_c]
               * fptr[b + (numdirs_c+2)*blocksize_c];

      if( numdims_c==3)
      {
          usq += fptr[b + (numdirs_c+3)*blocksize_c]
               * fptr[b + (numdirs_c+3)*blocksize_c];
      }

      // Calculate the collision operator and add to f resulting from first
      // streaming (TODO: Why do we pass f_mem_d to this function?)
      for( a=0; a<numdirs_c; a++)
      {
        calc_f_tilde_d( f_mem_d, subs, a, b, blocksize_c, fptr, usq);
      }
    }
#else
#if 0
          if( i==1 && j==1)
          {
            fptr[threadIdx.x + C *blockDim.x] = 0.12345678;
            fptr[threadIdx.x + E *blockDim.x] = 1.0;
            fptr[threadIdx.x + W *blockDim.x] = 2.0;
            fptr[threadIdx.x + N *blockDim.x] = 3.0;
            fptr[threadIdx.x + S *blockDim.x] = 4.0;
            fptr[threadIdx.x + NE*blockDim.x] = 5.0;
            fptr[threadIdx.x + SW*blockDim.x] = 6.0;
            fptr[threadIdx.x + NW*blockDim.x] = 7.0;
            fptr[threadIdx.x + SE*blockDim.x] = 8.0;
          }
          else
          {
            fptr[threadIdx.x + C *blockDim.x] = 0.0;
            fptr[threadIdx.x + E *blockDim.x] = 0.0;
            fptr[threadIdx.x + W *blockDim.x] = 0.0;
            fptr[threadIdx.x + N *blockDim.x] = 0.0;
            fptr[threadIdx.x + S *blockDim.x] = 0.0;
            fptr[threadIdx.x + NE*blockDim.x] = 0.0;
            fptr[threadIdx.x + SW*blockDim.x] = 0.0;
            fptr[threadIdx.x + NW*blockDim.x] = 0.0;
            fptr[threadIdx.x + SE*blockDim.x] = 0.0;
          }
#endif
#endif

    // Finally, save results back to global memory in adjacent nodes.  This is
    // a second streaming operation.  Note that the 'swap' that occurs in the
    // CPU code is combined with this operation without penalty.  This utilizes
    // the rearrangement of the v vectors into opposite pairs.
    // Populate shared memory for a given node with global memory values from
    // surrounding nodes.  This is a streaming operation. Splitting to odd and
    // even parts is necessary for the boundary condition implementation.
#if 1
    set_f1d_d( f_mem_d, solids_mem_d
                   , subs
                   , i,j,k,n
                   , vx_c[0],vy_c[0],vz_c[0]
                   , 0, fptr[b + 0*blocksize_c]);

    for( a=1; a<numdirs_c; a+=2)
    {
      set_f1d_d( f_mem_d, solids_mem_d
                   , subs
                   , i,j,k,n
                   , vx_c[a],vy_c[a],vz_c[a]
                   , a+1, fptr[b + a*blocksize_c]);
    }

    for( a=2; a<numdirs_c; a+=2)
    {
      set_f1d_d( f_mem_d, solids_mem_d
                   , subs
                   , i,j,k,n
                   , vx_c[a],vy_c[a],vz_c[a]
                   , a-1, fptr[b + a*blocksize_c]);
    }
#endif
//    #if !(COMPUTE_ON_SOLIDS)
//    }
//    #endif


//  __syncthreads();
  #if __CUDA_ARCH__ < 200
  }  /*for( klc=0; klc < kloop_c; klc++)*/
  #endif
  }  /* subs loop */
}
