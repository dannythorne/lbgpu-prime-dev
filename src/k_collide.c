extern __shared__ real fptr[];

__global__
void k_collide(
  real* f_mem_d
, real* mv_mem_d
, unsigned char* solids_mem_d
, int* is_end_of_frame_mem_d
)
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
    b = threadIdx.x + threadIdx.y*blockDim.y
        + threadIdx.z*bixbj_c;

    // Populate shared memory for a given node with global memory values from
    // that node.  Note that in global memory, each distribution function is
    // swapped with its opposite partner.  The correct ordering will be
    // restored while loading into shared memory using the fact that opposite
    // pairs are stored adjacently in v.  Hopefully this will occur without
    // bank conflicts.  For Compute 1.0/1.1 devices, this is preferable to
    // restoring the order after the calc.  For Compute 1.3 and greater, it
    // probably doesn't matter.

    fptr[b] = get_f1d_d( f_mem_d, solids_mem_d
                                 , subs
                                 , i,j,k,n
                                 , 0,0,0
                                 , 0, 0);

    for( a=1; a<numdirs_c; a+=2)
    {
      fptr[b + a*blockDim.x]
        = get_f1d_d( f_mem_d, solids_mem_d
                   , subs
                   , i,j,k,n
                   , 0,0,0
                   , a,1 );
    }

    for( a=2; a<numdirs_c; a+=2)
    {
      fptr[b + a*blockDim.x]
        = get_f1d_d( f_mem_d, solids_mem_d
                   , subs
                   , i,j,k,n
                   , 0,0,0
                   , a,-1 );
    }

    // Initialize shared memory values for calculating macro vars.
    fptr[b + (numdirs_c+0)*blocksize] = 0.;
    fptr[b + (numdirs_c+1)*blocksize] = 0.;
    fptr[b + (numdirs_c+2)*blocksize] = 0.;
   if( numdims_c==3)
   {
    fptr[b + (numdirs_c+3)*blocksize] = 0.;
   }

    // Calculate macroscopic variables.
    for( a=0; a<numdirs_c; a++)
    {
      fptr[b + (numdirs_c+0)*blocksize]
        += fptr[b + a*blocksize];

      if( /*debug*/0)
      {
        fptr[b + (numdirs_c+0)*blocksize] = 9.;
      }

      fptr[b + (numdirs_c+1)*blocksize]
        += vx_c[a]*fptr[b + a*blocksize];

      fptr[b + (numdirs_c+2)*blocksize]
        += vy_c[a]*fptr[b + a*blocksize];

     if( numdims_c==3)
     {
      fptr[b + (numdirs_c+3)*blocksize]
        += vz_c[a]*fptr[b + a*blocksize];
     }
    }

    fptr[b + (numdirs_c+1)*blocksize] /=
      fptr[b + (numdirs_c+0)*blocksize];

    fptr[b + (numdirs_c+2)*blocksize] /=
      fptr[b + (numdirs_c+0)*blocksize];

   if( numdims_c==3)
   {
    fptr[b + (numdirs_c+3)*blocksize] /=
      fptr[b + (numdirs_c+0)*blocksize];
   }

    if( !d_skip_updating_macrovars())
    {
      // Store macroscopic variables in global memory.
      for( a=0; a<=numdims_c; a++)
      {
        set_mv_d( mv_mem_d
                , subs, i, j, k, a
                , fptr[b + (numdirs_c + a)*blocksize]);

        if( /*debug*/0)
        {
          set_mv_d( mv_mem_d, subs, i, j, k, a, 7.);
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
          apply_accel_mv( subs, a, b, blocksize, fptr);
        }
      }

      // Calculate u-squared since it is used many times
      real usq = fptr[b + (numdirs_c+1)*blocksize]
               * fptr[b + (numdirs_c+1)*blocksize]

               + fptr[b + (numdirs_c+2)*blocksize]
               * fptr[b + (numdirs_c+2)*blocksize];

      if( numdims_c==3)
      {
          usq += fptr[b + (numdirs_c+3)*blocksize]
               * fptr[b + (numdirs_c+3)*blocksize];
      }

      // Calculate the collision operator and add to f resulting from first
      // streaming
      for( a=0; a<numdirs_c; a++)
      {
        calc_f_tilde_d( f_mem_d, subs, a, b, blocksize, fptr, usq);
      }
    }

    // Finally, save results back to global memory in the local node.  The
    // ordering was already corrected in the first step, so nothing to worry
    // about here.
    for( a=0; a<numdirs_c; a++)
    {
      set_f1d_d( f_mem_d, subs, i, j, k, a, fptr[b + a*blocksize]);
    }

    // I do not understand the purpose of the following section

#if 0

    if( *is_end_of_frame_mem_d)
    {
      // Initialize shared memory values for calculating macro vars.
      fptr[b + (numdirs_c+0)*blocksize] = 0.;
      fptr[b + (numdirs_c+1)*blocksize] = 0.;
      fptr[b + (numdirs_c+2)*blocksize] = 0.;
     if( numdims_c==3)
     {
      fptr[b + (numdirs_c+3)*blocksize] = 0.;
     }

      // Calculate macroscopic variables.
      for( a=0; a<numdirs_c; a++)
      {
        fptr[b + (numdirs_c+0)*blocksize]
          += fptr[b + a*blocksize];

        if( /*debug*/0)
        {
          fptr[b + (numdirs_c+0)*blocksize] = 9.;
        }

        fptr[b + (numdirs_c+1)*blocksize]
          += vx_c[a]*fptr[b + a*blocksize];

        fptr[b + (numdirs_c+2)*blocksize]
          += vy_c[a]*fptr[b + a*blocksize];

       if( numdims_c==3)
       {
        fptr[b + (numdirs_c+3)*blocksize]
          += vz_c[a]*fptr[b + a*blocksize];
       }
      }

      fptr[b + (numdirs_c+1)*blocksize] /=
        fptr[b + (numdirs_c+0)*blocksize];

      fptr[b + (numdirs_c+2)*blocksize] /=
        fptr[b + (numdirs_c+0)*blocksize];

     if( numdims_c==3)
     {
      fptr[b + (numdirs_c+3)*blocksize] /=
        fptr[b + (numdirs_c+0)*blocksize];
     }

      if( !d_skip_updating_macrovars())
      {
        // Store macroscopic variables in global memory.
        for( a=0; a<=numdims_c; a++)
        {
          set_mv_d( mv_mem_d
                  , subs, i, j, k, a
                  , fptr[b + (numdirs_c + a)*blocksize]);

          if( /*debug*/0)
          {
            set_mv_d( mv_mem_d, subs, i, j, k, a, 7.);
          }
        }
      }
    }
#endif
#if __CUDA_ARCH__ < 200
  }  /*for( klc=0; klc < kloop_c; klc++)*/
#endif
  }
}
