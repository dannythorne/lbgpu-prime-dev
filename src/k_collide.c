extern __shared__ real fptr[];

__global__
void k_collide(
  real* f_mem_d
, real* mv_mem_d
, unsigned char* solids_mem_d )
{
  int n = threadIdx.x + blockIdx.x*blockDim.x;

#if 1
  int i = n % ni_c;
  int j = (n % (nixnj_c)) / ni_c;
  int k = n / (nixnj_c);
#else
  // If ni and nj are both powers of two, these bitwise operations will
  // evaluate more quickly than the above more general arithmetic.
  // TODO: Figure out how to implement this automatically.
  int i = n & (ni_c-1);
  int j = (n & (nixnj_c-1)) >> log2(ni_c);
  int k = n >> log2(nixnj_c);
#endif

  int a, subs;

#if 1
  for( subs=0; subs<numsubs_c; subs++)
  {
    // Populate shared memory for a given node with global memory values from
    // that node.  Note that in global memory, each distribution function is
    // swapped with its opposite partner.  The correct ordering will be
    // restored while loading into shared memory using the fact that opposite
    // pairs are stored adjacently in v.  Hopefully this will occur without
    // bank conflicts.  For Compute 1.0/1.1 devices, this is preferable to
    // restoring the order after the calc.  For Compute 1.3 and greater, it
    // probably doesn't matter.

    fptr[threadIdx.x] = get_f1d_d( f_mem_d, solids_mem_d
                                 , subs
                                 , i,j,k
                                 , 0,0,0
                                 , 0 );

    for( a=1; a<numdirs_c; a+=2)
    {
      fptr[threadIdx.x + a*blockDim.x]
        = get_f1d_d( f_mem_d, solids_mem_d
                   , subs
                   , i,j,k
                   , 0,0,0
                   , a+1 );
    }

    for( a=2; a<numdirs_c; a+=2)
    {
      fptr[threadIdx.x + a*blockDim.x]
        = get_f1d_d( f_mem_d, solids_mem_d
                   , subs
                   , i,j,k
                   , 0,0,0
                   , a-1 );
    }

    // Initialize shared memory values for calculating macro vars.
    fptr[threadIdx.x + (numdirs_c+0)*blockDim.x] = 0.;
    fptr[threadIdx.x + (numdirs_c+1)*blockDim.x] = 0.;
    fptr[threadIdx.x + (numdirs_c+2)*blockDim.x] = 0.;
   if( numdims_c==3)
   {
    fptr[threadIdx.x + (numdirs_c+3)*blockDim.x] = 0.;
   }

    // Calculate macroscopic variables.
    for( a=0; a<numdirs_c; a++)
    {
      fptr[threadIdx.x + (numdirs_c+0)*blockDim.x]
        += fptr[threadIdx.x + a*blockDim.x];

      if( /*debug*/0)
      {
        fptr[threadIdx.x + (numdirs_c+0)*blockDim.x] = 9.;
      }

      fptr[threadIdx.x + (numdirs_c+1)*blockDim.x]
        += vx_c[a]*fptr[threadIdx.x + a*blockDim.x];

      fptr[threadIdx.x + (numdirs_c+2)*blockDim.x]
        += vy_c[a]*fptr[threadIdx.x + a*blockDim.x];

     if( numdims_c==3)
     {
      fptr[threadIdx.x + (numdirs_c+3)*blockDim.x]
        += vz_c[a]*fptr[threadIdx.x + a*blockDim.x];
     }
    }

    fptr[threadIdx.x + (numdirs_c+1)*blockDim.x] /=
      fptr[threadIdx.x + (numdirs_c+0)*blockDim.x];

    fptr[threadIdx.x + (numdirs_c+2)*blockDim.x] /=
      fptr[threadIdx.x + (numdirs_c+0)*blockDim.x];

   if( numdims_c==3)
   {
    fptr[threadIdx.x + (numdirs_c+3)*blockDim.x] /=
      fptr[threadIdx.x + (numdirs_c+0)*blockDim.x];
   }

    if( !d_skip_updating_macrovars())
    {
      // Store macroscopic variables in global memory.
      for( a=0; a<=numdims_c; a++)
      {
        set_mv_d( mv_mem_d
                , subs, i, j, k, a
                , fptr[threadIdx.x + (numdirs_c + a)*blockDim.x]);

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
          apply_accel_mv( subs, a, threadIdx.x, blockDim.x, fptr);
        }
      }

      // Calculate u-squared since it is used many times
      real usq = fptr[threadIdx.x + (numdirs_c+1)*blockDim.x]
               * fptr[threadIdx.x + (numdirs_c+1)*blockDim.x]

               + fptr[threadIdx.x + (numdirs_c+2)*blockDim.x]
               * fptr[threadIdx.x + (numdirs_c+2)*blockDim.x];

      if( numdims_c==3)
      {
          usq += fptr[threadIdx.x + (numdirs_c+3)*blockDim.x]
               * fptr[threadIdx.x + (numdirs_c+3)*blockDim.x];
      }

      // Calculate the collision operator and add to f resulting from first
      // streaming
      for( a=0; a<numdirs_c; a++)
      {
        calc_f_tilde_d( f_mem_d, subs, a, threadIdx.x, blockDim.x, fptr, usq);
      }
    }

    // Finally, save results back to global memory in the local node.  The
    // ordering was already corrected in the first step, so nothing to worry
    // about here.
    for( a=0; a<numdirs_c; a++)
    {
      set_f1d_d( f_mem_d, subs, i, j, k, a, fptr[threadIdx.x + a*blockDim.x]);
    }

#if 0 // TODO: Should recompute macrovars before outputting them.

    // Initialize shared memory values for calculating macro vars.
    fptr[threadIdx.x + (numdirs_c+0)*blockDim.x] = 0.;
    fptr[threadIdx.x + (numdirs_c+1)*blockDim.x] = 0.;
    fptr[threadIdx.x + (numdirs_c+2)*blockDim.x] = 0.;
   if( numdims_c==3)
   {
    fptr[threadIdx.x + (numdirs_c+3)*blockDim.x] = 0.;
   }

    // Calculate macroscopic variables.
    for( a=0; a<numdirs_c; a++)
    {
      fptr[threadIdx.x + (numdirs_c+0)*blockDim.x]
        += fptr[threadIdx.x + a*blockDim.x];

      if( /*debug*/0)
      {
        fptr[threadIdx.x + (numdirs_c+0)*blockDim.x] = 9.;
      }

      fptr[threadIdx.x + (numdirs_c+1)*blockDim.x]
        += vx_c[a]*fptr[threadIdx.x + a*blockDim.x];

      fptr[threadIdx.x + (numdirs_c+2)*blockDim.x]
        += vy_c[a]*fptr[threadIdx.x + a*blockDim.x];

     if( numdims_c==3)
     {
      fptr[threadIdx.x + (numdirs_c+3)*blockDim.x]
        += vz_c[a]*fptr[threadIdx.x + a*blockDim.x];
     }
    }

    fptr[threadIdx.x + (numdirs_c+1)*blockDim.x] /=
      fptr[threadIdx.x + (numdirs_c+0)*blockDim.x];

    fptr[threadIdx.x + (numdirs_c+2)*blockDim.x] /=
      fptr[threadIdx.x + (numdirs_c+0)*blockDim.x];

   if( numdims_c==3)
   {
    fptr[threadIdx.x + (numdirs_c+3)*blockDim.x] /=
      fptr[threadIdx.x + (numdirs_c+0)*blockDim.x];
   }

    if( !d_skip_updating_macrovars())
    {
      // Store macroscopic variables in global memory.
      for( a=0; a<=numdims_c; a++)
      {
        set_mv_d( mv_mem_d
                , subs, i, j, k, a
                , fptr[threadIdx.x + (numdirs_c + a)*blockDim.x]);

        if( /*debug*/0)
        {
          set_mv_d( mv_mem_d, subs, i, j, k, a, 7.);
        }
      }
    }
#endif

//  __syncthreads();
  }
#else
  //mv_mem_d[n] = 9.;
  //set_mv_d( mv_mem_d, 0, i, j, k, /*a*/0, 9.);
#endif
}
