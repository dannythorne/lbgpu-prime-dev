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

#if !(IGNORE_SOLIDS) && !(COMPUTE_ON_SOLIDS)
      if( d_is_not_solid( solids_mem_d, n + end_bound_c))
      {
#endif

        b = threadIdx.x + threadIdx.y*blockDim.x
          + threadIdx.z*bixbj_c;

        // Populate shared memory for a given node with global memory values from
        // that node.  Note that in global memory between stream_collide_stream
        // and collide, each distribution function is swapped with its opposite
        // partner.  The correct ordering is restored while loading into shared
        // memory using the fact that opposite pairs are stored adjacently in v.

        // If fptr = a for all nodes, then the correct values are displayed, therefore
        // fptr working, set_f1d_d working, problem must be with get_f1d_d;
#if 1
        fptr[b] = get_f1d_d( f_mem_d, solids_mem_d
            , subs
            , i,j,k,n
            , 0,0,0
            , 0, 0);

        for( a=1; a<numdirs_c; a+=2)
        {
          fptr[b + a*blocksize_c]
            = get_f1d_d( f_mem_d, solids_mem_d
                , subs
                , i,j,k,n
                , 0,0,0
                , a+1, 0);
        }

        for( a=2; a<numdirs_c; a+=2)
        {
          fptr[b + a*blocksize_c]
            = get_f1d_d( f_mem_d, solids_mem_d
                , subs
                , i,j,k,n
                , 0,0,0
                , a-1, 0);
        }
#endif
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
            fptr[b + (numdirs_c+0)*blocksize_c] = 9.;
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
                , fptr[b + (numdirs_c + a)*blocksize_c]);

            if( /*debug*/0)
            {
              set_mv_d( mv_mem_d, subs, n, a, 7.);
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
          // streaming
          for( a=0; a<numdirs_c; a++)
          {
            calc_f_tilde_d( f_mem_d, subs, a, b, blocksize_c, fptr, usq);
          }
        }
        // Finally, save results back to global memory in the local node.  The
        // ordering was already corrected in the first step, so nothing to worry
        // about here.
#if 1
        for( a=0; a<numdirs_c; a++)
        {
          set_f1d_d( f_mem_d, solids_mem_d
              , subs
              , i,j,k,n
              , 0, 0, 0
              , a, fptr[b + a*blocksize_c]);
        }
#endif
        // Calculate macroscopic variables after the collision step, for the
        // purpose of writing these to host arrays and output files.  Note that
        // for the purpose of efficiency, this (and possibly between sc and s in
        // k_scs) should be the only place in the code at which macroscopic variables
        // are either stored in device global memory or transferred to the host.
#if 1

        if( is_end_of_frame_mem_c)
        {
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
              fptr[b + (numdirs_c+0)*blocksize_c] = 9.;
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

          //if( !d_skip_updating_macrovars())
          // {
          // Store macroscopic variables in global memory.
          for( a=0; a<=numdims_c; a++)
          {
            set_mv_d( mv_mem_d
                , subs, n, a
                , fptr[b + (numdirs_c + a)*blocksize_c]);

            if( /*debug*/0)
            {
              set_mv_d( mv_mem_d, subs, n, a, (real) n);
            }
          }
          //}
        }
#endif
#if !(IGNORE_SOLIDS) && !(COMPUTE_ON_SOLIDS)
      }  /*if( d_is_not_solid(solids_mem_d, n))*/
#endif
#if __CUDA_ARCH__ < 200
    }  /*for( klc=0; klc < kloop_c; klc++)*/
#endif
  }  /*for( subs=0; subs<numsubs_c; subs++)*/
}
