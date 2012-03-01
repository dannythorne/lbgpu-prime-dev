extern __shared__ real fptr[];

__global__
void k_collide(
    real* f_mem_d
    , real* mv_mem_d
    , unsigned char* solids_mem_d
    , real* ns_mem_d
    )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = threadIdx.y + blockIdx.y*blockDim.y;

  int a, subs, k, n, b;

#if __CUDA_ARCH__ < 200
  int klc;
  for( klc=0; klc < kloop_c; klc++)
  {
    k = threadIdx.z + klc*blockDim.z;
#else
    k = threadIdx.z + blockIdx.z*blockDim.z;
#endif

    n = i + j * ni_c + k * nixnj_c;

#if !(COMPUTE_ON_SOLIDS)
    if( d_is_not_solid( solids_mem_d, n + end_bound_c))
    {
#endif

      b = threadIdx.x + threadIdx.y*blockDim.x
        + threadIdx.z*bixbj_c;

      for( subs=0; subs<numsubs_c; subs++)
      {

#if (INAMURO_SIGMA_COMPONENT)
        if( subs != 1 || (time_c >= sigma_t_on_c && time_c <= sigma_t_off_c))
        {
#endif

          // Populate shared memory for a given node with global memory values from
          // that node.  Note that in global memory between stream_collide_stream
          // and collide, each distribution function is swapped with its opposite
          // partner.  The correct ordering is restored while loading into shared
          // memory using the fact that opposite pairs are stored adjacently in v.

          // If fptr = a for all nodes, then the correct values are displayed, therefore
          // fptr working, set_f1d_d working, problem must be with get_f1d_d;

          fptr[b] = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
              , subs
              , i,j,k,n
              , 0,0,0
              , 0, 0);

          for( a=1; a<numdirs_c; a+=2)
          {
            fptr[b + a*blocksize_c]
              = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
                  , subs
                  , i,j,k,n
                  , 0,0,0
                  , a+1, 0);
          }

          for( a=2; a<numdirs_c; a+=2)
          {
            fptr[b + a*blocksize_c]
              = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
                  , subs
                  , i,j,k,n
                  , 0,0,0
                  , a-1, 0);
          }

#if INAMURO_SIGMA_COMPONENT 
          if( subs == 1)
          {
            // Initialize shared memory values for calculating rho.
            fptr[b + (numdirs_c+0)*blocksize_c] = 0.;

            // Note that the velocity macroscopic variables are already
            // set from substance 0, so we don't need to touch them

            // Calculate rho
            for( a=0; a<numdirs_c; a++)
            {
              fptr[b + (numdirs_c+0)*blocksize_c]
                += fptr[b + a*blocksize_c];
            }
#if 0
            if( is_end_of_frame_mem_c)
            {
              for( a=1; a<=numdims_c; a++)
              {
                // Velocity of subs 1 is identical to that of subs 0
                fptr[b + (numdirs_c + a)*blocksize_c] = get_mv_d( mv_mem_d
                    , 0, n, a);
              }
            }
#endif
          }
          else
          {
#endif  // INAMURO_SIGMA_COMPONENT

            // Initialize shared memory values for calculating macro vars.
            fptr[b + (numdirs_c+0)*blocksize_c] = 0.;

            // Calculate macroscopic variables.
            for( a=0; a<numdirs_c; a++)
            {
              fptr[b + (numdirs_c+0)*blocksize_c]
                += fptr[b + a*blocksize_c];

              if( /*debug*/0)
              {
                fptr[b + (numdirs_c+0)*blocksize_c] = 8.;
              }
            }

            if( numdims_c == 2)
            {
              fptr[b + (numdirs_c+1)*blocksize_c] = 0.;
              fptr[b + (numdirs_c+2)*blocksize_c] = 0.;

              if( fptr[b + (numdirs_c+0)*blocksize_c] > EPSILON)
              {
                for( a=0; a<numdirs_c; a++)
                {
                  fptr[b + (numdirs_c+1)*blocksize_c]
                    += vx_c[a]*fptr[b + a*blocksize_c];

                  fptr[b + (numdirs_c+2)*blocksize_c]
                    += vy_c[a]*fptr[b + a*blocksize_c];
                }
#if WALSH_NS_ON
                fptr[b + (numdirs_c+1)*blocksize_c] *=
                  1. - ns_mem_d[n + end_bound_c];

                fptr[b + (numdirs_c+2)*blocksize_c] *=
                  1. - ns_mem_d[n + end_bound_c];
#endif

                fptr[b + (numdirs_c+1)*blocksize_c] /=
                  fptr[b + (numdirs_c+0)*blocksize_c];

                fptr[b + (numdirs_c+2)*blocksize_c] /=
                  fptr[b + (numdirs_c+0)*blocksize_c];
              }   //rho > EPSILON
            }
            else  // numdims_c == 3
            {
              fptr[b + (numdirs_c+1)*blocksize_c] = 0.;
              fptr[b + (numdirs_c+2)*blocksize_c] = 0.;
              fptr[b + (numdirs_c+3)*blocksize_c] = 0.;

              if( fptr[b + (numdirs_c+0)*blocksize_c] > EPSILON)
              {
                for( a=0; a<numdirs_c; a++)
                {
                  fptr[b + (numdirs_c+1)*blocksize_c]
                    += vx_c[a]*fptr[b + a*blocksize_c];

                  fptr[b + (numdirs_c+2)*blocksize_c]
                    += vy_c[a]*fptr[b + a*blocksize_c];

                  fptr[b + (numdirs_c+3)*blocksize_c]
                    += vz_c[a]*fptr[b + a*blocksize_c];

                }
#if WALSH_NS_ON
                fptr[b + (numdirs_c+1)*blocksize_c] *=
                  1. - ns_mem_d[n + end_bound_c];

                fptr[b + (numdirs_c+2)*blocksize_c] *=
                  1. - ns_mem_d[n + end_bound_c];

                fptr[b + (numdirs_c+3)*blocksize_c] *=
                  1. - ns_mem_d[n + end_bound_c];

#endif

                fptr[b + (numdirs_c+1)*blocksize_c] /=
                  fptr[b + (numdirs_c+0)*blocksize_c];

                fptr[b + (numdirs_c+2)*blocksize_c] /=
                  fptr[b + (numdirs_c+0)*blocksize_c];

                fptr[b + (numdirs_c+3)*blocksize_c] /=
                  fptr[b + (numdirs_c+0)*blocksize_c];

              }   //rho > EPSILON
            }     // numdims_c == 2

#if INAMURO_SIGMA_COMPONENT
          }
#endif  // INAMURO_SIGMA_COMPONENT
#if 0
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
#endif
#if 0 //WALSH_NS_ON
          if( d_is_solid( solids_mem_d, n + end_bound_c))
          {
            real temp;
            for( a=1; a<numdirs_c; a+=2)
            { 
              temp = fptr[b + (a+1)*blocksize_c];
            }
            for( a=2; a<numdirs_c; a+=2)
            { 
              fptr[b + a*blocksize_c] = fptr[b + (a-1)*blocksize_c];
            }
            for( a=1; a<numdirs_c; a+=2)
            { 
              fptr[b + a*blocksize_c]=  temp;
            }
          }
          else
          {
#endif

            if( !d_skip_collision_step())
            {
              if( !d_skip_body_force_term())
              {
#if INAMURO_SIGMA_COMPONENT
                if( subs != 1)
                {
#endif
                  // Modify macroscopic variables with a body force
                  for( a=1; a<=numdims_c; a++)
                  {
                    apply_accel_mv( subs, a, b, blocksize_c, n, fptr, ns_mem_d);
                  }
#if INAMURO_SIGMA_COMPONENT
                }
#endif
              }

              // Calculate u-squared since it is used many times
              real usq;
#if INAMURO_SIGMA_COMPONENT
              if( subs != 1)
              {
#endif
                usq = fptr[b + (numdirs_c+1)*blocksize_c]
                  * fptr[b + (numdirs_c+1)*blocksize_c]

                  + fptr[b + (numdirs_c+2)*blocksize_c]
                  * fptr[b + (numdirs_c+2)*blocksize_c];

                if( numdims_c==3)
                {
                  usq += fptr[b + (numdirs_c+3)*blocksize_c]
                    * fptr[b + (numdirs_c+3)*blocksize_c];
                }
#if INAMURO_SIGMA_COMPONENT
              }
#endif

              // Calculate the collision operator and add to f resulting from first
              // streaming

#if WALSH_NS_ON
              real temp1, temp2;
              
              temp1 = fptr[b];
              calc_f_tilde_d( f_mem_d, subs, 0, b, blocksize_c, fptr, usq);
              fptr[b] *= 1. - ns_mem_d[ n + end_bound_c];
              fptr[b] += ns_mem_d[ n + end_bound_c] * temp1;

              for( a=1; a<numdirs_c; a+=2)
              {
                temp1 = fptr[b + a*blocksize_c];
                temp2 = fptr[b + (a+1)*blocksize_c];
                calc_f_tilde_d( f_mem_d, subs, a, b, blocksize_c, fptr, usq);
                calc_f_tilde_d( f_mem_d, subs, a+1, b, blocksize_c, fptr, usq);

                fptr[b + a*blocksize_c] *= 1. - ns_mem_d[ n + end_bound_c];
                fptr[b + a*blocksize_c] += ns_mem_d[ n + end_bound_c] * temp2;
                fptr[b + (a+1)*blocksize_c] *= 1. - ns_mem_d[ n + end_bound_c];
                fptr[b + (a+1)*blocksize_c] += ns_mem_d[ n + end_bound_c] * temp1;

              }

#else
              for( a=0; a<numdirs_c; a++)
              {
                calc_f_tilde_d( f_mem_d, subs, a, b, blocksize_c, fptr, usq);
              }
#endif
            }



#if 0 //WALSH_NS_ON
          } 
#endif


          // Finally, save results back to global memory in the local node.  The
          // ordering was already corrected in the first step, so nothing to worry
          // about here.

          for( a=0; a<numdirs_c; a++)
          {
            set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
                , subs
                , i,j,k,n
                , 0, 0, 0
                , a, fptr[b + a*blocksize_c]);
          }

#if INAMURO_SIGMA_COMPONENT
        }
#endif

      }  /*for( subs=0; subs<numsubs_c; subs++)*/


      for( subs=0; subs<numsubs_c; subs++)
      {


#if (INAMURO_SIGMA_COMPONENT)
        if( subs != 1 || (time_c >= sigma_t_on_c && time_c <= sigma_t_off_c))
        {
#endif

          // Calculate macroscopic variables after the collision step, for the
          // purpose of writing these to host arrays and output files.  Note that
          // for the purpose of efficiency, this (and possibly between sc and s in
          // k_scs) should be the only place in the code at which macroscopic variables
          // are either stored in device global memory or transferred to the host.

          if( is_end_of_frame_mem_c)
          {
            for( a=0; a<numdirs_c; a++)
            {
              fptr[b + a*blocksize_c]
                = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
                    , subs
                    , i,j,k,n
                    , 0,0,0
                    , a, 0);
            }


#if INAMURO_SIGMA_COMPONENT 
            if( subs == 1)
            {
              // Initialize shared memory values for calculating rho.
              fptr[b + (numdirs_c+0)*blocksize_c] = 0.;

              // Note that the velocity macroscopic variables are already
              // set from substance 0, so we don't need to touch them

              // Calculate rho
              for( a=0; a<numdirs_c; a++)
              {
                fptr[b + (numdirs_c+0)*blocksize_c]
                  += fptr[b + a*blocksize_c];
              }
#if 0
              for( a=1; a<=numdims_c; a++)
              {
                // Velocity of subs 1 is identical to that of subs 0
                fptr[b + (numdirs_c + a)*blocksize_c] = get_mv_d( mv_mem_d
                    , 0, n, a);
              }
#endif
            }
            else
            {
#endif  // INAMURO_SIGMA_COMPONENT

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

#if INAMURO_SIGMA_COMPONENT
            }
#endif  // INAMURO_SIGMA_COMPONENT
#if 1
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
                set_mv_d( mv_mem_d, subs, n, a, f_mem_d[cumul_stride_c[1]+n]);
              }
            }
            //}
#endif
          }

#if (INAMURO_SIGMA_COMPONENT)
        }
#endif

      }         // subs loop

#if !(COMPUTE_ON_SOLIDS)
    }  /*if( d_is_not_solid(solids_mem_d, n))*/
#endif


#if __CUDA_ARCH__ < 200
  }  /*for( klc=0; klc < kloop_c; klc++)*/
#endif

}
