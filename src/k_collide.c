// Should be already defined in ./src/system_boundary_kernels.c
//extern __shared__ real fptr[];

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

  for( subs=0; subs<numsubs_c; subs++)
  {
#if (INAMURO_SIGMA_COMPONENT)
    if( subs != 1 || (time_c >= sigma_t_on_c && time_c <= sigma_t_off_c))
    {
#endif

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
            fptr[b + (numdirs_c+0)*blocksize_c]
              *= (1.0 - sink_c * tau_c[1]);

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
                    += ((real) vx_c[a])*fptr[b + a*blocksize_c];

                  fptr[b + (numdirs_c+2)*blocksize_c]
                    += ((real) vy_c[a])*fptr[b + a*blocksize_c];
                }
#if WALSH_NS_ON
                fptr[b + (numdirs_c+1)*blocksize_c] *=
                  (1. - ns_mem_d[n + end_bound_c]);

                fptr[b + (numdirs_c+2)*blocksize_c] *=
                  (1. - ns_mem_d[n + end_bound_c]);
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
                    += ((real) vx_c[a])*fptr[b + a*blocksize_c];

                  fptr[b + (numdirs_c+2)*blocksize_c]
                    += ((real) vy_c[a])*fptr[b + a*blocksize_c];

                  fptr[b + (numdirs_c+3)*blocksize_c]
                    += ((real) vz_c[a])*fptr[b + a*blocksize_c];

                }
#if WALSH_NS_ON
                fptr[b + (numdirs_c+1)*blocksize_c] *=
                  (1. - ns_mem_d[n + end_bound_c]);

                fptr[b + (numdirs_c+2)*blocksize_c] *=
                  (1. - ns_mem_d[n + end_bound_c]);

                fptr[b + (numdirs_c+3)*blocksize_c] *=
                  (1. - ns_mem_d[n + end_bound_c]);

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

          if( !d_skip_body_force_term())
          {
            // Modify macroscopic variables with a body force
#if INAMURO_SIGMA_COMPONENT
            if( subs != 1)
            {
#endif
              for( a=1; a<=numdims_c; a++)
              {
                apply_accel_mv( subs, a, b, blocksize_c, n, fptr, ns_mem_d);
              }
#if INAMURO_SIGMA_COMPONENT
            }
            else
            {
              fptr[b + (numdirs_c+0)*blocksize_c]
                *= (1.0 - sink_c * tau_c[1]);

            }
#endif

          }


#if 1
          if( is_end_of_frame_mem_c)
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
          else
          {
            if( pest_output_flag_c && subs == 1)
            {
              set_mv_d( mv_mem_d
                  , subs, n, 0
                  , fptr[ b + (numdirs_c + 0)*blocksize_c] );
            }
          }

#endif

          if( !d_skip_collision_step())
          {

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


#if !(COMPUTE_ON_SOLIDS)
        }  /*if( d_is_not_solid(solids_mem_d, n))*/
#endif


#if __CUDA_ARCH__ < 200
      }  /*for( klc=0; klc < kloop_c; klc++)*/
#endif

#if (INAMURO_SIGMA_COMPONENT)
    }
#endif

  }         // subs loop


}
