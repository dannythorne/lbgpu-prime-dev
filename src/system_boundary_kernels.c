
extern __shared__ real fptr[];
#if 1
__global__
void k_sysbound_pressure_n_1(
    real* f_mem_d
    , unsigned char* solids_mem_d
    , real* ns_mem_d
    )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = nj_c - 1;

  int a, k, n, b;

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
      b = threadIdx.x + threadIdx.z*bixbk_c;

      // Load variables into shared memory. It remains to be seen
      // whether this is necessary or not.

      fptr[b] = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[0],-vy_c[0],-vz_c[0]
          , 0, 0);

      for( a=1; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
            , subs_c
            , i,j,k,n
            , -vx_c[a],-vy_c[a],-vz_c[a]
            , a, 1);
      }
      for( a=2; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
            , subs_c
            , i,j,k,n
            , -vx_c[a],-vy_c[a],-vz_c[a]
            , a, -1);
      }

      fptr[b + numdirs_c*bixbk_c] =
        fptr[b + C*bixbk_c]
        + fptr[b + E*bixbk_c]
        + fptr[b + W*bixbk_c]
        + fptr[b + N*bixbk_c]
        + fptr[b + N*bixbk_c]
        + fptr[b + NE*bixbk_c]
        + fptr[b + NE*bixbk_c]
        + fptr[b + NW*bixbk_c] 
        + fptr[b + NW*bixbk_c];

      fptr[b + numdirs_c*bixbk_c] -= fixed_bound_var_c;
#if ACCURATE
      fptr[b + numdirs_c*bixbk_c] += rho_A_c[subs_c];
#endif

      fptr[b + S*bixbk_c] = fptr[b + N*bixbk_c] 
        - (2./3.) * fptr[b + numdirs_c*bixbk_c];

      fptr[b + SW*bixbk_c] = fptr[b + NE*bixbk_c] 
        + (1./2.) * (fptr[b + E*bixbk_c]-fptr[b + W*bixbk_c])
        - (1./6.) * fptr[b + numdirs_c*bixbk_c];

      fptr[b + SE*bixbk_c] = fptr[b + NW*bixbk_c] 
        + (1./2.) * (fptr[b + W*bixbk_c]-fptr[b + E*bixbk_c])
        - (1./6.) * fptr[b + numdirs_c*bixbk_c];

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[S],-vy_c[S],-vz_c[S]
          , S, 0, fptr[b + S*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[SE],-vy_c[SE],-vz_c[SE]
          , SE, 0, fptr[b + SE*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[SW],-vy_c[SW],-vz_c[SW]
          , SW, 0, fptr[b + SW*bixbk_c]);

#if !(COMPUTE_ON_SOLIDS)
    }
#endif

#if __CUDA_ARCH__ < 200
  }  /*for( klc=0; klc < kloop_c; klc++)*/
#endif

}

__global__
void k_sysbound_pressure_n_2(
    real* f_mem_d
    , unsigned char* solids_mem_d
    , real* ns_mem_d
    )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = nj_c - 1;

  int a, k, n, b;

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
      b = threadIdx.x + threadIdx.z*bixbk_c;

      // Load variables into shared memory. It remains to be seen
      // whether this is necessary or not.

      fptr[b + 0*bixbk_c] 
        = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
            , subs_c
            , i,j,k,n
            , 0,0,0
            , 0, 0);

      for( a=1; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] 
          = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
              , subs_c
              , i,j,k,n
              , 0,0,0
              , a+1, 0);
      }

      for( a=2; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] 
          = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
              , subs_c
              , i,j,k,n
              , 0,0,0
              , a-1, 0);
      }

      fptr[b + numdirs_c*bixbk_c] =
        fptr[b + C*bixbk_c]
        + fptr[b + E*bixbk_c]
        + fptr[b + W*bixbk_c]
        + fptr[b + N*bixbk_c]
        + fptr[b + N*bixbk_c]
        + fptr[b + NE*bixbk_c]
        + fptr[b + NE*bixbk_c]
        + fptr[b + NW*bixbk_c] 
        + fptr[b + NW*bixbk_c];

      fptr[b + numdirs_c*bixbk_c] -= fixed_bound_var_c;
#if ACCURATE
      fptr[b + numdirs_c*bixbk_c] += rho_A_c[subs_c];
#endif


      fptr[b + S*bixbk_c] = fptr[b + N*bixbk_c] 
        - (2./3.) * fptr[b + numdirs_c*bixbk_c];

      fptr[b + SW*bixbk_c] = fptr[b + NE*bixbk_c] 
        + (1./2.) * (fptr[b + E*bixbk_c]-fptr[b + W*bixbk_c])
        - (1./6.) * fptr[b + numdirs_c*bixbk_c];

      fptr[b + SE*bixbk_c] = fptr[b + NW*bixbk_c] 
        + (1./2.) * (fptr[b + W*bixbk_c]-fptr[b + E*bixbk_c])
        - (1./6.) * fptr[b + numdirs_c*bixbk_c];

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , 0,0,0
          , N, 0, fptr[b + S*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , 0,0,0
          , NW, 0, fptr[b + SE*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , 0,0,0
          , NE, 0, fptr[b + SW*bixbk_c]);

#if !(COMPUTE_ON_SOLIDS)
    }
#endif


#if __CUDA_ARCH__ < 200
  }  /*for( klc=0; klc < kloop_c; klc++)*/
#endif

}

__global__
void k_sysbound_pressure_s_1(
    real* f_mem_d
    , unsigned char* solids_mem_d
    , real* ns_mem_d
    )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = 0;

  int a, k, n, b;

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
      b = threadIdx.x + threadIdx.z*bixbk_c;

      // Load variables into shared memory. It remains to be seen
      // whether this is necessary or not.

      fptr[b] = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[0],-vy_c[0],-vz_c[0]
          , 0, 0);

      for( a=1; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
            , subs_c
            , i,j,k,n
            , -vx_c[a],-vy_c[a],-vz_c[a]
            , a, 1);
      }
      for( a=2; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
            , subs_c
            , i,j,k,n
            , -vx_c[a],-vy_c[a],-vz_c[a]
            , a, -1);
      }

      fptr[b + numdirs_c*bixbk_c] =
        - fptr[b + C*bixbk_c]
        - fptr[b + E*bixbk_c]
        - fptr[b + W*bixbk_c]
        - fptr[b + S*bixbk_c]
        - fptr[b + S*bixbk_c]
        - fptr[b + SE*bixbk_c]
        - fptr[b + SE*bixbk_c]
        - fptr[b + SW*bixbk_c] 
        - fptr[b + SW*bixbk_c];

      fptr[b + numdirs_c*bixbk_c] += fixed_bound_var_c;
#if ACCURATE
      fptr[b + numdirs_c*bixbk_c] -= rho_A_c[subs_c];
#endif

      fptr[b + N*bixbk_c] = fptr[b + S*bixbk_c] 
        + (2./3.) * fptr[b + numdirs_c*bixbk_c];

      fptr[b + NE*bixbk_c] = fptr[b + SW*bixbk_c] 
        + (1./2.) * (fptr[b + W*bixbk_c]-fptr[b + E*bixbk_c])
        + (1./6.) * fptr[b + numdirs_c*bixbk_c];

      fptr[b + NW*bixbk_c] = fptr[b + SE*bixbk_c] 
        + (1./2.) * (fptr[b + E*bixbk_c]-fptr[b + W*bixbk_c])
        + (1./6.) * fptr[b + numdirs_c*bixbk_c];

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[N],-vy_c[N],-vz_c[N]
          , N, 0, fptr[b + N*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[NW],-vy_c[NW],-vz_c[NW]
          , NW, 0, fptr[b + NW*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[NE],-vy_c[NE],-vz_c[NE]
          , NE, 0, fptr[b + NE*bixbk_c]);

#if !(COMPUTE_ON_SOLIDS)
    }
#endif


#if __CUDA_ARCH__ < 200
  }  /*for( klc=0; klc < kloop_c; klc++)*/
#endif

}

__global__
void k_sysbound_pressure_s_2(
    real* f_mem_d
    , unsigned char* solids_mem_d
    , real* ns_mem_d
    )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = 0;

  int a, k, n, b;

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
      b = threadIdx.x + threadIdx.z*bixbk_c;

      // Load variables into shared memory. It remains to be seen
      // whether this is necessary or not.

      fptr[b + 0*bixbk_c] 
        = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
            , subs_c
            , i,j,k,n
            , 0,0,0
            , 0, 0);

      for( a=1; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] 
          = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
              , subs_c
              , i,j,k,n
              , 0,0,0
              , a+1, 0);
      }

      for( a=2; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] 
          = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
              , subs_c
              , i,j,k,n
              , 0,0,0
              , a-1, 0);
      }

      fptr[b + numdirs_c*bixbk_c] =
        - fptr[b + C*bixbk_c]
        - fptr[b + E*bixbk_c]
        - fptr[b + W*bixbk_c]
        - fptr[b + S*bixbk_c]
        - fptr[b + S*bixbk_c]
        - fptr[b + SE*bixbk_c]
        - fptr[b + SE*bixbk_c]
        - fptr[b + SW*bixbk_c] 
        - fptr[b + SW*bixbk_c];

      fptr[b + numdirs_c*bixbk_c] += fixed_bound_var_c;
#if ACCURATE
      fptr[b + numdirs_c*bixbk_c] -= rho_A_c[subs_c];
#endif


      fptr[b + N*bixbk_c] = fptr[b + S*bixbk_c] 
        + (2./3.) * fptr[b + numdirs_c*bixbk_c];

      fptr[b + NE*bixbk_c] = fptr[b + SW*bixbk_c] 
        + (1./2.) * (fptr[b + W*bixbk_c]-fptr[b + E*bixbk_c])
        + (1./6.) * fptr[b + numdirs_c*bixbk_c];

      fptr[b + NW*bixbk_c] = fptr[b + SE*bixbk_c] 
        + (1./2.) * (fptr[b + E*bixbk_c]-fptr[b + W*bixbk_c])
        + (1./6.) * fptr[b + numdirs_c*bixbk_c];

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , 0,0,0
          , S, 0, fptr[b + N*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , 0,0,0
          , SE, 0, fptr[b + NW*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , 0,0,0
          , SW, 0, fptr[b + NE*bixbk_c]);

#if !(COMPUTE_ON_SOLIDS)
    }
#endif


#if __CUDA_ARCH__ < 200
  }  /*for( klc=0; klc < kloop_c; klc++)*/
#endif

}


__global__
void k_sysbound_zeroconcgrad_n_1(
    real* f_mem_d
    , unsigned char* solids_mem_d
    , real* ns_mem_d
    )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = nj_c - 1;

  int a, k, n, b;

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
      b = threadIdx.x + threadIdx.z*bixbk_c;

      // Load variables into shared memory. It remains to be seen
      // whether this is necessary or not.

      fptr[b] = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[0],-vy_c[0],-vz_c[0]
          , 0, 0);

      for( a=1; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
            , subs_c
            , i,j,k,n
            , -vx_c[a],-vy_c[a],-vz_c[a]
            , a, 1);
      }
      for( a=2; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
            , subs_c
            , i,j,k,n
            , -vx_c[a],-vy_c[a],-vz_c[a]
            , a, -1);
      }

      fptr[b + numdirs_c*bixbk_c] =
        (fptr[b + N*bixbk_c]
         + fptr[b + NE*bixbk_c]
         + fptr[b + NW*bixbk_c])
        / ( wt_c[S]+wt_c[SW]+wt_c[SE]);

      fptr[b + S*bixbk_c] = wt_c[S]*fptr[b + numdirs_c*bixbk_c];
      fptr[b + SW*bixbk_c] = wt_c[SW]*fptr[b + numdirs_c*bixbk_c];
      fptr[b + SE*bixbk_c] = wt_c[SE]*fptr[b + numdirs_c*bixbk_c];

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[S],-vy_c[S],-vz_c[S]
          , S, 0, fptr[b + S*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[SE],-vy_c[SE],-vz_c[SE]
          , SE, 0, fptr[b + SE*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[SW],-vy_c[SW],-vz_c[SW]
          , SW, 0, fptr[b + SW*bixbk_c]);

#if !(COMPUTE_ON_SOLIDS)
    }
#endif

#if __CUDA_ARCH__ < 200
  }  /*for( klc=0; klc < kloop_c; klc++)*/
#endif

}

__global__
void k_sysbound_zeroconcgrad_n_2(
    real* f_mem_d
    , unsigned char* solids_mem_d
    , real* ns_mem_d
    )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = nj_c - 1;

  int a, k, n, b;

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
      b = threadIdx.x + threadIdx.z*bixbk_c;

      // Load variables into shared memory. It remains to be seen
      // whether this is necessary or not.

      fptr[b + 0*bixbk_c] 
        = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
            , subs_c
            , i,j,k,n
            , 0,0,0
            , 0, 0);

      for( a=1; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] 
          = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
              , subs_c
              , i,j,k,n
              , 0,0,0
              , a+1, 0);
      }

      for( a=2; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] 
          = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
              , subs_c
              , i,j,k,n
              , 0,0,0
              , a-1, 0);
      }

      fptr[b + numdirs_c*bixbk_c] =
        (fptr[b + N*bixbk_c]
         + fptr[b + NE*bixbk_c]
         + fptr[b + NW*bixbk_c])
        / ( wt_c[S]+wt_c[SW]+wt_c[SE]);

      fptr[b + S*bixbk_c] = wt_c[S]*fptr[b + numdirs_c*bixbk_c];
      fptr[b + SW*bixbk_c] = wt_c[SW]*fptr[b + numdirs_c*bixbk_c];
      fptr[b + SE*bixbk_c] = wt_c[SE]*fptr[b + numdirs_c*bixbk_c];

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , 0,0,0
          , N, 0, fptr[b + S*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , 0,0,0
          , NW, 0, fptr[b + SE*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , 0,0,0
          , NE, 0, fptr[b + SW*bixbk_c]);

#if !(COMPUTE_ON_SOLIDS)
    }
#endif


#if __CUDA_ARCH__ < 200
  }  /*for( klc=0; klc < kloop_c; klc++)*/
#endif

}

__global__
void k_sysbound_zeroconcgrad_s_1(
    real* f_mem_d
    , unsigned char* solids_mem_d
    , real* ns_mem_d
    )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = 0;

  int a, k, n, b;

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
      b = threadIdx.x + threadIdx.z*bixbk_c;

      // Load variables into shared memory. It remains to be seen
      // whether this is necessary or not.

      fptr[b] = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[0],-vy_c[0],-vz_c[0]
          , 0, 0);

      for( a=1; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
            , subs_c
            , i,j,k,n
            , -vx_c[a],-vy_c[a],-vz_c[a]
            , a, 1);
      }
      for( a=2; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
            , subs_c
            , i,j,k,n
            , -vx_c[a],-vy_c[a],-vz_c[a]
            , a, -1);
      }

      fptr[b + numdirs_c*bixbk_c] =
        (fptr[b + S*bixbk_c]
         + fptr[b + SW*bixbk_c]
         + fptr[b + SE*bixbk_c])
        / ( wt_c[N]+wt_c[NE]+wt_c[NW]);

      fptr[b + N*bixbk_c] = wt_c[N]*fptr[b + numdirs_c*bixbk_c];
      fptr[b + NE*bixbk_c] = wt_c[NE]*fptr[b + numdirs_c*bixbk_c];
      fptr[b + NW*bixbk_c] = wt_c[NW]*fptr[b + numdirs_c*bixbk_c];

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[N],-vy_c[N],-vz_c[N]
          , N, 0, fptr[b + N*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[NW],-vy_c[NW],-vz_c[NW]
          , NW, 0, fptr[b + NW*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , -vx_c[NE],-vy_c[NE],-vz_c[NE]
          , NE, 0, fptr[b + NE*bixbk_c]);

#if !(COMPUTE_ON_SOLIDS)
    }
#endif


#if __CUDA_ARCH__ < 200
  }  /*for( klc=0; klc < kloop_c; klc++)*/
#endif

}

__global__
void k_sysbound_zeroconcgrad_s_2(
    real* f_mem_d
    , unsigned char* solids_mem_d
    , real* ns_mem_d
    )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = 0;

  int a, k, n, b;

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
      b = threadIdx.x + threadIdx.z*bixbk_c;

      // Load variables into shared memory. It remains to be seen
      // whether this is necessary or not.

      fptr[b + 0*bixbk_c] 
        = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
            , subs_c
            , i,j,k,n
            , 0,0,0
            , 0, 0);

      for( a=1; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] 
          = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
              , subs_c
              , i,j,k,n
              , 0,0,0
              , a+1, 0);
      }

      for( a=2; a<numdirs_c; a+=2)
      { 
        fptr[b + a*bixbk_c] 
          = get_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
              , subs_c
              , i,j,k,n
              , 0,0,0
              , a-1, 0);
      }

      fptr[b + numdirs_c*bixbk_c] =
        (fptr[b + S*bixbk_c]
         + fptr[b + SW*bixbk_c]
         + fptr[b + SE*bixbk_c])
        / ( wt_c[N]+wt_c[NE]+wt_c[NW]);

      fptr[b + N*bixbk_c] = wt_c[N]*fptr[b + numdirs_c*bixbk_c];
      fptr[b + NE*bixbk_c] = wt_c[NE]*fptr[b + numdirs_c*bixbk_c];
      fptr[b + NW*bixbk_c] = wt_c[NW]*fptr[b + numdirs_c*bixbk_c];

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , 0,0,0
          , S, 0, fptr[b + N*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , 0,0,0
          , SE, 0, fptr[b + NW*bixbk_c]);

      set_f1d_d( f_mem_d, solids_mem_d, ns_mem_d
          , subs_c
          , i,j,k,n
          , 0,0,0
          , SW, 0, fptr[b + NE*bixbk_c]);

#if !(COMPUTE_ON_SOLIDS)
    }
#endif


#if __CUDA_ARCH__ < 200
  }  /*for( klc=0; klc < kloop_c; klc++)*/
#endif

}

#endif
