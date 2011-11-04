extern __shared__ real fptr[];

__global__
void k_stream_collide_stream(
  real* f_mem_d
, real* mv_mem_d
, real* solids_mem_d )
{
  int n = threadIdx.x + blockIdx.x*blockDim.x;

//If ni and nj are both powers of two, the bitwise operations in the second
//part of the following loop will evaluate more quickly.
//TODO: Figure out how to implement this automatically.
#if 1
  int i = n % ni_c;
  int j = (n % (nixnj_c)) / ni_c;
  int k = n / (nixnj_c);
#else
  int i = n & (ni_c-1);
  int j = (n & (nixnj_c-1)) >> log2(ni_c);
  int k = n >> log2(nixnj_c);
#endif

  int a, subs;

#if 1
  for( subs=0; subs<numsubs_c; subs++)
  {

//Populate shared memory for a given node with global memory values from
//surrounding nodes.  This is a streaming operation.
    for( a=0; a<numdirs_c; a++)
    {
      fptr[threadIdx.x + a*blockDim.x]
        = get_f1d_d( f_mem_d, solids_mem_d
                   , subs, i,j,k, -vx_c[a],-vy_c[a],-vz_c[a], a );
    }

//Initialize shared memory values for calculating macro vars.
    fptr[threadIdx.x + numdirs_c*blockDim.x] = 0.;
    fptr[threadIdx.x + (numdirs_c+1)*blockDim.x] = 0.;
    fptr[threadIdx.x + (numdirs_c+2)*blockDim.x] = 0.;
    if( numdims_c==3)
    {
      fptr[threadIdx.x + (numdirs_c+3)*blockDim.x] = 0.;
    }

//Calculate macroscopic variables.
    for( a=0; a<numdirs_c; a++)
    {
      fptr[threadIdx.x + numdirs_c*blockDim.x]
        += fptr[threadIdx.x + a*blockDim.x];

      if( /*debug*/0)
      {
        fptr[threadIdx.x + numdirs_c*blockDim.x] = 8.;
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

//Store macroscopic variables in global memory.
  for( a=0; a<=numdims_c; a++)
  {
    set_mv_d( mv_mem_d
            , subs, i, j, k, a
            , fptr[threadIdx.x + (numdirs_c + a)*blockDim.x]);
    if( /*debug*/0)
    {
      set_mv_d( mv_mem_d, subs, i, j, k, a, 6.);
    }
  }

#if 1
//Modify macroscopic variables with a body force
  for( a=1; a<=numdims_c; a++)
  {
    apply_accel_mv( subs, a, threadIdx.x, blockDim.x, fptr);
  }
#endif

//Calculate u-squared since it is used many times
    real usq = fptr[threadIdx.x + (numdirs_c+1)*blockDim.x]
             * fptr[threadIdx.x + (numdirs_c+1)*blockDim.x]

             + fptr[threadIdx.x + (numdirs_c+2)*blockDim.x]
             * fptr[threadIdx.x + (numdirs_c+2)*blockDim.x];

    if( numdims_c==3)
    {
      usq += fptr[threadIdx.x + (numdirs_c+3)*blockDim.x]
           * fptr[threadIdx.x + (numdirs_c+3)*blockDim.x];
    }

//Calculate the collision operator and add to f resulting from first streaming
    for( a=0; a<numdirs_c; a++)
    {
      calc_f_tilde_d( f_mem_d, subs, a, threadIdx.x, blockDim.x, fptr, usq);
    }

//Finally, save results back to global memory in adjacent nodes.
//This is a second streaming operation.  Note that the 'swap' that
//occurs in the CPU code is combined with this operation without penalty.
//This utilizes the rearrangement of the v vectors into opposite pairs.
    set_f1d_d( f_mem_d
             , subs, i+vx_c[0], j+vy_c[0], k+vz_c[0], 0
             , fptr[threadIdx.x + 0*blockDim.x]);

    for( a=1; a<numdirs_c; a+=2)
    {
      set_f1d_d( f_mem_d
               , subs, i+vx_c[a], j+vy_c[a], k+vz_c[a], a+1
               , fptr[threadIdx.x + a*blockDim.x]);
    }

    for( a=2; a<numdirs_c; a+=2)
    {
      set_f1d_d( f_mem_d
               , subs, i+vx_c[a], j+vy_c[a], k+vz_c[a], a-1
               , fptr[threadIdx.x + a*blockDim.x]);
    }
//  __syncthreads();
  }
#else
    set_mv_d( mv_mem_d, 0, i, j, k, /*a*/0, 8.);
#endif
}
