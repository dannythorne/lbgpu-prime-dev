__global__
void k_bound_DtH_1(
    real* f_mem_d
    , real* pos_dir_send_ptr_d
    , real* neg_dir_send_ptr_d
    )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = threadIdx.y + blockIdx.y*blockDim.y;

  int n = i + j * ni_c;

  int a, subs;

  for( subs=0; subs<numsubs_c; subs++)
  {
    // This assumes that numunbounddirs_c is odd
    for( a=numunbounddirs_c; a<numdirs_c; a+=2)
    {
      pos_dir_send_ptr_d[subs * numbounddirs_c * end_bound_c / 2  
        + (a - numunbounddirs_c) * end_bound_c / 2 + n]
        = f_mem_d[subs * cumul_stride_c[numdirs_c] 
        + cumul_stride_c[a] - end_bound_c + numnodes_c + n];
    }

    for( a=numunbounddirs_c+1; a<numdirs_c; a+=2)
    {
      neg_dir_send_ptr_d[subs * numbounddirs_c * end_bound_c / 2 
        + (a - numunbounddirs_c-1) * end_bound_c / 2 + n]
        = f_mem_d[subs * cumul_stride_c[numdirs_c] 
        + cumul_stride_c[a] /*+ numnodes_c*/ + n];
    }
  }
}

__global__
void k_bound_DtH_2(
    real* f_mem_d
    , real* pos_dir_send_ptr_d
    , real* neg_dir_send_ptr_d
    )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = threadIdx.y + blockIdx.y*blockDim.y;

  int n = i + j * ni_c;

  int a, subs;

  for( subs=0; subs<numsubs_c; subs++)
  {
    // This assumes that numunbounddirs_c is odd
    for( a=numunbounddirs_c; a<numdirs_c; a+=2)
    {
      neg_dir_send_ptr_d[subs * numbounddirs_c * end_bound_c / 2  
        + (a - numunbounddirs_c) * end_bound_c / 2 + n]
        = f_mem_d[subs * cumul_stride_c[numdirs_c] 
        + cumul_stride_c[a] - end_bound_c + n];
    }

    for( a=numunbounddirs_c+1; a<numdirs_c; a+=2)
    {
      pos_dir_send_ptr_d[subs * numbounddirs_c * end_bound_c / 2  
        + (a - numunbounddirs_c-1) * end_bound_c / 2 + n]
        = f_mem_d[subs * cumul_stride_c[numdirs_c] 
        + cumul_stride_c[a] + numnodes_c + n];
    }
  }
}

__global__
void k_bound_HtD_1(
    real* f_mem_d
    , real* pos_dir_recv_ptr_d
    , real* neg_dir_recv_ptr_d
    )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = threadIdx.y + blockIdx.y*blockDim.y;

  int n = i + j * ni_c;

  int a, subs;

  for( subs=0; subs<numsubs_c; subs++)
  {
    // This assumes that numunbounddirs_c is odd
    for( a=numunbounddirs_c; a<numdirs_c; a+=2)
    {
      f_mem_d[subs * cumul_stride_c[numdirs_c] 
        + cumul_stride_c[a] - end_bound_c + n] = 
        pos_dir_recv_ptr_d[subs * numbounddirs_c * end_bound_c / 2  
        + (a - numunbounddirs_c) * end_bound_c / 2 + n];
    }

    for( a=numunbounddirs_c+1; a<numdirs_c; a+=2)
    {
      f_mem_d[subs * cumul_stride_c[numdirs_c] 
        + cumul_stride_c[a] + numnodes_c + n] = 
        neg_dir_recv_ptr_d[subs * numbounddirs_c * end_bound_c / 2  
        + (a - numunbounddirs_c-1) * end_bound_c / 2 + n];
    }
  }
}

__global__
void k_bound_HtD_2(
    real* f_mem_d
    , real* pos_dir_recv_ptr_d
    , real* neg_dir_recv_ptr_d
    )
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = threadIdx.y + blockIdx.y*blockDim.y;

  int n = i + j * ni_c;

  int a, subs;

  for( subs=0; subs<numsubs_c; subs++)
  {
    // This assumes that numunbounddirs_c is odd
    for( a=numunbounddirs_c; a<numdirs_c; a+=2)
    {
      f_mem_d[subs * cumul_stride_c[numdirs_c] 
        + cumul_stride_c[a] - end_bound_c + numnodes_c + n] = 
        neg_dir_recv_ptr_d[subs * numbounddirs_c * end_bound_c / 2  
        + (a - numunbounddirs_c) * end_bound_c / 2 + n];
    }

    for( a=numunbounddirs_c+1; a<numdirs_c; a+=2)
    {
      f_mem_d[subs * cumul_stride_c[numdirs_c] 
        + cumul_stride_c[a] + n] =  
        pos_dir_recv_ptr_d[subs * numbounddirs_c * end_bound_c / 2  
        + (a - numunbounddirs_c-1) * end_bound_c / 2 + n];
    }
  }
}
