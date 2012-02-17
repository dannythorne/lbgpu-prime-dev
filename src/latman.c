//##############################################################################
//
// all_nodes_latman.c
//
//  - Lattice Manager.
//
//  - Routines for managing a lattice:
//
//    - construct
//    - init
//    - destruct
//
//  - This file, with prefix "all_nodes_", is for the version of the code
//    that stores all nodes of the domain even if they are interior solid
//    nodes that are not involved in any computations.  This is more
//    efficient, in spite of storing unused nodes, if the ratio of
//    interior solid nodes to total nodes is sufficiently low.  How
//    low is sufficient? is a difficult question.  If storage is the
//    only consideration, then the two approaches balance at somewhere
//    around .75 .  But more observations need to be made to characterize
//    the trade-off in terms of computational efficiency.
//

void process_matrix( struct lattice_struct *lattice, int **matrix)
{
  // Variable declarations.
  int i,  j;
  int ei, ej;
  int n;

  // Ending indices.
  ei = get_LX(lattice)-1;
  ej = get_LY(lattice)-1;

  for( j=0; j<=ej; j++)
  {
    n = j*get_LX(lattice);

    for( i=0; i<=ei; i++, n++)
    {
      set_is_solid( lattice, n, matrix[j][i]);

    } /* for( i=0; i<=ei; i++, n++) */

  } /* for( j=0; j<=ej; j++) */

} /* void process_matrix( struct lattice_struct *lattice, int **matrix) */

// void construct_lattice( struct lattice_struct *lattice)
//##############################################################################
//
// C O N S T R U C T   L A T T I C E
//
//  - Construct lattice.
//
void construct_lattice( lattice_ptr *lattice, int argc, char **argv)
{
  // Variable declarations
  int    i,
         j;
  int    subs;
  char   filename[1024];
  char   dirname[1024];

  // Allocate the lattice structure.
  *lattice = ( struct lattice_struct*)malloc( sizeof(struct lattice_struct));

  assert(*lattice!=NULL);

  process_init( *lattice, argc, argv);

  // Read problem parameters
  if( argc == 2)
  {
    sprintf( dirname, "in_%s", argv[1]);
  }
  else
  {
    sprintf( dirname, "in");
  }
  sprintf(filename, "./%s/%s", dirname, "params.in");
  printf("%s %d %04d >> Reading params from file \"%s\"\n",
      __FILE__, __LINE__, get_proc_id( *lattice), filename);

  printf("Before read_params: rho_A = %f\n",get_rho_A(*lattice,0));

  read_params( *lattice, filename);

  printf("After read_params: rho_A = %f\n",get_rho_A(*lattice,0));
  printf("After read_params: rho_B = %f\n",get_rho_B(*lattice,0));

  set_NumSubs( *lattice, NUM_FLUID_COMPONENTS);
  set_NumDims( *lattice, NUM_DIMENSIONS);
  set_NumVelDirs( *lattice, 2*NUM_DIMENSIONS*NUM_DIMENSIONS + 1);

  if( get_NumDims( *lattice) == 2)
  {
    set_NumBoundDirs( *lattice, 6);
  }

  if( get_NumDims( *lattice) == 3)
  {
    set_NumBoundDirs( *lattice, 10);
  }

  set_NumUnboundDirs( *lattice, get_NumVelDirs( *lattice)
      - get_NumBoundDirs( *lattice));

  process_compute_local_params( *lattice);


#ifdef __CUDACC__
  // Copy certain parameters to __constant__ memory on the gpu device

  cudaMemcpyToSymbol( vx_c, vx, 19*sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");

  cudaMemcpyToSymbol( vy_c, vy, 19*sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");

  cudaMemcpyToSymbol( vz_c, vz, 19*sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");

  int a;

  if( get_NumDims( *lattice)==2)
  {
    real W0 = 4./9.;
    real W1 = 1./9.;
    real W2 = 1./36.;
    // wt = { W0
    //      , W1, W1, W1, W1
    //      , W2, W2, W2, W2
    //      , 0., 0.
    //      , 0., 0., 0., 0., 0., 0., 0., 0.};
    wt[0] = W0;
    for( a=1; a<5 ; a++) { wt[a] = W1;}
    for(    ; a<9 ; a++) { wt[a] = W2;}
    for(    ; a<19; a++) { wt[a] = 0.;}

  }

  if( get_NumDims( *lattice) == 3)
  {
    real W0 = 1./3.;
    real W1 = 1./18.;
    real W2 = 1./36.;
    // wt = { W0
    //      , W1, W1, W1, W1
    //      , W2, W2, W2, W2
    //      , W1, W1
    //      , W2, W2, W2, W2, W2, W2, W2, W2};
    wt[0] = W0;
    for( a=1; a<5 ; a++) { wt[a] = W1;}
    for(    ; a<9 ; a++) { wt[a] = W2;}
    for(    ; a<11; a++) { wt[a] = W1;}
    for(    ; a<19; a++) { wt[a] = W2;}
  }

  // For direction a, cumul_stride[a] is the 'cumulative stride',
  // i.e. the first non-boundary memory address in direction a is 
  // f + cumul_stride[a].  A final (additional) entry is the size
  // of a single substance's worth of pdfs.

  cumul_stride[0] = 0;

  for( a=1; a<get_NumUnboundDirs( *lattice); a++)
  {
    cumul_stride[a] = cumul_stride[a-1] + get_NumNodes( *lattice);
  }

  cumul_stride[a++] = cumul_stride[a-1] 
    + get_NumNodes( *lattice) + get_EndBoundSize( *lattice);
  cumul_stride[a++] = cumul_stride[a-1] 
    + get_NumNodes( *lattice);

  for(    ; a<get_NumVelDirs( *lattice); a+=2) 
  { 
    cumul_stride[a] = cumul_stride[a-2] + 2 * get_NumNodes( *lattice) 
      + 2 * get_EndBoundSize( *lattice);
    cumul_stride[a+1] = cumul_stride[a-1] + 2 * get_NumNodes( *lattice) 
      + 2 * get_EndBoundSize( *lattice);
  }

  cumul_stride[a] = cumul_stride[a-1] + get_NumNodes( *lattice) 
    + get_EndBoundSize( *lattice);
  cudaMemcpyToSymbol( wt_c, wt, 19*sizeof(real));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  cudaMemcpyToSymbol( cumul_stride_c, cumul_stride, 20*sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");

  cudaMemcpyToSymbol( tau_c
      , get_tau_ptr( *lattice)
      , get_NumSubs(*lattice)*sizeof(real) );
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");

  /*  printf(" \n\n %f %f %f %f %f %f \n\n", *((*lattice)->param.gforce[0])
      , *((*lattice)->param.gforce[0]+1)
      , *((*lattice)->param.gforce[0]+2)
      , *((*lattice)->param.gforce[0]+3)
      , *((*lattice)->param.gforce[0]+4)
      , *((*lattice)->param.gforce[0]+5));*/




  for( subs = 0; subs < get_NumSubs(*lattice); subs++)
  {
    cudaMemcpyToSymbol(
        gaccel_c
        , get_gaccel_ptr( *lattice, subs)
        , get_NumDims(*lattice)*sizeof(real)
        , subs*3*sizeof(real) // offset from gaccel_c
        , cudaMemcpyHostToDevice);
    checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  }

  // cudaMemcpyToSymbol( gaccel_c
  //                   , get_gaccel_ptr( *lattice, 0)
  //                   , get_NumDims(*lattice)*sizeof(real));
  // if(get_NumSubs(*lattice) == 2)
  // {
  //   cudaMemcpyToSymbol( gaccel_c
  //                     , get_gaccel_ptr( *lattice, 1)
  //                     , get_NumDims(*lattice)*sizeof(real)
  //                     , get_NumDims(*lattice)*sizeof(real)
  //                     , cudaMemcpyHostToDevice);
  // }

  cudaMemcpyToSymbol( numsubs_c, get_NumSubs_ptr( *lattice), sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  cudaMemcpyToSymbol( numdims_c, get_NumDims_ptr( *lattice), sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  cudaMemcpyToSymbol( numdirs_c, get_NumVelDirs_ptr( *lattice), sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  cudaMemcpyToSymbol( numbounddirs_c, get_NumBoundDirs_ptr( *lattice), sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  cudaMemcpyToSymbol( numunbounddirs_c, get_NumUnboundDirs_ptr( *lattice), sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  cudaMemcpyToSymbol( end_bound_c, get_EndBoundSize_ptr( *lattice), sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  cudaMemcpyToSymbol( proc_id_c, get_proc_id_ptr( *lattice), sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  cudaMemcpyToSymbol( num_procs_c, get_num_procs_ptr( *lattice), sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");

  cudaMemcpyToSymbol( ni_c, get_ni_ptr( *lattice), sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  cudaMemcpyToSymbol( nj_c, get_nj_ptr( *lattice), sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  cudaMemcpyToSymbol( nk_c, get_nk_ptr( *lattice), sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");

  int temp;
#ifdef __CUDACC__
#if __CUDACC__ >= 200
  temp = 1;
  cudaMemcpyToSymbol( kloop_c, &temp, sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
#else
  temp = get_LZ( *lattice) / get_BZ( *lattice);
  cudaMemcpyToSymbol( kloop_c, &temp, sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
#endif
#endif

  temp = get_BX( *lattice) * get_BY( *lattice);
  cudaMemcpyToSymbol( bixbj_c, &temp, sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  temp = get_BX( *lattice) * get_BZ( *lattice);
  cudaMemcpyToSymbol( bixbk_c, &temp, sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  temp = get_BY( *lattice) * get_BZ( *lattice);
  cudaMemcpyToSymbol( bjxbk_c, &temp, sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  temp = get_ni( *lattice) * get_nj( *lattice);
  cudaMemcpyToSymbol( nixnj_c, &temp, sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  temp = get_BX( *lattice)*get_BY( *lattice)*get_BZ( *lattice);
  cudaMemcpyToSymbol( blocksize_c, &temp, sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
  cudaMemcpyToSymbol( numnodes_c, get_NumNodes_ptr( *lattice), sizeof(int));
  checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
#endif

  //The following has been moved above the cudaMemcpyToSymbol section
  //process_compute_local_params( *lattice);

  // Allocate space for solids.
  (*lattice)->solids_memblock =
    (unsigned char*)malloc((*lattice)->NumNodes*sizeof(unsigned char));

  if( get_NumDims(*lattice)==2)
  {
    // Allocate matrix for storing information from bmp file.
    int** matrix = (int**)malloc( get_LY(*lattice)*sizeof(int*));
    for( j=0; j<get_LY(*lattice); j++)
    {
      matrix[j] = (int*)malloc( get_LX(*lattice)*sizeof(int));
    }

    // Initialize matrix[][].
    for( j=0; j<get_LY(*lattice); j++)
    {
      for( i=0; i<get_LX(*lattice); i++)
      {
        matrix[j][i] = 0;
      }
    }

    // Read solids from .bmp file.
    sprintf(filename, "./in/%dx%d.bmp",
        get_g_LX(*lattice),
        get_g_LY(*lattice) );
    printf("%s %d >> Reading solids from file \"%s\"\n",
        __FILE__, __LINE__, filename);
    spy_bmp( *lattice, filename, matrix);
    process_matrix( *lattice, matrix);
  }
  else
  {
    // Read solids from .raw file.
    sprintf(filename, "./in/%dx%dx%d.raw",
        get_g_LX(*lattice),
        get_g_LY(*lattice),
        get_g_LZ(*lattice) );
    printf("%s %d >> Reading solids from file \"%s\"\n",
        __FILE__, __LINE__, filename);
    read_solids( *lattice, filename);
  }

  // Allocate vars struct.
  (*lattice)->vars =
    ( struct vars_struct*)malloc(
        get_NumSubs(*lattice)*sizeof(struct vars_struct));

  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
    // Allocate NumNodes particle distribution functions.
    (*lattice)->vars[subs].f_memblock =
      (real*)malloc(
          sizeof(real)
          * get_NumNodes(*lattice)
          * get_NumVelDirs(*lattice)
          );
    for( i=0; i<get_NumNodes(*lattice)* get_NumVelDirs(*lattice); i++)
    {
      (*lattice)->vars[subs].f_memblock[i] = 0.;
    }
    if( (*lattice)->vars[subs].f_memblock==NULL)
    {
      printf(
          "%s %d %04d >> "
          "construct_lattice() -- ERROR:  "
          "Attempt to allocate %d struct pdf_struct types failed.  "
          "Exiting!\n",
          __FILE__,__LINE__,get_proc_id(*lattice),
          (*lattice)->NumNodes
          );
      process_exit(1);
    }

#ifdef __CUDACC__
    // Allocate mem block for f on the gpu device
    //printf(" Allocating memory on the device \n");

    f_mem_size = get_NumSubs( *lattice) * cumul_stride[get_NumVelDirs( *lattice)];
    // printf(" \n\n\n XXXXXXXXX %d \n\n\n\n\n", f_mem_size);
    cudaMalloc( (void**)&f_mem_d, f_mem_size*sizeof(real));
    checkCUDAError( __FILE__, __LINE__, "cudaMalloc");

    // Allocate mem block for macro vars on the gpu device
    mv_mem_size = get_NumSubs( *lattice)
      * get_NumNodes( *lattice)
      * (1 + get_NumDims( *lattice));

    cudaMalloc( (void**)&mv_mem_d, mv_mem_size*sizeof(real));
    checkCUDAError( __FILE__, __LINE__, "cudaMalloc");

    solids_mem_size = get_NumNodes( *lattice)
      + 2 * get_EndBoundSize( *lattice);

    cudaMalloc( (void**)&solids_mem_d, solids_mem_size*sizeof(unsigned char));
    checkCUDAError( __FILE__, __LINE__, "cudaMalloc");
#if 0
    cudaMalloc( (void**)&is_end_of_frame_mem_d, sizeof(int));
    checkCUDAError( __FILE__, __LINE__, "cudaMalloc");
#endif
#endif

    // Allocate pointers to the individual f arrays.
    (*lattice)->vars[subs].f1d =
      (real**)malloc( sizeof(real*)*get_NumVelDirs(*lattice));
    for( i=0; i<get_NumVelDirs(*lattice); i++)
    {
      (*lattice)->vars[subs].f1d[i] =
        (*lattice)->vars[subs].f_memblock + i*get_NumNodes(*lattice);
    }

    // Allocate NumNodes macroscopic variables.
    printf("numDims = %d\n",get_NumDims(*lattice));
    printf("numNodes = %d\n",get_NumNodes(*lattice));
    (*lattice)->vars[subs].macrovars_memblock =
      (real*)malloc(
          sizeof(real)*get_NumNodes(*lattice)*( get_NumDims(*lattice)+1));
    if( (*lattice)->vars[subs].macrovars_memblock==NULL)
    {
      printf(
          "%s %d %04d >> "
          "construct_lattice() -- ERROR:  "
          "Attempt to allocate %d real types failed.  "
          "Exiting!\n",
          __FILE__,__LINE__,get_proc_id(*lattice),
          (*lattice)->NumNodes
          );
      process_exit(1);
    }

    // Allocate pointers to the individual macrovar arrays.
    (*lattice)->vars[subs].macrovars1d =
      (real**)malloc( sizeof(real*)*( get_NumDims(*lattice)+1));
    for( i=0; i<(get_NumDims(*lattice)+1); i++)
    {
      (*lattice)->vars[subs].macrovars1d[i] =
        (*lattice)->vars[subs].macrovars_memblock + i*get_NumNodes(*lattice);
    }

#if 0 // TODO: Revise for new data structures or delete.
#if NON_LOCAL_FORCES
    // Allocate NumNodes elements for force.
    (*lattice)->force[subs] =
      ( struct force_struct*)malloc(
          (*lattice)->NumNodes*sizeof( struct force_struct));
    if( (*lattice)->force[subs]==NULL)
    {
      printf(
          "construct_lattice() -- ERROR:  "
          "Attempt to allocate %d struct force_struct types failed.  "
          "Exiting!\n",
          (*lattice)->NumNodes
          );
      process_exit(1);
    }
#endif /* NON_LOCAL_FORCES */
#endif

  } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

#if 0 // TODO: Revise for new data structures or delete.
#if STORE_UEQ
  // Allocate NumNodes elements for ueq.
  (*lattice)->ueq =
    ( struct ueq_struct*)malloc(
        (*lattice)->NumNodes*sizeof( struct ueq_struct));
  if( (*lattice)->ueq==NULL)
  {
    printf(
        "construct_lattice() -- ERROR:  "
        "Attempt to allocate %d struct ueq_struct types failed.  "
        "Exiting!\n",
        (*lattice)->NumNodes
        );
    process_exit(1);
  }
#endif /* STORE_UEQ */

#if POROUS_MEDIA
  switch( (*lattice)->param.ns_flag)
  {
    case 0:
      {
        if( (*lattice)->param.ns > 1. || (*lattice)->param.ns < 0.)
        {
          printf(
              "latman.c: construct_lattice() -- "
              "ERROR: ns = %f. "
              "Should have 0 <= ns <=1. "
              "Exiting!\n", (*lattice)->param.ns
              );
          process_exit(1);
        }
        break;
      }

    case 1:
      {
        // Allocate space for ns values.
        (*lattice)->ns =
          (struct ns_struct*)
          malloc((*lattice)->NumNodes*sizeof(struct ns_struct));
        if( (*lattice)->ns==NULL)
        {
          printf(
              "construct_lattice() -- ERROR:  "
              "Attempt to allocate %d struct ns_struct types failed.  "
              "Exiting!\n",
              (*lattice)->NumNodes
              );
          process_exit(1);
        }

#if 0
        // Try to read ns<LX>x<LY>.bmp file.
        sprintf( filename, "./in/ns%dx%d.bmp",
            get_LX(*lattice),
            get_LY(*lattice));
        if( in = fopen( filename, "r+"))
        {
          printf("%s %d >> Reading file \"%s\".\n",__FILE__,__LINE__,filename);
          bmp_read_header( in, &bmih);
          int n;
          for( n=0; n<(*lattice)->NumNodes; n++)
          {
            bmp_read_entry( in, &bmih, &r, &g, &b);

            // Verify grayscale.
            if(    (real)r != (real)g
                || (real)g != (real)b
                || (real)r != (real)b)
            {
              printf(
                  "%s %d >> latman.c: construct_lattice() -- "
                  "n=%d:  [ r g b] = [ %3u %3u %3u]\n",__FILE__,__LINE__,
                  n, (unsigned int)r%256, (unsigned int)g%256, (unsigned int)b%256);
              printf(
                  "%s %d >> latman.c: construct_lattice() -- "
                  "ERROR: File %s needs to be grayscale. "
                  "Exiting!\n",__FILE__,__LINE__, filename);
              process_exit(1);
            }

            // Assign ns value.
            (*lattice)->ns[n].ns = ((real)((unsigned int)r%256))/255.;
#if 0 && VERBOSITY_LEVEL>0
            printf("%s %d >> n=%d, ns=%f\n",
                __FILE__, __LINE__, n, (*lattice)->ns[n].ns);
#endif /* 1 && VERBOSITY_LEVEL>0 */

          } /* for( n=0; n<(*lattice)->NumNodes; n++) */

          fclose( in);

        } /* if( in = fopen( filename, "r+")) */

        else /* !( in = fopen( filename, "r+")) */
        {
          // Can't read ns.bmp file, so use default values.
          printf("%s %d >> WARNING: Can't read \"%s\". "
              "Using default ns values.\n",__FILE__,__LINE__,filename);
        } /* if( in = fopen( filename, "r+")) else */


#else
        // Read solids from .raw file.
        sprintf(filename, "./in/ns%dx%dx%d.raw",
            get_g_LX(*lattice),
            get_g_LY(*lattice),
            get_g_LZ(*lattice) );
        printf("%s %d >> Reading ns values from file \"%s\"\n",
            __FILE__, __LINE__, filename);
        read_ns( *lattice, filename);
#endif
        break;
      }

    case 2:
      {
#if 0
        // Allocate space for ns values.
        (*lattice)->ns =
          (struct ns_struct*)malloc( (*lattice)->NumNodes*sizeof(struct ns_struct));
        if( (*lattice)->ns==NULL)
        {
          printf(
              "construct_lattice() -- ERROR:  "
              "Attempt to allocate %d struct ns_struct types failed.  "
              "Exiting!\n",
              (*lattice)->NumNodes
              );
          process_exit(1);
        }

        // Try to read ns<LX>x<LY>.bmp file.
        sprintf( filename, "./in/ns%dx%d.bmp",
            get_LX(*lattice),
            get_LY(*lattice));
        if( in = fopen( filename, "r+"))
        {
          printf("%s %d >> Reading file \"%s\".\n",__FILE__,__LINE__,filename);
          bmp_read_header( in, &bmih);
          int n;
          for( n=0; n<(*lattice)->NumNodes; n++)
          {
            bmp_read_entry( in, &bmih, &r, &g, &b);

            // Verify grayscale.
            if(    (real)r != (real)g
                || (real)g != (real)b
                || (real)r != (real)b)
            {
              printf(
                  "%s %d >> latman.c: construct_lattice() -- "
                  "n=%d:  [ r g b] = [ %3u %3u %3u]\n",__FILE__,__LINE__,
                  n, (unsigned int)r%256, (unsigned int)g%256, (unsigned int)b%256);
              printf(
                  "%s %d >> latman.c: construct_lattice() -- "
                  "ERROR: File %s needs to be grayscale. "
                  "Exiting!\n",__FILE__,__LINE__, filename);
              process_exit(1);
            }

            if( ((unsigned int)r%256) != 0 && ((unsigned int)r%256) != 255 )
            {
              printf(
                  "%s %d >> latman.c: construct_lattice() -- "
                  "ERROR: File %s needs to be black and white. "
                  "Exiting!\n",__FILE__,__LINE__, filename);
              process_exit(1);
            }

            // Assign ns value.
            if( ((unsigned int)r%256) == 0)
            {
              (*lattice)->ns[n].ns = (*lattice)->param.ns;
            }
            else
            {
              (*lattice)->ns[n].ns = 0.;
            }
#if 0 && VERBOSITY_LEVEL>0
            printf("%s %d >> n=%d, ns=%f\n",
                __FILE__, __LINE__, n, (*lattice)->ns[n].ns);
#endif /* 1 && VERBOSITY_LEVEL>0 */

          } /* for( n=0; n<(*lattice)->NumNodes; n++) */

          fclose( in);

        } /* if( in = fopen( filename, "r+")) */

        else /* !( in = fopen( filename, "r+")) */
        {
          // Can't read ns.bmp file, so use default values.
          printf("%s %d >> WARNING: Can't read \"%s\". "
              "Using default ns values.\n",__FILE__,__LINE__,filename);
        } /* if( in = fopen( filename, "r+")) else */
        break;
#else
        printf("%s %d >> Case ns_flag==2 pending. (Exiting)\n",__FILE__,__LINE__);
        process_exit(1);
#endif
      }

    default:
      {
        printf("%s %d >> construct_lattice() -- Unhandled case: "
            "ns_flag = %d . (Exiting!)\n",
            __FILE__,__LINE__,(*lattice)->param.ns_flag < 0.);
        process_exit(1);
        break;
      }

  } /* switch( (*lattice)->param.ns_flag) */
  (*lattice)->nsterm = (real*)malloc( (*lattice)->NumNodes*Q*sizeof(real));
  if( (*lattice)->nsterm==NULL)
  {
    printf(
        "construct_lattice() -- ERROR:  "
        "Attempt to allocate %d real types for nsterm failed.  "
        "Exiting!\n",
        (*lattice)->NumNodes
        );
    process_exit(1);
  }
#endif /* POROUS_MEDIA */
#endif

  dump_params( *lattice);

} /* void construct_lattice( struct lattice_struct **lattice) */

// void init_problem( struct lattice_struct *lattice)
//##############################################################################
//
// I N I T   P R O B L E M
//
//  - Initialize the problem on the lattice.
//
//    - Set the initial density and velocity.
//
//    - Compute the initial feq.
//
void init_problem( lattice_ptr lattice)
{
  int    a, n, i, j, k, ni, nj, nk;
  real *rho;

#if STORE_UEQ
  real *ueq_x, *ueq_y, *ueq_z;
#endif /* STORE_UEQ */
#if NON_LOCAL_FORCES
  real *force;
#endif /* NON_LOCAL_FORCES */
  int    subs;

#if VERBOSITY_LEVEL > 0
  printf("init_problem() -- Initilizing problem...\n");
#endif /* VERBOSITY_LEVEL > 0 */

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  lattice->time = 0;


  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
    for( n=0; n<get_NumNodes(lattice); n++)
    {
      i = N2I(n,ni,nj,nk);//n%lattice->param.LX;
      j = N2J(n,ni,nj,nk);//n/lattice->param.LX;
      k = N2K(n,ni,nj,nk);

      rho = get_rho_ptr(lattice,subs,n);
      //u_x = get_ux_ptr(lattice,subs,n);
      //u_y = get_uy_ptr(lattice,subs,n);
      //u_z = get_uz_ptr(lattice,subs,n);

#if 0 // ============================ TODO ====================================
#if STORE_UEQ
      ueq_x = &( lattice->ueq[n].u[0]);
      ueq_y = &( lattice->ueq[n].u[1]);
      ueq_z = &( lattice->ueq[n].u[2]);
#endif /* STORE_UEQ */
#endif// ======================================================================


      // Set initial density.
      if( ( 1 || is_not_solid(lattice,n)) )
      {
        switch( lattice->param.initial_condition)
        {
          case 0: // Uniform
            {
              // WHY switch( NUM_FLUID_COMPONENTS)
              // WHY {
              // WHY   case 1:
              // WHY   {
              // WHY     set_rho( lattice, subs, n, get_rho_A(lattice,subs));
              // WHY     break;
              // WHY   }
              // WHY   case 2:
              // WHY   {
              // WHY     set_rho( lattice, subs, n, get_rho_B(lattice,subs));
              // WHY     break;
              // WHY   }
              // WHY }
              set_rho( lattice, subs, n, get_rho_A(lattice,subs));
              break;
            }

#if 0 // ============================ TODO ====================================

          case 1: // Random
            {
              switch( NUM_FLUID_COMPONENTS)
              {
                case 1:
                  {
                    *rho =
                      lattice->param.rho_A[0] + 10.*(real)rand()/(real)RAND_MAX;
                    break;
                  }

                case 2:
                  {
                    if( (real)rand()/(real)RAND_MAX < lattice->param.cut)
                    {
                      *rho = lattice->param.rho_A[subs];
                    }
                    else
                    {
                      *rho = lattice->param.rho_B[subs];
                    }
                    break;
                  }

                default:
                  {
                    printf("%s %d >> ERROR: Unhandled case %d.\n",
                        __FILE__,__LINE__,
                        NUM_FLUID_COMPONENTS);
                    process_exit(1);
                    break;
                  }
              }
              break;
            }
#endif
          case 2: // Cube
            {
              switch( NUM_FLUID_COMPONENTS)
              {
                case 1:
                  {
                    if( i >= lattice->param.x1 && i <= lattice->param.x2
#if PARALLEL
#if NUM_DIMENSIONS == 2
                        && j+ get_proc_id(lattice)*nj >= lattice->param.y1 
                        && j+ get_proc_id(lattice)*nj <= lattice->param.y2
#else
                        && j >= lattice->param.y1 
                        && j <= lattice->param.y2
                        && k+ get_proc_id(lattice)*nk >= lattice->param.z1
                        && k+ get_proc_id(lattice)*nk <= lattice->param.z2
#endif
#else   // !(PARALLEL)
#if NUM_DIMENSIONS == 2
                        && j >= lattice->param.y1 
                        && j <= lattice->param.y2
#else
                        && j >= lattice->param.y1 
                        && j <= lattice->param.y2
                        && k >= lattice->param.z1
                        && k <= lattice->param.z2
#endif
#endif
                       )
                    {
                      *rho = lattice->param.rho_A[0];
                    }
                    else
                    {
                      *rho = lattice->param.rho_B[0];
                    }
                    break;
                  }

                case 2:
                  {
#if INAMURO_SIGMA_COMPONENT
                    if( subs == 0)
                    {
                      *rho = lattice->param.rho_A[0];
                    }
                    else
                    {
                      if( i >= lattice->param.x1 && i <= lattice->param.x2
#if PARALLEL
#if NUM_DIMENSIONS == 2
                        && j+ get_proc_id(lattice)*nj >= lattice->param.y1 
                        && j+ get_proc_id(lattice)*nj <= lattice->param.y2
#else
                        && j >= lattice->param.y1 
                        && j <= lattice->param.y2
                        && k+ get_proc_id(lattice)*nk >= lattice->param.z1
                        && k+ get_proc_id(lattice)*nk <= lattice->param.z2
#endif
#else   // !(PARALLEL)
#if NUM_DIMENSIONS == 2
                        && j >= lattice->param.y1 
                        && j <= lattice->param.y2
#else
                        && j >= lattice->param.y1 
                        && j <= lattice->param.y2
                        && k >= lattice->param.z1
                        && k <= lattice->param.z2
#endif
#endif
                        )
                      {
                        *rho = lattice->param.rho_sigma;
                      }
                      else
                      {
                        *rho = 0.;
                      }
                    }
                    break;

#else
                    if( i >= lattice->param.x1 && i <= lattice->param.x2
#if PARALLEL
#if NUM_DIMENSIONS == 2
                        && j+ get_proc_id(lattice)*nj >= lattice->param.y1 
                        && j+ get_proc_id(lattice)*nj <= lattice->param.y2
#else
                        && j >= lattice->param.y1 
                        && j <= lattice->param.y2
                        && k+ get_proc_id(lattice)*nk >= lattice->param.z1
                        && k+ get_proc_id(lattice)*nk <= lattice->param.z2
#endif
#else   // !(PARALLEL)
#if NUM_DIMENSIONS == 2
                        && j >= lattice->param.y1 
                        && j <= lattice->param.y2
#else
                        && j >= lattice->param.y1 
                        && j <= lattice->param.y2
                        && k >= lattice->param.z1
                        && k <= lattice->param.z2
#endif
#endif
                      )
                    {
                      *rho = lattice->param.rho_A[subs];
                    }
                    else
                    {
                      *rho = lattice->param.rho_B[subs];
                    }
                    break;

#endif
                  }

                default:
                  {
                    printf("%s %d >> ERROR: Unhandled case %d.\n",
                        __FILE__,__LINE__,
                        NUM_FLUID_COMPONENTS);
                    process_exit(1);
                    break;
                  }
              }

              break;
            }
#if 0
          case 3: // Sphere
            {
              switch( NUM_FLUID_COMPONENTS)
              {

                //#if PARALLEL
                //   k+ get_proc_id(lattice)*nk >= lattice->param.z1
                //&& k+ get_proc_id(lattice)*nk <= lattice->param.z2
                ////if 0      k >=k1 && k<=k2
                //#else
                //k >= lattice->param.z1 && k <= lattice->param.z2
                //#endif

                case 1:
                  {
#if PARALLEL
                    if( (i-lattice->param.x0)*(i-lattice->param.x0)
                        + (j-lattice->param.y0)*(j-lattice->param.y0)
                        + (k+get_proc_id(lattice)*nk-lattice->param.z0)
                        *(k+get_proc_id(lattice)*nk-lattice->param.z0)
                        <= lattice->param.r0*lattice->param.r0 )
#else
                      if( (i-lattice->param.x0)*(i-lattice->param.x0)
                          + (j-lattice->param.y0)*(j-lattice->param.y0)
                          + (k-lattice->param.z0)*(k-lattice->param.z0)
                          <= lattice->param.r0*lattice->param.r0 )
#endif
                      {
                        *rho = lattice->param.rho_A[0];
                      }
                      else
                      {
                        *rho = lattice->param.rho_B[0];
                      }
                    break;
                  }
                case 2:
                  {
#if PARALLEL
                    if( (i-lattice->param.x0)*(i-lattice->param.x0)
                        + (j-lattice->param.y0)*(j-lattice->param.y0)
                        + (k+ get_proc_id(lattice)*nk-lattice->param.z0)
                        *(k+ get_proc_id(lattice)*nk-lattice->param.z0)
                        <= lattice->param.r0*lattice->param.r0 )
#else
                      if( (i-lattice->param.x0)*(i-lattice->param.x0)
                          + (j-lattice->param.y0)*(j-lattice->param.y0)
                          + (k-lattice->param.z0)*(k-lattice->param.z0)
                          <= lattice->param.r0*lattice->param.r0 )
#endif
                      {
                        *rho = lattice->param.rho_A[subs];
                      }
                      else
                      {
                        *rho = lattice->param.rho_B[subs];
                      }
                    break;
                  }
                default:
                  {
                    printf(
                        "%s (%d) >> init_problem() -- Unhandled case  "
                        "NUM_FLUID_COMPONENTS = %d.  "
                        "Exiting!\n", __FILE__, __LINE__,
                        NUM_FLUID_COMPONENTS);
                    process_exit(1);
                  }
              }
              break;
            }
          case 4: // Linear Density Gradient
            {
              switch( NUM_FLUID_COMPONENTS)
              {
                case 1:
                  {

                    *rho = lattice->param.rho_in
                      + k*(lattice->param.rho_out-lattice->param.rho_in)/nk;
                    //lattice->param.rho_A[0];
                    break;
                  }

                case 2:
                  {
                    if( i > 20 && i < 32
                        && j > 10 && j < 22
                        && k >  8 && k < 18
                      )
                    {
                      *rho = lattice->param.rho_A[subs];
                    }
                    else
                    {
                      *rho = lattice->param.rho_B[subs];
                    }
                    break;
                  }

                default:
                  {
                    printf("%s %d >> ERROR: Unhandled case %d.\n",
                        __FILE__,__LINE__,
                        NUM_FLUID_COMPONENTS);
                    process_exit(1);
                    break;
                  }
              }

              break;
            }
#endif // =====================================================================

          default:
            {
              printf(
                  "%s (%d) >> init_problem() -- Unhandled case  "
                  "lattice->param.initial_condition = %d.  "
                  "Exiting!\n", __FILE__, __LINE__,
                  lattice->param.initial_condition );
              process_exit(1);
              break;
            }
        } /* switch( lattice->param.initial_condition) */
      } /* if( ( 1 || is_not_solid(lattice,n)) ) */
      else
      {
        //if( is_solid(lattice,n))
        //{
        //  *macro_var_ptr++ = lattice->param.rho_A[subs];
        //  *macro_var_ptr++ = lattice->param.rho_in;
        //}
        //else
        //{
        set_rho(lattice,subs,n,0.);
        //}
      } /* if( ( 1 || is_not_solid(lattice,n)) ) else */

#if 0 // ============================ TODO ====================================

#if SPONGE
      if(// i >= lattice->param.x1 && i <= lattice->param.x2 &&
          // j >= lattice->param.y1 && j <= lattice->param.y2 &&
#if PARALLEL
          k+ get_proc_id(lattice)*nk < lattice->param.z1 || k+ get_proc_id(lattice)*nk > lattice->param.z2)
        //      k >=k1 && k<=k2 )
#else
        k < lattice->param.z1 || k > lattice->param.z2  )
#endif
        {
          *rho = lattice->param.rho_A[subs];
        }
#endif
#endif// ======================================================================

      // Set initial velocity.
#if INITIALIZE_WITH_UX_IN
      set_ux( lattice, subs, n, lattice->param.ux_in);
#else
      set_ux( lattice, subs, n, 0.);
#endif
#if INITIALIZE_WITH_UY_IN
      set_uy( lattice, subs, n, lattice->param.uy_in);
#else
      set_uy( lattice, subs, n, 0.);
#endif
#if INITIALIZE_WITH_UZ_IN
      set_uz( lattice, subs, n, lattice->param.uz_in);
#else
      set_uz( lattice, subs, n, 0.);
#endif

      if( do_diagnostic_init_of_rho(lattice))
      {
        set_rho( lattice, subs, n, get_rho_A(lattice,subs)*(n+1)/10);
      }
      if( do_diagnostic_init_of_u(lattice))
      {
        //set_ux( lattice, subs, n, (real)(n+1)/get_NumNodes(lattice)/100);
        //set_ux( lattice, subs, n, 0.12345);
        set_ux( lattice, subs, n, 0.);

        //set_uy( lattice, subs, n, (real)(n+1)/get_NumNodes(lattice)/100);
        set_uy( lattice, subs, n, 0.00000);

        //set_uz( lattice, subs, n, (real)(n+1)/get_NumNodes(lattice)/100);
        set_uz( lattice, subs, n, 0.);
      }

#if 0 // ============================ TODO ====================================
#if STORE_UEQ
      *ueq_x = 0.;
      *ueq_y = 0.;
      *ueq_z = 0.;
#endif /* STORE_UEQ */
#endif // =====================================================================

    } /* for( n=0; n<lattice->NumNodes; n++) */

  } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

  printf("%s %d >> Before compute_feq...\n",__FILE__,__LINE__);
  for( subs=0; subs<get_NumSubs(lattice); subs++)
  {

    ni = get_ni(lattice);
    nj = get_nj(lattice);
    nk = get_nk(lattice);

    real W0;
    real W1;
    real W2;

    if( get_NumDims(lattice)==2)
    {
      W0 = 4./9.;
      W1 = 1./9.;
      W2 = 1./36.;
    }
    else
    {
      W0 = 1./3.;
      W1 = 1./18.;
      W2 = 1./36.;
    }

    real usq, udotx;

    real** fptr;

    real rho_val;
    real ux_val;
    real uy_val;
    real uz_val;

    fptr = (real**)malloc( sizeof(real*)*get_NumVelDirs(lattice));

    if( !do_diagnostic_init_of_f( lattice))
    {
      n=0;
      for( k=0; k<nk; k++)
      {
        for( j=0; j<nj; j++)
        {
          for( i=0; i<ni; i++)
          {
            fptr[C ] = get_fptr(lattice,subs, i,j,k, 0,0,0, C );
            fptr[E ] = get_fptr(lattice,subs, i,j,k, 0,0,0, E );
            fptr[W ] = get_fptr(lattice,subs, i,j,k, 0,0,0, W );
            fptr[N ] = get_fptr(lattice,subs, i,j,k, 0,0,0, N );
            fptr[S ] = get_fptr(lattice,subs, i,j,k, 0,0,0, S );
            fptr[NE] = get_fptr(lattice,subs, i,j,k, 0,0,0, NE);
            fptr[SW] = get_fptr(lattice,subs, i,j,k, 0,0,0, SW);
            fptr[NW] = get_fptr(lattice,subs, i,j,k, 0,0,0, NW);
            fptr[SE] = get_fptr(lattice,subs, i,j,k, 0,0,0, SE);
            if( get_NumDims(lattice)==3)
            {
              fptr[T ] = get_fptr(lattice,subs, i,j,k, 0,0,0, T );
              fptr[B ] = get_fptr(lattice,subs, i,j,k, 0,0,0, B );
              fptr[TE] = get_fptr(lattice,subs, i,j,k, 0,0,0, TE);
              fptr[BW] = get_fptr(lattice,subs, i,j,k, 0,0,0, BW);
              fptr[TW] = get_fptr(lattice,subs, i,j,k, 0,0,0, TW);
              fptr[BE] = get_fptr(lattice,subs, i,j,k, 0,0,0, BE);
              fptr[TN] = get_fptr(lattice,subs, i,j,k, 0,0,0, TN);
              fptr[BS] = get_fptr(lattice,subs, i,j,k, 0,0,0, BS);
              fptr[TS] = get_fptr(lattice,subs, i,j,k, 0,0,0, TS);
              fptr[BN] = get_fptr(lattice,subs, i,j,k, 0,0,0, BN);
            }

            rho_val = get_rho(lattice,subs,n);
            ux_val  = get_ux (lattice,subs,n);
            uy_val  = get_uy (lattice,subs,n);
            uz_val  = get_uz (lattice,subs,n);

            usq = ux_val*ux_val + uy_val*uy_val + uz_val*uz_val;

#if 1
            *(fptr[0]) = ( /*feq[a]*/
                W0*rho_val*(1. - 1.5*usq)
                ) / get_tau(lattice,subs);

            for( a=1; a<=4; a++)
            {
              udotx = ((real)vx[a]*ux_val
                  +(real)vy[a]*uy_val
                  +(real)vz[a]*uz_val);

              *(fptr[a]) = ( /*feq[a]*/
                  W1*rho_val*(1. + 3.*udotx + 4.5 *udotx*udotx - 1.5*usq)
                  ) / get_tau(lattice,subs);
            }

            for( ; a<=8; a++)
            {
              udotx = ((real)vx[a]*ux_val
                  +(real)vy[a]*uy_val
                  +(real)vz[a]*uz_val);

              *(fptr[a]) = ( /*feq[a]*/
                  W2*rho_val*(1. + 3.*udotx + 4.5 *udotx*udotx - 1.5*usq)
                  ) / get_tau(lattice,subs);
            }

            if( get_NumDims(lattice)==3)
            {
              for( ; a<=10; a++)
              {
                udotx = ((real)vx[a]*ux_val
                    +(real)vy[a]*uy_val
                    +(real)vz[a]*uz_val);

                *(fptr[a]) = ( /*feq[a]*/
                    W1*rho_val*(1. + 3.*udotx + 4.5 *udotx*udotx - 1.5*usq)
                    ) / get_tau(lattice,subs);
              }

              for( ; a<get_NumVelDirs(lattice); a++)
              {
                udotx = ((real)vx[a]*ux_val
                    +(real)vy[a]*uy_val
                    +(real)vz[a]*uz_val);

                *(fptr[a]) = ( /*feq[a]*/
                    W2*rho_val*(1. + 3.*udotx + 4.5 *udotx*udotx - 1.5*usq)
                    ) / get_tau(lattice,subs);
              }
            }
#else
            // Just assign the weighted rho_val for debugging.
            *(fptr[0]) = W0*rho_val;

            for( a=1; a<=4; a++)
            {
              udotx = ((real)vx[a+((a%2)?(1):(-1))]*ux_val
                  +(real)vy[a+((a%2)?(1):(-1))]*uy_val
                  +(real)vz[a+((a%2)?(1):(-1))]*uz_val);

              *(fptr[a]) = W1*rho_val;
            }

            for( ; a<=8; a++)
            {
              udotx = ((real)vx[a+((a%2)?(1):(-1))]*ux_val
                  +(real)vy[a+((a%2)?(1):(-1))]*uy_val
                  +(real)vz[a+((a%2)?(1):(-1))]*uz_val);

              *(fptr[a]) = W2*rho_val;
            }

            if( get_NumDims(lattice)==3)
            {
              for( ; a<=10; a++)
              {
                udotx = ((real)vx[a+((a%2)?(1):(-1))]*ux_val
                    +(real)vy[a+((a%2)?(1):(-1))]*uy_val
                    +(real)vz[a+((a%2)?(1):(-1))]*uz_val);

                *(fptr[a]) = W1*rho_val;
              }

              for( ; a<get_NumVelDirs(lattice); a++)
              {
                udotx = ((real)vx[a+((a%2)?(1):(-1))]*ux_val
                    +(real)vy[a+((a%2)?(1):(-1))]*uy_val
                    +(real)vz[a+((a%2)?(1):(-1))]*uz_val);

                *(fptr[a]) = W2*rho_val;
              }
            }
#endif
            rho_val = 0.;
            ux_val = 0.;
            uy_val = 0.;
            uz_val = 0.;
            for( a=0; a<get_NumVelDirs(lattice); a++)
            {
              rho_val+= (*(fptr[a]));
              ux_val += (*(fptr[a]))*vx[a];
              uy_val += (*(fptr[a]))*vy[a];
              if( get_NumDims(lattice)==3)
              {
                uz_val += (*(fptr[a]))*vz[a];
              }
            }
            ux_val /= rho_val;
            uy_val /= rho_val;
            uz_val /= rho_val;

            n++;
          } // i
        } // j
      } // k
    }
    else
    {
      // Assign canned values for diagnostics purposes. This should be done
      // only with very small domains to help with troubleshooting things like
      // memory management, boundary conditions, ...
      n=0;
      for( k=0; k<nk; k++)
      {
        for( j=0; j<nj; j++)
        {
          for( i=0; i<ni; i++)
          {
            fptr[C ] = get_fptr(lattice,subs, i,j,k, 0,0,0, C );
            fptr[E ] = get_fptr(lattice,subs, i,j,k, 0,0,0, E );
            fptr[W ] = get_fptr(lattice,subs, i,j,k, 0,0,0, W );
            fptr[N ] = get_fptr(lattice,subs, i,j,k, 0,0,0, N );
            fptr[S ] = get_fptr(lattice,subs, i,j,k, 0,0,0, S );
            fptr[NE] = get_fptr(lattice,subs, i,j,k, 0,0,0, NE);
            fptr[SW] = get_fptr(lattice,subs, i,j,k, 0,0,0, SW);
            fptr[NW] = get_fptr(lattice,subs, i,j,k, 0,0,0, NW);
            fptr[SE] = get_fptr(lattice,subs, i,j,k, 0,0,0, SE);
            if( get_NumDims(lattice)==3)
            {
              fptr[T ] = get_fptr(lattice,subs, i,j,k, 0,0,0, T );
              fptr[B ] = get_fptr(lattice,subs, i,j,k, 0,0,0, B );
              fptr[TE] = get_fptr(lattice,subs, i,j,k, 0,0,0, TE);
              fptr[BW] = get_fptr(lattice,subs, i,j,k, 0,0,0, BW);
              fptr[TW] = get_fptr(lattice,subs, i,j,k, 0,0,0, TW);
              fptr[BE] = get_fptr(lattice,subs, i,j,k, 0,0,0, BE);
              fptr[TN] = get_fptr(lattice,subs, i,j,k, 0,0,0, TN);
              fptr[BS] = get_fptr(lattice,subs, i,j,k, 0,0,0, BS);
              fptr[TS] = get_fptr(lattice,subs, i,j,k, 0,0,0, TS);
              fptr[BN] = get_fptr(lattice,subs, i,j,k, 0,0,0, BN);
            }

            *(fptr[C ]) = i/100. + j/10000. + k/1000000. + 0.00000000;
            *(fptr[E ]) = i/100. + j/10000. + k/1000000. + 0.00000001;
            *(fptr[W ]) = i/100. + j/10000. + k/1000000. + 0.00000003;
            *(fptr[N ]) = i/100. + j/10000. + k/1000000. + 0.00000002;
            *(fptr[S ]) = i/100. + j/10000. + k/1000000. + 0.00000004;
            *(fptr[NE]) = i/100. + j/10000. + k/1000000. + 0.00000005;
            *(fptr[SW]) = i/100. + j/10000. + k/1000000. + 0.00000007;
            *(fptr[NW]) = i/100. + j/10000. + k/1000000. + 0.00000006;
            *(fptr[SE]) = i/100. + j/10000. + k/1000000. + 0.00000008;
            if( get_NumDims(lattice)==3)
            {
              // TODO: Assign debug values to 3D directions.
              *(fptr[T ]) = 0.00000000;
              *(fptr[B ]) = 0.00000000;
              *(fptr[TE]) = 0.00000000;
              *(fptr[BW]) = 0.00000000;
              *(fptr[TW]) = 0.00000000;
              *(fptr[BE]) = 0.00000000;
              *(fptr[TN]) = 0.00000000;
              *(fptr[BS]) = 0.00000000;
              *(fptr[TS]) = 0.00000000;
              *(fptr[BN]) = 0.00000000;
            }

#if 0 // TEMP
            if( i==1 && j==1)
            {
              *(fptr[C ]) = 0.99999999;
              *(fptr[E ]) = 1.0;
              *(fptr[W ]) = 2.0;
              *(fptr[N ]) = 3.0;
              *(fptr[S ]) = 4.0;
              *(fptr[NE]) = 5.0;
              *(fptr[SW]) = 6.0;
              *(fptr[NW]) = 7.0;
              *(fptr[SE]) = 8.0;
            }
            else
            {
              *(fptr[C ]) = 0.0;
              *(fptr[E ]) = 0.0;
              *(fptr[W ]) = 0.0;
              *(fptr[N ]) = 0.0;
              *(fptr[S ]) = 0.0;
              *(fptr[NE]) = 0.0;
              *(fptr[SW]) = 0.0;
              *(fptr[NW]) = 0.0;
              *(fptr[SE]) = 0.0;
            }
#endif

            n++;
          } // i
        } // j
      } // k
    }

    free(fptr);

#if 0 // ============================ TODO ====================================
#if NON_LOCAL_FORCES
    for( n=0; n<lattice->NumNodes; n++)
    {
      lattice->force[subs][n].force[0] = 0.;
      lattice->force[subs][n].force[1] = 0.;
      lattice->force[subs][n].force[2] = 0.;
      lattice->force[subs][n].sforce[0] = 0.;
      lattice->force[subs][n].sforce[1] = 0.;
      lattice->force[subs][n].sforce[2] = 0.;
    }
#endif /* NON_LOCAL_FORCES */
#endif// =======================================================================

  } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */
  printf("%s %d >> After compute_feq!\n",__FILE__,__LINE__);

#if VERBOSITY_LEVEL > 0
  printf("init_problem() -- Problem initialized.\n");
#endif /* VERBOSITY_LEVEL > 0 */

} /* void init_problem( lattice_ptr *lattice) */

// void destruct_lattice( struct lattice_struct *lattice)
//##############################################################################
//
// D E S T R U C T   L A T T I C E
//
//  - Destruct lattice.
//
void destruct_lattice( lattice_ptr lattice)
{
#if 0
  int subs;

  assert( lattice!=NULL);

  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
    assert( lattice->pdf[subs]!=NULL);
    free(   lattice->pdf[subs]);

    assert( lattice->macro_vars[subs]!=NULL);
    free(   lattice->macro_vars[subs]);

    assert( lattice->solids[subs]!=NULL);
    free(   lattice->solids[subs]);
#if NON_LOCAL_FORCES
    assert( lattice->force[subs]!=NULL);
    free(   lattice->force[subs]);
#endif /* NON_LOCAL_FORCES */
  }

#if STORE_UEQ
  free(   lattice->ueq);
#endif /* STORE_UEQ */

  process_finalize();

  free(   lattice);

#else
  int subs;
  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
    // Deallocate NumNodes particle distribution functions.
    free( lattice->vars[subs].f_memblock);

    // Deallocate pointers to the individual f arrays.
    free( lattice->vars[subs].f1d);

    // Deallocate NumNodes macroscopic variables.
    free( lattice->vars[subs].macrovars_memblock);

    // Deallocate pointers to the individual macrovar arrays.
    free( lattice->vars[subs].macrovars1d);
  }

  // Deallocate vars struct.
  free( lattice->vars);

  // Deallocate space for solids.
  free( lattice->solids_memblock);

  // Deallocate the lattice structure.
  free( lattice);

  // Free arrays on gpu device

#ifdef __CUDACC__
  cudaFree(f_mem_d);
  cudaFree(mv_mem_d);
  cudaFree(solids_mem_d);
  //cudaFree(is_end_of_frame_mem_d);
#if PAGE_LOCKED
  cudaFreeHost(lattice->process.pos_dir_pdf_to_send);
  cudaFreeHost(lattice->process.pos_dir_pdf_to_recv);
  cudaFreeHost(lattice->process.neg_dir_pdf_to_send);
  cudaFreeHost(lattice->process.neg_dir_pdf_to_recv);

#if BOUNDARY_KERNEL && !(POINTER_MAPPING)
  cudaFree( pos_dir_send_ptr_d);
  cudaFree( pos_dir_recv_ptr_d);
  cudaFree( neg_dir_send_ptr_d);
  cudaFree( neg_dir_recv_ptr_d);
#endif

#else

#if BOUNDARY_KERNEL
  cudaFree( pos_dir_send_ptr_d);
  cudaFree( pos_dir_recv_ptr_d);
  cudaFree( neg_dir_send_ptr_d);
  cudaFree( neg_dir_recv_ptr_d);
#endif

#endif
#endif

  process_finalize();
#endif
} /* void destruct_lattice( lattice_ptr lattice) */

//##############################################################################
void dump_north_pointing_pdfs(
    lattice_ptr lattice,
    const int subs,
    const int z_slice,
    char *comment_str,
    const int which_pdf)
{
#if 0
  int n, p;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  real pdfs[5*get_LX(lattice)*get_LY(lattice)];

  if( domain_is_too_big_to_display( lattice)) { return;}

  for( p=0; p<get_num_procs( lattice); p++)
  {
    process_barrier();
    if( p == get_proc_id(lattice))
    {
      if( z_slice >= 0)
      {
        k = z_slice;

        gather_north_pointing_pdfs( lattice, pdfs, subs, k, which_pdf);

        printf("\n\n// Proc %d, k=%d >> North pointing, \"%s\".",
            get_proc_id( lattice), k, comment_str);
        printf("\n ");
        for( i=0; i<ni; i++)
        {
          printf("+");
          printf("---");
          printf("---");
          printf("---");
          printf("-");
        }
        printf("+");

        for( j=0; j<nj; j++)
        {
          // 0 1 2 3 4
          // O W E N S

          // South
          n = 5*j*ni + 4;
          printf("\n ");
          for( i=0; i<ni; i++)
          {
            printf("|");
            printf("   ");
            printf(" %2.0f", pdfs[n]);
            printf("   ");
            printf(" ");
            n+=5;
          }
          printf("|");

          // West/O/East
          n = 5*j*ni + 0;
          printf("\n ");
          for( i=0; i<ni; i++)
          {
            printf("|");
            printf(" %2.0f", pdfs[n+1]);
            printf(" %2.0f", pdfs[n]);
            printf(" %2.0f", pdfs[n+2]);
            printf(" ");
            n+=5;
          }
          printf("|");

          // North
          n = 5*j*ni + 3;
          printf("\n ");
          for( i=0; i<ni; i++)
          {
            printf("|");
            printf("   ");
            printf(" %2.0f", pdfs[n]);
            printf("   ");
            printf(" ");
            n+=5;
          }
          printf("|");

          printf("\n ");
          for( i=0; i<ni; i++)
          {
            printf("+");
            printf("---");
            printf("---");
            printf("---");
            printf("-");
          }
          printf("+");

        } /* if( j=0; j<nj; j++) */

      } /* if( z_slice >= 0) */

      else
      {
        for( k=0; k<get_LZ(lattice); k++)
        {

          gather_north_pointing_pdfs( lattice, pdfs, subs, k, which_pdf);

          printf("\n\n// Proc %d, k=%d >> North pointing, \"%s\".",
              get_proc_id( lattice), k, comment_str);
          printf("\n ");
          for( i=0; i<ni; i++)
          {
            printf("+");
            printf("---");
            printf("---");
            printf("---");
            printf("-");
          }
          printf("+");

          for( j=0; j<nj; j++)
          {
            // 0 1 2 3 4
            // O W E N S

            // South
            n = 5*j*ni + 4;
            printf("\n ");
            for( i=0; i<ni; i++)
            {
              printf("|");
              printf("   ");
              printf(" %2.0f", pdfs[n]);
              printf("   ");
              printf(" ");
              n+=5;
            }
            printf("|");

            // West/O/East
            n = 5*j*ni + 0;
            printf("\n ");
            for( i=0; i<ni; i++)
            {
              printf("|");
              printf(" %2.0f", pdfs[n+1]);
              printf(" %2.0f", pdfs[n]);
              printf(" %2.0f", pdfs[n+2]);
              printf(" ");
              n+=5;
            }
            printf("|");

            // North
            n = 5*j*ni + 3;
            printf("\n ");
            for( i=0; i<ni; i++)
            {
              printf("|");
              printf("   ");
              printf(" %2.0f", pdfs[n]);
              printf("   ");
              printf(" ");
              n+=5;
            }
            printf("|");

            printf("\n ");
            for( i=0; i<ni; i++)
            {
              printf("+");
              printf("---");
              printf("---");
              printf("---");
              printf("-");
            }
            printf("+");

          } /* if( j=0; j<nj; j++) */

        }

      } /* if( z_slice >= 0) else */

    }

  } /* for( p=0; p<get_num_procs( lattice); p++) */
  process_barrier();

#else
  // TODO
#endif
} /* void dump_north_pointing_pdfs( lattice_ptr lattice, const int subs, ... */

//##############################################################################
void dump_south_pointing_pdfs(
    lattice_ptr lattice,
    const int subs,
    const int z_slice,
    char *comment_str,
    const int which_pdf)
{
#if 0
  int n, p;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  real pdfs[5*get_LX(lattice)*get_LY(lattice)];

  if( domain_is_too_big_to_display( lattice)) { return;}

  for( p=0; p<get_num_procs( lattice); p++)
  {
    process_barrier();
    if( p == get_proc_id(lattice))
    {
      if( z_slice >= 0)
      {
        k = z_slice;

        gather_south_pointing_pdfs( lattice, pdfs, subs, k, which_pdf);

        printf("\n\n// Proc %d, k=%d >> South pointing, \"%s\".",
            get_proc_id( lattice), k, comment_str);
        printf("\n ");
        for( i=0; i<ni; i++)
        {
          printf("+");
          printf("---");
          printf("---");
          printf("---");
          printf("-");
        }
        printf("+");

        for( j=0; j<nj; j++)
        {
          // 0 1 2 3 4
          // O W E N S

          // South
          n = 5*j*ni + 4;
          printf("\n ");
          for( i=0; i<ni; i++)
          {
            printf("|");
            printf("   ");
            printf(" %2.0f", pdfs[n]);
            printf("   ");
            printf(" ");
            n+=5;
          }
          printf("|");

          // West/O/East
          n = 5*j*ni + 0;
          printf("\n ");
          for( i=0; i<ni; i++)
          {
            printf("|");
            printf(" %2.0f", pdfs[n+1]);
            printf(" %2.0f", pdfs[n]);
            printf(" %2.0f", pdfs[n+2]);
            printf(" ");
            n+=5;
          }
          printf("|");

          // North
          n = 5*j*ni + 3;
          printf("\n ");
          for( i=0; i<ni; i++)
          {
            printf("|");
            printf("   ");
            printf(" %2.0f", pdfs[n]);
            printf("   ");
            printf(" ");
            n+=5;
          }
          printf("|");

          printf("\n ");
          for( i=0; i<ni; i++)
          {
            printf("+");
            printf("---");
            printf("---");
            printf("---");
            printf("-");
          }
          printf("+");

        } /* if( j=0; j<nj; j++) */

      } /* if( z_slice >= 0) */

      else
      {
        for( k=0; k<get_LZ(lattice); k++)
        {

          gather_south_pointing_pdfs( lattice, pdfs, subs, k, which_pdf);

          printf("\n\n// Proc %d, k=%d >> South pointing, \"%s\".",
              get_proc_id( lattice), k, comment_str);
          printf("\n ");
          for( i=0; i<ni; i++)
          {
            printf("+");
            printf("---");
            printf("---");
            printf("---");
            printf("-");
          }
          printf("+");

          for( j=0; j<nj; j++)
          {
            // 0 1 2 3 4
            // O W E N S

            // South
            n = 5*j*ni + 4;
            printf("\n ");
            for( i=0; i<ni; i++)
            {
              printf("|");
              printf("   ");
              printf(" %2.0f", pdfs[n]);
              printf("   ");
              printf(" ");
              n+=5;
            }
            printf("|");

            // West/O/East
            n = 5*j*ni + 0;
            printf("\n ");
            for( i=0; i<ni; i++)
            {
              printf("|");
              printf(" %2.0f", pdfs[n+1]);
              printf(" %2.0f", pdfs[n]);
              printf(" %2.0f", pdfs[n+2]);
              printf(" ");
              n+=5;
            }
            printf("|");

            // North
            n = 5*j*ni + 3;
            printf("\n ");
            for( i=0; i<ni; i++)
            {
              printf("|");
              printf("   ");
              printf(" %2.0f", pdfs[n]);
              printf("   ");
              printf(" ");
              n+=5;
            }
            printf("|");

            printf("\n ");
            for( i=0; i<ni; i++)
            {
              printf("+");
              printf("---");
              printf("---");
              printf("---");
              printf("-");
            }
            printf("+");

          } /* if( j=0; j<nj; j++) */

        }

      } /* if( z_slice >= 0) else */

    }

  } /* for( p=0; p<get_num_procs( lattice); p++) */
  process_barrier();

#else
  // TODO
#endif
} /* void dump_south_pointing_pdfs( lattice_ptr lattice, const int subs, ... */

int domain_is_too_big_to_display( lattice_ptr lattice)
{
  if(   get_LX(lattice) > 12
      || get_LY(lattice) > 12
      || get_LZ(lattice) > 12 ) { return 1;} else { return 0;}
}

int domain_is_not_too_big_to_display( lattice_ptr lattice)
{
  return !domain_is_too_big_to_display(lattice);
}

void display_warning_about_contrived_data( lattice_ptr lattice)
{
  printf("\n");
  printf("\n");
  printf("DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG");
  printf("DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG");
  printf("\n");
  printf("\n");
  printf(" W A R N I N G:");
  printf("\n");
  printf("\n");
  printf("      Contrived data! This run good for DEBUG purposes only.");
  printf("\n");
  printf("\n");
  printf("DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG");
  printf("DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG");
  printf("\n");
  printf("\n");
}

// int get_sizeof_lattice_structure( lattice_ptr lattice)
//##############################################################################
//
// G E T _ S I Z E O F _ L A T T I C E _ S T R U C T U R E
//
//  - Return size of struct lattice_struct in bytes.
//
int get_sizeof_lattice_structure( lattice_ptr lattice)
{
#if 0
  return sizeof( struct lattice_struct);

#else
  // TODO
  return -1;
#endif
} /* int get_sizeof_lattice_structure( lattice_ptr lattice) */

// int get_sizeof_lattice( lattice_ptr lattice)
//##############################################################################
//
// G E T _ S I Z E O F _ L A T T I C E
//
//  - Return size of lattice in bytes.
//
int get_sizeof_lattice( lattice_ptr lattice)
{
#if 0
  return
    sizeof(int) /* NumNodes*/
    + sizeof(int) /* NumTimeSteps */
    + sizeof(int) /* time */
    + sizeof(int)*NUM_FLUID_COMPONENTS /* periodic_x */
    + sizeof(int)*NUM_FLUID_COMPONENTS /* periodic_y */

#if INAMURO_SIGMA_COMPONENT
    + sizeof(int) /* SizeBTC */
    + sizeof(int) /* FlowDir */
#endif

    + sizeof(struct param_struct)

    + lattice->NumNodes
    * (
        NUM_FLUID_COMPONENTS*sizeof(struct pdf_struct)
        + NUM_FLUID_COMPONENTS*sizeof(struct macro_vars_struct)
        + NUM_FLUID_COMPONENTS*sizeof(struct solids_struct)
#if NON_LOCAL_FORCES
        + NUM_FLUID_COMPONENTS*sizeof(struct force_struct)
#endif /* NON_LOCAL_FORCES */
#if STORE_UEQ
        + sizeof(struct ueq_struct)
#endif /* STORE_UEQ */
#if POROUS_MEDIA
        + sizeof(struct ns_struct)
        + Q*sizeof(real) /* nsterm */
#endif /* POROUS_MEDIA */
        + sizeof(struct process_struct)
      );

#else
  // TODO
  return -1;
#endif
} /* int get_sizeof_lattice( lattice_ptr lattice) */

