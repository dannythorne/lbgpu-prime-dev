//##############################################################################
//
// process.c
//
//##############################################################################
void process_init( lattice_ptr lattice, int argc, char **argv)
{
#if PARALLEL
  // Routine inititialization calls.
  MPI_Init( &argc, &argv);
  MPI_Comm_size( MPI_COMM_WORLD, &(lattice->process.num_procs));
  MPI_Comm_rank( MPI_COMM_WORLD, &(lattice->process.id));
#else
  lattice->process.id = 0;
  lattice->process.num_procs = 1;
#endif
#if VERBOSITY_LEVEL > 0
  // Say hi.
  printf("Hello >>  ProcID = %d, NumProcs = %d.\n",
      get_proc_id( lattice),
      get_num_procs( lattice) );
#endif
} /* void process_init( lattice_ptr lattice, int argc, char **argv) */

//##############################################################################
void process_compute_local_params( lattice_ptr lattice)
{
#if PARALLEL
  int NumLayersOnRoot;
  int NumLayersPerProc;

  // Save a copy of global dimensions.
  set_g_LX( lattice, get_LX( lattice));
  set_g_LY( lattice, get_LY( lattice));
  set_g_LZ( lattice, get_LZ( lattice));
  set_g_SX( lattice, 0);
  set_g_SY( lattice, 0);
  set_g_SZ( lattice, 0);
  set_g_EX( lattice, get_LX( lattice) - 1);
  set_g_EY( lattice, get_LY( lattice) - 1);
  set_g_EZ( lattice, get_LZ( lattice) - 1);
  set_g_NumNodes( lattice, get_NumNodes( lattice));

  // Adjust local z-dimension according to local subdomain.
  // Currently only supports partitioning in y-direction
  // for 2D and the z-direction for 3D.  Note however that
  // this setup may be optimal for the GPU - memory locations
  // are adjacent for the x-direction, leading to minimal
  // number of memcpys and in the case of calculation, memory
  // coalescense.  The minimum number of communications is only
  // 1/3 less than that used here, and the cost of inefficient
  // memcpys/non-coalescence is presumably higher.

#if NUM_DIMENSIONS == 2
  NumLayersOnRoot = get_g_LY( lattice) % get_num_procs( lattice);
  if( NumLayersOnRoot != 0)
  {
    NumLayersPerProc =
      ( get_g_LY( lattice) - NumLayersOnRoot) / (get_num_procs(lattice)-1);
  }
  else
  {
    NumLayersPerProc =
      ( get_g_LY( lattice) - NumLayersOnRoot) / (get_num_procs(lattice));
    NumLayersOnRoot = NumLayersPerProc;
  }

  if( is_on_root_proc( lattice))
  {
    // Assign the left-over (modulus) layers.
    set_LY( lattice, NumLayersOnRoot);
    set_g_SY( lattice, 0);
    set_g_EY( lattice, 0 + NumLayersOnRoot - 1);
    set_g_StartNode( lattice, 0);
  }
  else
  {
    set_LY( lattice, NumLayersPerProc);
    set_g_SY( lattice,
	NumLayersOnRoot
	+ NumLayersPerProc*(get_proc_id(lattice)-1) );
    set_g_EY( lattice,
	NumLayersOnRoot
	+ NumLayersPerProc*(get_proc_id(lattice)-1)
	+ NumLayersPerProc
	- 1);

    set_g_StartNode( lattice,
	get_g_SY( lattice)*( get_LX(lattice)) );
  }

  int pdf_buffer_size = 3 * get_LX( lattice);
  int rho_buffer_size = 2 * get_LX( lattice);
  int sol_buffer_size = 1 * get_LX( lattice); 

#else
  NumLayersOnRoot = get_g_LZ( lattice) % get_num_procs( lattice);
  if( NumLayersOnRoot != 0)
  {
    NumLayersPerProc =
      ( get_g_LZ( lattice) - NumLayersOnRoot) / (get_num_procs(lattice)-1);
  }
  else
  {
    NumLayersPerProc =
      ( get_g_LZ( lattice) - NumLayersOnRoot) / (get_num_procs(lattice));
    NumLayersOnRoot = NumLayersPerProc;
  }

  if( is_on_root_proc( lattice))
  {
    // Assign the left-over (modulus) layers.
    set_LZ( lattice, NumLayersOnRoot);
    set_g_SZ( lattice, 0);
    set_g_EZ( lattice, 0 + NumLayersOnRoot - 1);
    set_g_StartNode( lattice, 0);
  }
  else
  {
    set_LZ( lattice, NumLayersPerProc);
    set_g_SZ( lattice,
	NumLayersOnRoot
	+ NumLayersPerProc*(get_proc_id(lattice)-1) );
    set_g_EZ( lattice,
	NumLayersOnRoot
	+ NumLayersPerProc*(get_proc_id(lattice)-1)
	+ NumLayersPerProc
	- 1);

    set_g_StartNode( lattice,
	get_g_SZ( lattice)*( get_LX(lattice)*get_LY(lattice)) );
  }

  int pdf_buffer_size = 5 * get_LX( lattice) * get_LY( lattice);
  int rho_buffer_size = 2 * get_LX( lattice) * get_LY( lattice);
  int sol_buffer_size = 1 * get_LX( lattice) * get_LY( lattice); 

#endif


  set_NumNodes( lattice);



  lattice->process.pos_dir_pdf_to_send =
    (real*)malloc( pdf_buffer_size*sizeof(real));
  lattice->process.pos_dir_pdf_to_recv =
    (real*)malloc( pdf_buffer_size*sizeof(real));
  lattice->process.neg_dir_pdf_to_send =
    (real*)malloc( pdf_buffer_size*sizeof(real));
  lattice->process.neg_dir_pdf_to_recv =
    (real*)malloc( pdf_buffer_size*sizeof(real));

  lattice->process.pos_dir_rho_to_send =
    (real*)malloc( rho_buffer_size*sizeof(real));
  lattice->process.pos_dir_rho_to_recv =
    (real*)malloc( rho_buffer_size*sizeof(real));
  lattice->process.neg_dir_rho_to_send =
    (real*)malloc( rho_buffer_size*sizeof(real));
  lattice->process.neg_dir_rho_to_recv =
    (real*)malloc( rho_buffer_size*sizeof(real));

  lattice->process.pos_dir_solid_to_send =
    (int*)malloc( sol_buffer_size*sizeof(int));
  lattice->process.pos_dir_solid_to_recv =
    (int*)malloc( sol_buffer_size*sizeof(int));
  lattice->process.neg_dir_solid_to_send =
    (int*)malloc( sol_buffer_size*sizeof(int));
  lattice->process.neg_dir_solid_to_recv =
    (int*)malloc( sol_buffer_size*sizeof(int));

#endif

#if VERBOSITY_LEVEL > 0
#if PARALLEL
  printf(
      "Proc %04d"
      ", g_SX = %d"
      ", g_EX = %d"
      ", g_LX = %d"
      ", g_SY = %d"
      ", g_EY = %d"
      ", g_LY = %d"
      ", g_SZ = %d"
      ", g_EZ = %d"
      ", g_LZ = %d"
      ", g_StartNode = %d"
      ", g_NumNodes = %d"
      ".\n"
      ,get_proc_id( lattice)
      ,get_g_SX( lattice)
      ,get_g_EX( lattice)
      ,get_g_LX( lattice)
      ,get_g_SY( lattice)
      ,get_g_EY( lattice)
      ,get_g_LY( lattice)
      ,get_g_SZ( lattice)
      ,get_g_EZ( lattice)
      ,get_g_LZ( lattice)
      ,get_g_StartNode( lattice)
      ,get_g_NumNodes( lattice)
      );
#endif
#endif

} /* void process_init( process_ptr *process, int argc, char **argv) */

//##############################################################################
void process_send_recv_begin( lattice_ptr lattice, const int subs)
{
#if PARALLEL
  int n;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  int mpierr;

  // A C C U M U L A T E   P D F S   T O   S E N D
  //#########################################################################
  n = 0;
  k = get_LZ(lattice)-1;
  for( j=0; j<nj; j++)
  {
    for( i=0; i<ni; i++)
    {
#if SAVE_MEMO
#else
      lattice->process.pos_dir_pdf_to_send[n] =
	lattice->pdf[subs][ IJK2N( i , j , k, ni, nj)].f[ T];
      lattice->process.neg_dir_pdf_to_send[n] =
	lattice->pdf[subs][ IJK2N( i , j , 0, ni, nj)].f[ B];
      n++;
      lattice->process.pos_dir_pdf_to_send[n] =
	lattice->pdf[subs][ IJK2N( i , j , k, ni, nj)].f[TW];
      lattice->process.neg_dir_pdf_to_send[n] =
	lattice->pdf[subs][ IJK2N( i , j , 0, ni, nj)].f[BW];
      n++;
      lattice->process.pos_dir_pdf_to_send[n] =
	lattice->pdf[subs][ IJK2N( i , j , k, ni, nj)].f[TE];
      lattice->process.neg_dir_pdf_to_send[n] =
	lattice->pdf[subs][ IJK2N( i , j , 0, ni, nj)].f[BE];
      n++;
      lattice->process.pos_dir_pdf_to_send[n] =
	lattice->pdf[subs][ IJK2N( i , j , k, ni, nj)].f[TN];
      lattice->process.neg_dir_pdf_to_send[n] =
	lattice->pdf[subs][ IJK2N( i , j , 0, ni, nj)].f[BN];
      n++;
      lattice->process.pos_dir_pdf_to_send[n] =
	lattice->pdf[subs][ IJK2N( i , j , k, ni, nj)].f[TS];
      lattice->process.neg_dir_pdf_to_send[n] =
	lattice->pdf[subs][ IJK2N( i , j , 0, ni, nj)].f[BS];
      n++;
#endif
    } /* if( i=0; i<ni; i++) */
  } /* if( j=0; j<nj; j++) */


#if 0
  // Contrived debug data...
  display_warning_about_contrived_data( lattice);
  n = 0;
  k = 1;
  for( j=0; j<nj; j++)
  {
    for( i=0; i<ni; i++)
    {
#if 0
      lattice->process.pos_dir_pdf_to_send[n] = 1;
      lattice->process.neg_dir_pdf_to_send[n] = 1;
      n++;
      lattice->process.pos_dir_pdf_to_send[n] = 2;
      lattice->process.neg_dir_pdf_to_send[n] = 2;
      n++;
      lattice->process.pos_dir_pdf_to_send[n] = 3;
      lattice->process.neg_dir_pdf_to_send[n] = 3;
      n++;
      lattice->process.pos_dir_pdf_to_send[n] = 4;
      lattice->process.neg_dir_pdf_to_send[n] = 4;
      n++;
      lattice->process.pos_dir_pdf_to_send[n] = 5;
      lattice->process.neg_dir_pdf_to_send[n] = 5;
      n++;
#endif
#if 1
      lattice->process.pos_dir_pdf_to_send[n] = k;
      lattice->process.neg_dir_pdf_to_send[n] = k;
      n++;
      lattice->process.pos_dir_pdf_to_send[n] = k;
      lattice->process.neg_dir_pdf_to_send[n] = k;
      n++;
      lattice->process.pos_dir_pdf_to_send[n] = k;
      lattice->process.neg_dir_pdf_to_send[n] = k;
      n++;
      lattice->process.pos_dir_pdf_to_send[n] = k;
      lattice->process.neg_dir_pdf_to_send[n] = k;
      n++;
      lattice->process.pos_dir_pdf_to_send[n] = k;
      lattice->process.neg_dir_pdf_to_send[n] = k;
      n++;

      k++;
#endif
    } /* if( i=0; i<ni; i++) */
  } /* if( j=0; j<nj; j++) */
#endif

  //process_dump_pdfs_to_send( lattice, "Before Send/Recv");

  // S E N D   I N   P O S I T I V E   D I R E C T I O N
  //#########################################################################
#if VERBOSITY_LEVEL > 1
  printf( "%s %d %04d >> "
      "MPI_Isend( %04d)"
      "\n",
      __FILE__,__LINE__,get_proc_id(lattice),
      (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
#endif
  mpierr =
    MPI_Isend(
	/*void *buf*/          lattice->process.pos_dir_pdf_to_send,
	/*int count*/          pdf_buffer_size,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*int dest*/         ( get_proc_id(lattice)
	  + get_num_procs(lattice)+1)
	% get_num_procs(lattice),
	/*int tag*/            0,
	/*MPI_Comm comm*/      MPI_COMM_WORLD,
	/*MPI_Request *req*/   &(lattice->process.send_req_0)
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Isend( %04d)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice),
	mpierr,
	(get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
    process_finalize();
    exit(1);
  }
  // R E C V   F R O M   N E G A T I V E   D I R E C T I O N
  //#########################################################################
#if VERBOSITY_LEVEL > 1
  printf( "%s %d %04d >> "
      "MPI_Irecv( %04d)"
      "\n",
      __FILE__,__LINE__,get_proc_id(lattice),
      (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
#endif
  mpierr =
    MPI_Irecv(
	/*void *buf*/          lattice->process.pos_dir_pdf_to_recv,
	/*int count*/          pdf_buffer_size,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*int src*/          ( get_proc_id(lattice)
	  + get_num_procs(lattice)-1)
	% get_num_procs(lattice),
	/*int tag*/            0,
	/*MPI_Comm comm*/      MPI_COMM_WORLD,
	/*MPI_Request *req*/   &(lattice->process.recv_req_0)
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Irecv( %04d)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice),
	mpierr,
	(get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
    process_finalize();
    exit(1);
  }
  // S E N D   I N   N E G A T I V E   D I R E C T I O N
  //#########################################################################
#if VERBOSITY_LEVEL > 1
  printf( "%s %d %04d >> "
      "MPI_Isend( %04d)"
      "\n",
      __FILE__,__LINE__,get_proc_id(lattice),
      (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
#endif
  mpierr =
    MPI_Isend(
	/*void *buf*/          lattice->process.neg_dir_pdf_to_send,
	/*int count*/          pdf_buffer_size,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*int dest*/         ( get_proc_id(lattice)
	  + get_num_procs(lattice)-1)
	% get_num_procs(lattice),
	/*int tag*/            1,
	/*MPI_Comm comm*/      MPI_COMM_WORLD,
	/*MPI_Request *req*/   &(lattice->process.send_req_1)
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Isend( %04d)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice),
	mpierr,
	(get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
    process_finalize();
    exit(1);
  }
  // R E C V   F R O M   P O S I T I V E   D I R E C T I O N
  //#########################################################################
#if VERBOSITY_LEVEL > 1
  printf( "%s %d %04d >> "
      "MPI_Irecv( %04d)"
      "\n",
      __FILE__,__LINE__,get_proc_id(lattice),
      (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
#endif
  mpierr =
    MPI_Irecv(
	/*void *buf*/          lattice->process.neg_dir_pdf_to_recv,
	/*int count*/          pdf_buffer_size,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*int src*/          ( get_proc_id(lattice)
	  + get_num_procs(lattice)+1)
	% get_num_procs(lattice),
	/*int tag*/            1,
	/*MPI_Comm comm*/      MPI_COMM_WORLD,
	/*MPI_Request *req*/   &(lattice->process.recv_req_1)
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Irecv( %04d)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice),
	mpierr,
	(get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
    process_finalize();
    exit(1);
  }
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




#endif
} /* void process_send_recv_begin( lattice_ptr lattice, const int subs) */

//##############################################################################
void process_send_recv_end( lattice_ptr lattice, const int subs)
{
#if PARALLEL
  int n;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  int ip, in;
  int jp, jn;
  int mpierr;
  mpierr = MPI_Wait(
      /* MPI_Request *req */&(lattice->process.send_req_0),
      /* MPI_Status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
      /* MPI_Request *req */&(lattice->process.recv_req_0),
      /* MPI_Status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
      /* MPI_Request *req */&(lattice->process.send_req_1),
      /* MPI_Status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
      /* MPI_Request *req */&(lattice->process.recv_req_1),
      /* MPI_Status *stat */&(lattice->process.mpi_status));

  //process_dump_pdfs_to_recv( lattice, "After Send/Recv, Before Stream");
  //dump_north_pointing_pdfs( lattice, subs, -1,
  //    "After Send/Recv, Before Stream", 2);
  //dump_south_pointing_pdfs( lattice, subs, -1,
  //    "After Send/Recv, Before Stream", 2);

  // S T R E A M   I N   T H E   B O U N D A R I E S
  //###########################################################################
  n = 0;
  k = get_LZ(lattice)-1;
  for( j=0; j<nj; j++)
  {
    jp = ( j<nj-1)?( j+1):( 0   );
    jn = ( j>0   )?( j-1):( nj-1);

    for( i=0; i<ni; i++)
    {
      ip = ( i<ni-1)?( i+1):( 0   );
      in = ( i>0   )?( i-1):( ni-1);

      lattice->pdf[subs][ IJK2N( i , j , 0, ni, nj)].ftemp[ T] =
	lattice->process.pos_dir_pdf_to_recv[n];
      lattice->pdf[subs][ IJK2N( i , j , k, ni, nj)].ftemp[ B] =
	lattice->process.neg_dir_pdf_to_recv[n];
      n++;
      lattice->pdf[subs][ IJK2N( in, j , 0, ni, nj)].ftemp[TW] =
	lattice->process.pos_dir_pdf_to_recv[n];
      lattice->pdf[subs][ IJK2N( in, j , k, ni, nj)].ftemp[BW] =
	lattice->process.neg_dir_pdf_to_recv[n];
      n++;
      lattice->pdf[subs][ IJK2N( ip, j , 0, ni, nj)].ftemp[TE] =
	lattice->process.pos_dir_pdf_to_recv[n];
      lattice->pdf[subs][ IJK2N( ip, j , k, ni, nj)].ftemp[BE] =
	lattice->process.neg_dir_pdf_to_recv[n];
      n++;
      lattice->pdf[subs][ IJK2N( i , jp, 0, ni, nj)].ftemp[TN] =
	lattice->process.pos_dir_pdf_to_recv[n];
      lattice->pdf[subs][ IJK2N( i , jp, k, ni, nj)].ftemp[BN] =
	lattice->process.neg_dir_pdf_to_recv[n];
      n++;
      lattice->pdf[subs][ IJK2N( i , jn, 0, ni, nj)].ftemp[TS] =
	lattice->process.pos_dir_pdf_to_recv[n];
      lattice->pdf[subs][ IJK2N( i , jn, k, ni, nj)].ftemp[BS] =
	lattice->process.neg_dir_pdf_to_recv[n];
      n++;

    } /* if( i=0; i<ni; i++) */
  } /* if( j=0; j<nj; j++) */

#if 0
  // Copy back to check with a call to process_dump_pdfs...
#if 0
  n = 0;
  k = get_LZ(lattice)-1;
  for( j=0; j<nj; j++)
  {
    for( i=0; i<ni; i++)
    {
      lattice->process.pos_dir_pdf_to_recv[n] =
	lattice->pdf[subs][ IJK2N( i , j , k, ni, nj)].ftemp[ T];
      lattice->process.neg_dir_pdf_to_recv[n] =
	lattice->pdf[subs][ IJK2N( i , j , 0, ni, nj)].ftemp[ B];
      n++;
      lattice->process.pos_dir_pdf_to_recv[n] =
	lattice->pdf[subs][ IJK2N( i , j , k, ni, nj)].ftemp[TW];
      lattice->process.neg_dir_pdf_to_recv[n] =
	lattice->pdf[subs][ IJK2N( i , j , 0, ni, nj)].ftemp[BW];
      n++;
      lattice->process.pos_dir_pdf_to_recv[n] =
	lattice->pdf[subs][ IJK2N( i , j , k, ni, nj)].ftemp[TE];
      lattice->process.neg_dir_pdf_to_recv[n] =
	lattice->pdf[subs][ IJK2N( i , j , 0, ni, nj)].ftemp[BE];
      n++;
      lattice->process.pos_dir_pdf_to_recv[n] =
	lattice->pdf[subs][ IJK2N( i , j , k, ni, nj)].ftemp[TN];
      lattice->process.neg_dir_pdf_to_recv[n] =
	lattice->pdf[subs][ IJK2N( i , j , 0, ni, nj)].ftemp[BN];
      n++;
      lattice->process.pos_dir_pdf_to_recv[n] =
	lattice->pdf[subs][ IJK2N( i , j , k, ni, nj)].ftemp[TS];
      lattice->process.neg_dir_pdf_to_recv[n] =
	lattice->pdf[subs][ IJK2N( i , j , 0, ni, nj)].ftemp[BS];
      n++;

    } /* if( i=0; i<ni; i++) */
  } /* if( j=0; j<nj; j++) */
#else
  gather_north_pointing_pdfs( lattice,
      lattice->process.pos_dir_pdf_to_recv,
      subs, get_LZ(lattice)-1, 2);
  gather_south_pointing_pdfs( lattice,
      lattice->process.neg_dir_pdf_to_recv,
      subs, 0, 2);
#endif

  //process_dump_pdfs_to_recv( lattice, "After Stream");

#endif

#endif
} /* void process_send_recv_end( lattice_ptr lattice, const int subs) */



void rho_send_recv_begin( lattice_ptr lattice, const int subs)
{
#if PARALLEL
  int n, a;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  int mpierr;

  //rev Huang
  //To send "rho" and "solid" to neighbour blocks
  n =0;
  k = get_LZ(lattice)-1;
  // for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  for( a=0; a<NUM_FLUID_COMPONENTS; a++)
    //changed by Shadab on Danny's suggestion!!
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
#if SAVE_MEMO

#else

	lattice->process.pos_dir_rho_to_send[n] =
	  lattice->macro_vars[a][ IJK2N( i , j , k, ni, nj)].rho;
	lattice->process.neg_dir_rho_to_send[n] =
	  lattice->macro_vars[a][ IJK2N( i , j , 0, ni, nj)].rho;
	n++;
#endif
      } /* if( i=0; i<ni; i++) */
    } /* if( j=0; j<nj; j++) */
  } //for  subs=0; subs<NUM_FLUID_COMPONENTS; subs++

  // S E N D   I N   P O S I T I V E   D I R E C T I O N
  //#########################################################################
#if VERBOSITY_LEVEL > 1
  printf( "%s %d %04d >> "
      "MPI_Isend( %04d)"
      "\n",
      __FILE__,__LINE__,get_proc_id(lattice),
      (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
#endif
  mpierr =
    MPI_Isend(
	/*void *buf*/          lattice->process.pos_dir_rho_to_send,
	/*int count*/          rho_buffer_size,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*int dest*/         ( get_proc_id(lattice)
	  + get_num_procs(lattice)+1)
	% get_num_procs(lattice),
	/*int tag*/            2,
	/*MPI_Comm comm*/      MPI_COMM_WORLD,
	/*MPI_Request *req*/   &(lattice->process.send_req_2)
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Isend( %04d)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice),
	mpierr,
	(get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
    process_finalize();
    process_exit(1);
  }
  // R E C V   F R O M   N E G A T I V E   D I R E C T I O N
  //#########################################################################
#if VERBOSITY_LEVEL > 1
  printf( "%s %d %04d >> "
      "MPI_Irecv( %04d)"
      "\n",
      __FILE__,__LINE__,get_proc_id(lattice),
      (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
#endif
  mpierr =
    MPI_Irecv(
	/*void *buf*/          lattice->process.pos_dir_rho_to_recv,
	/*int count*/          rho_buffer_size,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*int src*/          ( get_proc_id(lattice)
	  + get_num_procs(lattice)-1)
	% get_num_procs(lattice),
	/*int tag*/            2,
	/*MPI_Comm comm*/      MPI_COMM_WORLD,
	/*MPI_Request *req*/   &(lattice->process.recv_req_2)
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Irecv( %04d)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice),
	mpierr,
	(get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
    process_finalize();
    process_exit(1);
  }
  // S E N D   I N   N E G A T I V E   D I R E C T I O N
  //#########################################################################
#if VERBOSITY_LEVEL > 1
  printf( "%s %d %04d >> "
      "MPI_Isend( %04d)"
      "\n",
      __FILE__,__LINE__,get_proc_id(lattice),
      (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
#endif
  mpierr =
    MPI_Isend(
	/*void *buf*/          lattice->process.neg_dir_rho_to_send,
	/*int count*/          rho_buffer_size,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*int dest*/         ( get_proc_id(lattice)
	  + get_num_procs(lattice)-1)
	% get_num_procs(lattice),
	/*int tag*/            3,
	/*MPI_Comm comm*/      MPI_COMM_WORLD,
	/*MPI_Request *req*/   &(lattice->process.send_req_3)
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Isend( %04d)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice),
	mpierr,
	(get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
    process_finalize();
    process_exit(1);
  }
  // R E C V   F R O M   P O S I T I V E   D I R E C T I O N
  //#########################################################################
#if VERBOSITY_LEVEL > 1
  printf( "%s %d %04d >> "
      "MPI_Irecv( %04d)"
      "\n",
      __FILE__,__LINE__,get_proc_id(lattice),
      (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
#endif
  mpierr =
    MPI_Irecv(
	/*void *buf*/          lattice->process.neg_dir_rho_to_recv,
	/*int count*/          rho_buffer_size,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*int src*/          ( get_proc_id(lattice)
	  + get_num_procs(lattice)+1)
	% get_num_procs(lattice),
	/*int tag*/            3,
	/*MPI_Comm comm*/      MPI_COMM_WORLD,
	/*MPI_Request *req*/   &(lattice->process.recv_req_3)
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Irecv( %04d)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice),
	mpierr,
	(get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
    process_finalize();
    process_exit(1);
  }



#endif
} /* void rho_send_recv_begin( lattice_ptr lattice, const int subs) */


void rho_send_recv_end( lattice_ptr lattice, const int subs)
{
#if PARALLEL
  int n;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  int ip, in;
  int jp, jn;
  int mpierr;

  mpierr = MPI_Wait(
      /* mpi_request *req */&(lattice->process.send_req_2),
      /* mpi_status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
      /* mpi_request *req */&(lattice->process.recv_req_2),
      /* mpi_status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
      /* mpi_request *req */&(lattice->process.send_req_3),
      /* mpi_status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
      /* mpi_request *req */&(lattice->process.recv_req_3),
      /* mpi_status *stat */&(lattice->process.mpi_status));

#endif
} /* void rho_send_recv_end( lattice_ptr lattice, const int subs) */
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void solid_send_recv_begin( lattice_ptr lattice, const int subs)
{
#if PARALLEL
  int n;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  int mpierr;

  //rev Huang
  //To send "rho" and "solid" to neighbour blocks
  n =0;
  k = get_LZ(lattice)-1;
  for( j=0; j<nj; j++)
  {
    for( i=0; i<ni; i++)
    {
#if SAVE_MEMO

#else
      lattice->process.pos_dir_solid_to_send[n] =
	lattice->solids[subs][ IJK2N( i , j , k, ni, nj)].is_solid;
      lattice->process.neg_dir_solid_to_send[n] =
	lattice->solids[subs][ IJK2N( i , j , 0, ni, nj)].is_solid;
      n++;
#endif
    } /* if( i=0; i<ni; i++) */
  } /* if( j=0; j<nj; j++) */


  // SOLID S E N D   I N   P O S I T I V E   D I R E C T I O N
  //#########################################################################
#if VERBOSITY_LEVEL > 1
  printf( "%s %d %04d >> "
      "MPI_Isend( %04d)"
      "\n",
      __FILE__,__LINE__,get_proc_id(lattice),
      (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
#endif
  mpierr =
    MPI_Isend(
	/*void *buf*/          lattice->process.pos_dir_solid_to_send,
	/*int count*/          sol_buffer_size,
	/*MPI_Datatype dtype*/ MPI_INT,
	/*int dest*/         ( get_proc_id(lattice)
	  + get_num_procs(lattice)+1)
	% get_num_procs(lattice),
	/*int tag*/            4,
	/*MPI_Comm comm*/      MPI_COMM_WORLD,
	/*MPI_Request *req*/   &(lattice->process.send_req_4)
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Isend( %04d)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice),
	mpierr,
	(get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
    process_finalize();
    process_exit(1);
  }
  // SOLID R E C V   F R O M   N E G A T I V E   D I R E C T I O N
  //#########################################################################
#if VERBOSITY_LEVEL > 1
  printf( "%s %d %04d >> "
      "MPI_Irecv( %04d)"
      "\n",
      __FILE__,__LINE__,get_proc_id(lattice),
      (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
#endif
  mpierr =
    MPI_Irecv(
	/*void *buf*/          lattice->process.pos_dir_solid_to_recv,
	/*int count*/          sol_buffer_size,
	/*MPI_Datatype dtype*/ MPI_INT,
	/*int src*/          ( get_proc_id(lattice)
	  + get_num_procs(lattice)-1)
	% get_num_procs(lattice),
	/*int tag*/            4,
	/*MPI_Comm comm*/      MPI_COMM_WORLD,
	/*MPI_Request *req*/   &(lattice->process.recv_req_4)
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Irecv( %04d)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice),
	mpierr,
	(get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
    process_finalize();
    process_exit(1);
  }
  // SOLID S E N D   I N   N E G A T I V E   D I R E C T I O N
  //#########################################################################
#if VERBOSITY_LEVEL > 1
  printf( "%s %d %04d >> "
      "MPI_Isend( %04d)"
      "\n",
      __FILE__,__LINE__,get_proc_id(lattice),
      (get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
#endif
  mpierr =
    MPI_Isend(
	/*void *buf*/          lattice->process.neg_dir_solid_to_send,
	/*int count*/          sol_buffer_size,
	/*MPI_Datatype dtype*/ MPI_INT,
	/*int dest*/         ( get_proc_id(lattice)
	  + get_num_procs(lattice)-1)
	% get_num_procs(lattice),
	/*int tag*/            5,
	/*MPI_Comm comm*/      MPI_COMM_WORLD,
	/*MPI_Request *req*/   &(lattice->process.send_req_5)
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Isend( %04d)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice),
	mpierr,
	(get_proc_id(lattice)+get_num_procs(lattice)-1)%get_num_procs(lattice));
    process_finalize();
    process_exit(1);
  }
  // SOLID R E C V   F R O M   P O S I T I V E   D I R E C T I O N
  //#########################################################################
#if VERBOSITY_LEVEL > 1
  printf( "%s %d %04d >> "
      "MPI_Irecv( %04d)"
      "\n",
      __FILE__,__LINE__,get_proc_id(lattice),
      (get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
#endif
  mpierr =
    MPI_Irecv(
	/*void *buf*/          lattice->process.neg_dir_solid_to_recv,
	/*int count*/          sol_buffer_size,
	/*MPI_Datatype dtype*/ MPI_INT,
	/*int src*/          ( get_proc_id(lattice)
	  + get_num_procs(lattice)+1)
	% get_num_procs(lattice),
	/*int tag*/            5,
	/*MPI_Comm comm*/      MPI_COMM_WORLD,
	/*MPI_Request *req*/   &(lattice->process.recv_req_5)
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Irecv( %04d)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice),
	mpierr,
	(get_proc_id(lattice)+get_num_procs(lattice)+1)%get_num_procs(lattice));
    process_finalize();
    process_exit(1);
  }



#endif
} /* void rho_send_recv_begin( lattice_ptr lattice, const int subs) */


void solid_send_recv_end( lattice_ptr lattice, const int subs)
{
#if PARALLEL
  int n;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  int ip, in;
  int jp, jn;
  int mpierr;

  mpierr = MPI_Wait(
      /* MPI_Request *req */&(lattice->process.send_req_4),
      /* MPI_Status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
      /* MPI_Request *req */&(lattice->process.recv_req_4),
      /* MPI_Status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
      /* MPI_Request *req */&(lattice->process.send_req_5),
      /* MPI_Status *stat */&(lattice->process.mpi_status));
  mpierr = MPI_Wait(
      /* MPI_Request *req */&(lattice->process.recv_req_5),
      /* MPI_Status *stat */&(lattice->process.mpi_status));


#endif
} /* void rho_send_recv_end( lattice_ptr lattice, const int subs) */
//#####################################################################
//#####################################################################




//##############################################################################
void process_dump_pdfs_to_recv( lattice_ptr lattice, char *comment_str)
{
#if PARALLEL
  char new_comment[1024];
  sprintf( new_comment, "pos pdfs to recv, %s", comment_str);
  process_dump_pdfs(
      lattice,
      new_comment,
      lattice->process.pos_dir_pdf_to_recv);
  sprintf( new_comment, "neg pdfs to recv, %s", comment_str);
  process_dump_pdfs(
      lattice,
      new_comment,
      lattice->process.neg_dir_pdf_to_recv);
#endif
} /* void process_dump_pdfs_to_recv( lattice_ptr lattice)  */

//##############################################################################
void process_dump_pdfs_to_send( lattice_ptr lattice, char *comment_str)
{
#if PARALLEL
  char new_comment[1024];
  sprintf( new_comment, "pos pdfs to send, %s", comment_str);
  process_dump_pdfs(
      lattice,
      new_comment,
      lattice->process.pos_dir_pdf_to_send);
  sprintf( new_comment, "neg pdfs to send, %s", comment_str);
  process_dump_pdfs(
      lattice,
      new_comment,
      lattice->process.neg_dir_pdf_to_send);
#endif
} /* void process_dump_pdfs_to_send( lattice_ptr lattice)  */

//##############################################################################
void process_dump_pdfs( lattice_ptr lattice, char *comment_str, real *pdfs)
{
#if PARALLEL
  int n, p;
  int i, j, k;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  int mpierr;

  for( p=0; p<get_num_procs( lattice); p++)
  {
    MPI_Barrier( MPI_COMM_WORLD);

    if( p == get_proc_id(lattice))
    {
      printf("\n\n// Proc %d, \"%s\".", get_proc_id( lattice), comment_str);
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
    } /* if( p == get_proc_id(lattice)) */
  } /* for( p=0; p<get_num_procs( lattice); p++) */

  MPI_Barrier( MPI_COMM_WORLD);

#endif
} /* void process_dump_pdfs_to_recv( lattice_ptr lattice)  */

//##############################################################################
void gather_north_pointing_pdfs(
    lattice_ptr lattice,
    real *north,
    const int subs,
    const int k,
    const int which_pdf)
{
#if PARALLEL
  int i, j, n;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);

  switch(which_pdf)
  {
    case 0:
      n = 0;
      for( j=0; j<nj; j++)
      {
	for( i=0; i<ni; i++)
	{
#if SAVE_MEMO
#else
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].feq[ T]; n++;
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].feq[TW]; n++;
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].feq[TE]; n++;
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].feq[TN]; n++;
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].feq[TS]; n++;
#endif
	} /* if( i=0; i<ni; i++) */
      } /* if( j=0; j<nj; j++) */
      break;
    case 1:
      n = 0;
      for( j=0; j<nj; j++)
      {
	for( i=0; i<ni; i++)
	{
#if SAVE_MEMO
#else
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].f[ T]; n++;
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].f[TW]; n++;
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].f[TE]; n++;
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].f[TN]; n++;
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].f[TS]; n++;
#endif
	} /* if( i=0; i<ni; i++) */
      } /* if( j=0; j<nj; j++) */
      break;
    case 2:
      n = 0;
      for( j=0; j<nj; j++)
      {
	for( i=0; i<ni; i++)
	{
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].ftemp[ T]; n++;
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].ftemp[TW]; n++;
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].ftemp[TE]; n++;
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].ftemp[TN]; n++;
	  north[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].ftemp[TS]; n++;

	} /* if( i=0; i<ni; i++) */
      } /* if( j=0; j<nj; j++) */
      break;
    default:
      printf("%s %d %04d >> ERROR: Unhandled case which_pdf=%d. Exiting!",
	  __FILE__,__LINE__,get_proc_id(lattice), which_pdf);
      exit(1);
      break;
  } /* switch(which_pdf) */
#endif
} /* void gather_north_pointing_pdfs( lattice_ptr lattice, real *north) */

//##############################################################################
void gather_south_pointing_pdfs(
    lattice_ptr lattice,
    real *south,
    const int subs,
    const int k,
    const int which_pdf)
{
#if PARALLEL
  int i, j, n;
  int ni = get_LX( lattice),
      nj = get_LY( lattice);
  switch(which_pdf)
  {
    case 0:
      n = 0;
      for( j=0; j<nj; j++)
      {
	for( i=0; i<ni; i++)
	{
#if SAVE_MEMO
#else
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].feq[ B]; n++;
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].feq[BW]; n++;
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].feq[BE]; n++;
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].feq[BN]; n++;
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].feq[BS]; n++;
#endif
	} /* if( i=0; i<ni; i++) */
      } /* if( j=0; j<nj; j++) */
      break;
    case 1:
      n = 0;
      for( j=0; j<nj; j++)
      {
	for( i=0; i<ni; i++)
	{
#if SAVE_MEMO
#else
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].f[ B]; n++;
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].f[BW]; n++;
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].f[BE]; n++;
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].f[BN]; n++;
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].f[BS]; n++;
#endif
	} /* if( i=0; i<ni; i++) */
      } /* if( j=0; j<nj; j++) */
      break;
    case 2:
      n = 0;
      for( j=0; j<nj; j++)
      {
	for( i=0; i<ni; i++)
	{
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].ftemp[ B]; n++;
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].ftemp[BW]; n++;
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].ftemp[BE]; n++;
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].ftemp[BN]; n++;
	  south[n] = lattice->pdf[subs][ IJK2N(i,j,k, ni, nj)].ftemp[BS]; n++;

	} /* if( i=0; i<ni; i++) */
      } /* if( j=0; j<nj; j++) */
      break;
    default:
      printf("%s %d %04d >> ERROR: Unhandled case which_pdf=%d. Exiting!",
	  __FILE__,__LINE__,get_proc_id(lattice), which_pdf);
      exit(1);
      break;
  }
#endif
} /* void gather_south_pointing_pdfs( lattice_ptr lattice, real *south) */

void process_reduce_real_sum( lattice_ptr lattice, real *arg_x)
{
#if PARALLEL
  real sum_x;
  int mpierr;

  //
  // INPUT PARAMETERS
  //        sbuf   - address of send buffer (choice)
  //        count  - number of elements in send buffer (integer)
  //        dtype  - data type of elements of send buffer (handle)
  //        op     - reduce operation (handle)
  //        root   - rank of root process (integer)
  //        comm   - communicator (handle)
  //
  // OUTPUT PARAMETER
  //        rbuf   - address of receive buffer (choice, sig't only at root )
  //
  mpierr =
    MPI_Reduce(
	/*void *sbuf*/         arg_x,
	/*void* rbuf*/        &sum_x,
	/*int count*/          1,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*MPI_Op op*/          MPI_SUM,
	/*int root*/           0,
	/*MPI_Comm comm*/      MPI_COMM_WORLD
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Reduce( ave_rho, MPI_DOUBLE, MPI_SUM)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice), mpierr);
    process_finalize();
    exit(1);
  }
  if( is_on_root_proc( lattice))
  {
    *arg_x = sum_x;
  }
#endif
} /* void process_reduce_real_sum( lattice_ptr lattice, real &arg_x) */

void process_allreduce_real_sum( lattice_ptr lattice, real *arg_x)
{
#if PARALLEL
  real sum_x;
  int mpierr;

  //
  // INPUT PARAMETERS
  //        sbuf   - address of send buffer (choice)
  //        count  - number of elements in send buffer (integer)
  //        dtype  - data type of elements of send buffer (handle)
  //        op     - reduce operation (handle)
  //        root   - rank of root process (integer)
  //        comm   - communicator (handle)
  //
  // OUTPUT PARAMETER
  //        rbuf   - address of receive buffer (choice, sig't only at root )
  //
  mpierr =
    MPI_Allreduce(
	/*void *sbuf*/         arg_x,
	/*void* rbuf*/        &sum_x,
	/*int count*/          1,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*MPI_Op op*/          MPI_SUM,
	/*MPI_Comm comm*/      MPI_COMM_WORLD
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Reduce( ave_rho, MPI_DOUBLE, MPI_SUM)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice), mpierr);
    process_finalize();
    exit(1);
  }

  *arg_x = sum_x;

#endif
} /* void process_reduce_real_sum( lattice_ptr lattice, real &arg_x) */

void process_reduce_int_sum( lattice_ptr lattice, int *arg_n)
{
#if PARALLEL
  real sum_n;
  int mpierr;

  //
  // INPUT PARAMETERS
  //        sbuf   - address of send buffer (choice)
  //        count  - number of elements in send buffer (integer)
  //        dtype  - data type of elements of send buffer (handle)
  //        op     - reduce operation (handle)
  //        root   - rank of root process (integer)
  //        comm   - communicator (handle)
  //
  // OUTPUT PARAMETER
  //        rbuf   - address of receive buffer (choice, sig't only at root )
  //
  mpierr =
    MPI_Reduce(
	/*void *sbuf*/         arg_n,
	/*void* rbuf*/        &sum_n,
	/*int count*/          1,
	/*MPI_Datatype dtype*/ MPI_INT,
	/*MPI_Op op*/          MPI_SUM,
	/*int root*/           0,
	/*MPI_Comm comm*/      MPI_COMM_WORLD
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Reduce( ave_rho, MPI_DOUBLE, MPI_SUM)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice), mpierr);
    process_finalize();
    exit(1);
  }
  if( is_on_root_proc( lattice))
  {
    *arg_n = sum_n;
  }
#endif
} /* void process_reduce_int_sum( lattice_ptr lattice, int *arg_n) */

void process_allreduce_int_sum( lattice_ptr lattice, int *arg_n)
{
#if PARALLEL
  real sum_n;
  int mpierr;

  //
  // INPUT PARAMETERS
  //        sbuf   - address of send buffer (choice)
  //        count  - number of elements in send buffer (integer)
  //        dtype  - data type of elements of send buffer (handle)
  //        op     - reduce operation (handle)
  //        root   - rank of root process (integer)
  //        comm   - communicator (handle)
  //
  // OUTPUT PARAMETER
  //        rbuf   - address of receive buffer (choice, sig't only at root )
  //
  mpierr =
    MPI_Allreduce(
	/*void *sbuf*/         arg_n,
	/*void* rbuf*/        &sum_n,
	/*int count*/          1,
	/*MPI_Datatype dtype*/ MPI_INT,
	/*MPI_Op op*/          MPI_SUM,
	/*MPI_Comm comm*/      MPI_COMM_WORLD
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Reduce( ave_rho, MPI_DOUBLE, MPI_SUM)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice), mpierr);
    process_finalize();
    exit(1);
  }

  *arg_n = sum_n;

#endif
} /* void process_reduce_int_sum( lattice_ptr lattice, int *arg_n) */

void process_reduce_real_max( lattice_ptr lattice, real *arg_x)
{
#if PARALLEL
  real max_x;
  int mpierr;

  //
  // INPUT PARAMETERS
  //        sbuf   - address of send buffer (choice)
  //        count  - number of elements in send buffer (integer)
  //        dtype  - data type of elements of send buffer (handle)
  //        op     - reduce operation (handle)
  //        root   - rank of root process (integer)
  //        comm   - communicator (handle)
  //
  // OUTPUT PARAMETER
  //        rbuf   - address of receive buffer (choice, sig't only at root )
  //
  mpierr =
    MPI_Reduce(
	/*void *sbuf*/         arg_x,
	/*void* rbuf*/        &max_x,
	/*int count*/          1,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*MPI_Op op*/          MPI_MAX,
	/*int root*/           0,
	/*MPI_Comm comm*/      MPI_COMM_WORLD
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Reduce( ave_rho, MPI_DOUBLE, MPI_SUM)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice), mpierr);
    process_finalize();
    exit(1);
  }
  if( is_on_root_proc( lattice))
  {
    *arg_x = max_x;
  }
#endif
} /* void process_reduce_real_sum( lattice_ptr lattice, real &arg_x) */

void process_allreduce_real_max( lattice_ptr lattice, real *arg_x)
{
#if PARALLEL
  real max_x;
  int mpierr;

  //
  // INPUT PARAMETERS
  //        sbuf   - address of send buffer (choice)
  //        count  - number of elements in send buffer (integer)
  //        dtype  - data type of elements of send buffer (handle)
  //        op     - reduce operation (handle)
  //        comm   - communicator (handle)
  //
  // OUTPUT PARAMETER
  //        rbuf   - address of receive buffer (choice, sig't only at root )
  //
  mpierr =
    MPI_Allreduce(
	/*void *sbuf*/         arg_x,
	/*void* rbuf*/        &max_x,
	/*int count*/          1,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*MPI_Op op*/          MPI_MAX,
	/*MPI_Comm comm*/      MPI_COMM_WORLD
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Allreduce( ave_rho, MPI_DOUBLE, MPI_SUM)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice), mpierr);
    process_finalize();
    exit(1);
  }

  *arg_x = max_x;

#endif
} /* void process_reduce_real_sum( lattice_ptr lattice, real &arg_x) */

void process_reduce_real_min( lattice_ptr lattice, real *arg_x)
{
#if PARALLEL
  real min_x;
  int mpierr;

  //
  // INPUT PARAMETERS
  //        sbuf   - address of send buffer (choice)
  //        count  - number of elements in send buffer (integer)
  //        dtype  - data type of elements of send buffer (handle)
  //        op     - reduce operation (handle)
  //        root   - rank of root process (integer)
  //        comm   - communicator (handle)
  //
  // OUTPUT PARAMETER
  //        rbuf   - address of receive buffer (choice, sig't only at root )
  //
  mpierr =
    MPI_Reduce(
	/*void *sbuf*/         arg_x,
	/*void* rbuf*/        &min_x,
	/*int count*/          1,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*MPI_Op op*/          MPI_MIN,
	/*int root*/           0,
	/*MPI_Comm comm*/      MPI_COMM_WORLD
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Reduce( ave_rho, MPI_DOUBLE, MPI_SUM)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice), mpierr);
    process_finalize();
    exit(1);
  }
  if( is_on_root_proc( lattice))
  {
    *arg_x = min_x;
  }
#endif
} /* void process_reduce_real_sum( lattice_ptr lattice, real &arg_x) */

void process_allreduce_real_min( lattice_ptr lattice, real *arg_x)
{
#if PARALLEL
  real min_x;
  int mpierr;

  //
  // INPUT PARAMETERS
  //        sbuf   - address of send buffer (choice)
  //        count  - number of elements in send buffer (integer)
  //        dtype  - data type of elements of send buffer (handle)
  //        op     - reduce operation (handle)
  //        comm   - communicator (handle)
  //
  // OUTPUT PARAMETER
  //        rbuf   - address of receive buffer (choice, sig't only at root )
  //
  mpierr =
    MPI_Allreduce(
	/*void *sbuf*/         arg_x,
	/*void* rbuf*/        &min_x,
	/*int count*/          1,
	/*MPI_Datatype dtype*/ MPI_DOUBLE,
	/*MPI_Op op*/          MPI_MIN,
	/*MPI_Comm comm*/      MPI_COMM_WORLD
	);
  if( mpierr != MPI_SUCCESS)
  {
    printf( "%s %d %04d >> "
	"ERROR: %d <-- MPI_Allreduce( ave_rho, MPI_DOUBLE, MPI_SUM)"
	"\n",
	__FILE__,__LINE__,get_proc_id(lattice), mpierr);
    process_finalize();
    exit(1);
  }

  *arg_x = min_x;

#endif
} /* void process_reduce_real_sum( lattice_ptr lattice, real &arg_x) */

//##############################################################################
void process_barrier()
{
#if PARALLEL
  MPI_Barrier( MPI_COMM_WORLD);
#endif
}

//##############################################################################
void process_tic( lattice_ptr lattice)
{
#if PARALLEL
  set_tic( lattice, MPI_Wtime());
#else
  set_tic( lattice, ((real)clock())/(real)CLK_TCK);
#endif
}

//##############################################################################
void process_toc( lattice_ptr lattice)
{
#if PARALLEL
  set_toc( lattice, MPI_Wtime());
#else
  set_toc( lattice, ((real)clock())/(real)CLK_TCK);
#endif
}

//##############################################################################
void process_finalize()
{
#if PARALLEL
  MPI_Finalize();
#endif
}

//##############################################################################
void process_exit( int exit_val)
{
  process_finalize();
  exit(exit_val);
}

//##############################################################################
//
// Accessor methods for the process struct.
//

int get_proc_id( lattice_ptr lattice) { return lattice->process.id;}
int get_num_procs( lattice_ptr lattice) { return lattice->process.num_procs;}
int is_on_root_proc( lattice_ptr lattice) { return !(lattice->process.id);}

#if PARALLEL
int get_g_LX( lattice_ptr lattice) { return lattice->process.g_LX;}
int get_g_LY( lattice_ptr lattice) { return lattice->process.g_LY;}
int get_g_LZ( lattice_ptr lattice) { return lattice->process.g_LZ;}
int get_g_SX( lattice_ptr lattice) { return lattice->process.g_SX;}
int get_g_SY( lattice_ptr lattice) { return lattice->process.g_SY;}
int get_g_SZ( lattice_ptr lattice) { return lattice->process.g_SZ;}
int get_g_EX( lattice_ptr lattice) { return lattice->process.g_EX;}
int get_g_EY( lattice_ptr lattice) { return lattice->process.g_EY;}
int get_g_EZ( lattice_ptr lattice) { return lattice->process.g_EZ;}

void set_g_LX( lattice_ptr lattice, const int arg_LX)
{
  lattice->process.g_LX = arg_LX;
}
void set_g_LY( lattice_ptr lattice, const int arg_LY)
{
  lattice->process.g_LY = arg_LY;
}
void set_g_LZ( lattice_ptr lattice, const int arg_LZ)
{
  lattice->process.g_LZ = arg_LZ;
}

void set_g_SX( lattice_ptr lattice, const int arg_SX)
{
  lattice->process.g_SX = arg_SX;
}
void set_g_SY( lattice_ptr lattice, const int arg_SY)
{
  lattice->process.g_SY = arg_SY;
}
void set_g_SZ( lattice_ptr lattice, const int arg_SZ)
{
  lattice->process.g_SZ = arg_SZ;
}

void set_g_EX( lattice_ptr lattice, const int arg_EX)
{
  lattice->process.g_EX = arg_EX;
}
void set_g_EY( lattice_ptr lattice, const int arg_EY)
{
  lattice->process.g_EY = arg_EY;
}
void set_g_EZ( lattice_ptr lattice, const int arg_EZ)
{
  lattice->process.g_EZ = arg_EZ;
}

int get_g_NumNodes( lattice_ptr lattice) { return lattice->process.g_NumNodes;}
void set_g_NumNodes( lattice_ptr lattice, const int arg_NumNodes)
{
  lattice->process.g_NumNodes = arg_NumNodes;
}
void set_g_StartNode( lattice_ptr lattice, const int arg_n)
{
  lattice->process.g_StartNode = arg_n;
}
int get_g_StartNode( lattice_ptr lattice)
{
  return lattice->process.g_StartNode;
}

#else
// Defaults for non-parallel runs.
int get_g_LX( lattice_ptr lattice) { return get_LX( lattice);}
int get_g_LY( lattice_ptr lattice) { return get_LY( lattice);}
int get_g_LZ( lattice_ptr lattice) { return get_LZ( lattice);}
int get_g_SX( lattice_ptr lattice) { return 0;}
int get_g_SY( lattice_ptr lattice) { return 0;}
int get_g_SZ( lattice_ptr lattice) { return 0;}
int get_g_EX( lattice_ptr lattice) { return get_LX( lattice)-1;}
int get_g_EY( lattice_ptr lattice) { return get_LY( lattice)-1;}
int get_g_EZ( lattice_ptr lattice) { return get_LZ( lattice)-1;}
int get_g_NumNodes( lattice_ptr lattice) { return get_NumNodes( lattice);}
int get_g_StartNode( lattice_ptr lattice) { return 0;}
#endif

