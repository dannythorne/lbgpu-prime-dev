//##############################################################################
//
// lbio.c
//
//  - Lattice Boltzmann I/O routines.
//
//  - Mainly, dump the data to files that can be read by Matlab.
//
//  - Also, output routines for facilitating debugging.
//
//  - Should have a routine that dumps a matlab script?
//

// Some compilers, e.g., VC++, don't have the usual round() function
// in their math library.  Alternatively, ROUND can be defined as
// ceil or floor or some other rounding function.  It is used in
// the below routines for converting the real number value of
// quantities at a lattice node into integer RGB values for writing
// to BMP files.
#define ROUND floor

//#if SWAP_BYTE_ORDER || OSTYPE==darwin
#if SWAP_BYTE_ORDER
// Swap byte order.
#define ENDIAN2(w) ((((w)&0x00ff)<<8)|(((w)&0xff00)>>8))
#define ENDIAN4(w) ((((w)&0x000000ff)<<24)|(((w)&0xff000000)>>24)|(((w)&0x0000ff00)<<8)|(((w)&0x00ff0000)>>8))
#else /* !( SWAP_BYTE_ORDER) */
#define ENDIAN2(w) (w)
#define ENDIAN4(w) (w)
#endif /* SWAP_BYTE_ORDER */

//void output_frame( lattice_ptr lattice)
//##############################################################################
//
// O U T P U T   F R A M E
//
void output_frame( lattice_ptr lattice)
{
  real s, u[3];
  real nu;
  real L;
  int    n, subs;

#if VERBOSITY_LEVEL > 0
  printf("\n");
  printf( "========================================"
          "========================================\n");
  printf("Begin file I/O at time = %d, frame = %d.\n",
      lattice->time, lattice->time/lattice->param.FrameRate);
  printf("\n");
#endif /* VERBOSITY_LEVEL > 0 */

  for( subs = 0; subs < NUM_FLUID_COMPONENTS; subs++)
  {
    compute_max_u( lattice, subs);
    compute_min_u( lattice, subs);
    compute_ave_u( lattice, subs);
    compute_max_rho( lattice, subs);
    compute_min_rho( lattice, subs);
    compute_ave_rho( lattice, subs);
    compute_flux( lattice, subs);
  }
  dump_frame_summary( lattice);

#if WRITE_MACRO_VAR_DAT_FILES || WRITE_MACRO_VAR_RAW_FILES  || WRITE_PLOT_FILE
  dump_macro_vars( lattice, lattice->time);
#endif /* WRITE_MACRO_VAR_DAT_FILES || WRITE_MACRO_VAR_RAW_FILES */

  if( write_debug_txt_files(lattice))
  {
    write_rho_txt(lattice);
  }

#if 0
#if WRITE_PDF_DAT_FILES
  dump_pdf( lattice, lattice->time);
#endif /* WRITE_PDF_DAT_FILES */

  slice( lattice);
#endif

#if VERBOSITY_LEVEL > 0
  printf("\n");
  printf("File I/O done.\n");
  printf("--\n");
#endif /* VERBOSITY_LEVEL > 0 */

#if 0
  nu = (1./3.)*(lattice->param.tau[0] - .5);
  L = lattice->param.length_scale;

  compute_ave_u( lattice, u, 0);
  s = sqrt( u[0]*u[0] + u[1]*u[1]);
  printf("subs 0: Re = ux_ave*L/nu = %f * %f / %f = %f\n",
    u[0], L, nu, u[0]*L/nu );
  printf("subs 0: Re = uy_ave*L/nu = %f * %f / %f = %f\n",
    u[1], L, nu, u[1]*L/nu );
  printf("subs 0: Re = u_ave*L/nu  = %f * %f / %f = %f\n",
    s, L, nu, s*L/nu );

#if NUM_FLUID_COMPONENTS == 2
  compute_ave_u( lattice, u, 1);
  s = sqrt( u[0]*u[0] + u[1]*u[1]);
  printf("subs 1: Re = ux_ave*L/nu = %f * %f / %f = %f\n",
    u[0], L, nu, u[0]*L/nu );
  printf("subs 1: Re = uy_ave*L/nu = %f * %f / %f = %f\n",
    u[1], L, nu, u[1]*L/nu );
  printf("subs 1: Re = u_ave*L/nu  = %f * %f / %f = %f\n",
    s, L, nu, s*L/nu );
#endif /* NUM_FLUID_COMPONENTS == 2 */

#if STORE_UEQ
  compute_ave_ueq( lattice, u);
  s = sqrt( u[0]*u[0] + u[1]*u[1]);
  printf("eq:     Re = ux_ave*L/nu = %f * %f / %f = %f\n",
    u[0], L, nu, u[0]*L/nu );
  printf("eq:     Re = uy_ave*L/nu = %f * %f / %f = %f\n",
    u[1], L, nu, u[1]*L/nu );
  printf("eq:     Re = u_ave*L/nu  = %f * %f / %f = %f\n",
    s, L, nu, s*L/nu );
#endif /* STORE_UEQ */
#endif

} /* void output_frame( lattice_ptr lattice) */

// void dump_frame_info( struct lattice_struct *lattice)
//##############################################################################
//
// D U M P   F R A M E   I N F O
//
void dump_frame_summary( struct lattice_struct *lattice)
{
  char   filename[1024];
  FILE   *o;
  real min_u[3], max_u[3],  ave_u[3];
  real flux[4];
  real min_rho, max_rho,   ave_rho;
  real rho_ratio, u_x_ratio, u_y_ratio, u_z_ratio;
  int    subs;

 for( subs = 0; subs < NUM_FLUID_COMPONENTS; subs++)
 {
  if( 1 || is_on_root_proc( lattice))
  {
    gen_filename( lattice, filename, "frames", -1, subs, ".dat");

  // On the first timestep, make sure we start with a new file.
  if( lattice->time==0)
  {
      if( !( o = fopen(filename,"w+")))
      {
        printf("ERROR: fopen(\"%s\",\"w+\") = NULL.  Bye, bye!\n", filename);
        process_exit(1);
      }
      else
      {
        // Put a header on the file.
        fprintf( o, "\n");
        fprintf( o,
        "        time "
        "         |j| "
        "         j_x "
        "         j_y "
        "         j_z "
        "    min_u[x] "
        "    min_u[y] "
        "    min_u[z] "
        "    max_u[x] "
        "    max_u[y] "
        "    max_u[z] "
        "    ave_u[x] "
        "    ave_u[y] "
        "    ave_u[z] "
        "  max/ave[x] "
        "  max/ave[y] "
        "  max/ave[z] "
        "     min_rho "
        "     max_rho "
        "     ave_rho "
        "     max/ave "
        "\n");fprintf( o,
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "------------ "
        "\n");
        fclose(o);
      }
  }
  if( !( o = fopen(filename,"a+")))
  {
    printf("ERROR: fopen(\"%s\",\"a+\") = NULL.  Bye, bye!\n", filename);
    process_exit(1);
  }
  }

  max_u[0] = get_max_ux( lattice, subs);
  min_u[0] = get_min_ux( lattice, subs);
  ave_u[0] = get_ave_ux( lattice, subs);
  max_u[1] = get_max_uy( lattice, subs);
  min_u[1] = get_min_uy( lattice, subs);
  ave_u[1] = get_ave_uy( lattice, subs);
  max_u[2] = get_max_uz( lattice, subs);
  min_u[2] = get_min_uz( lattice, subs);
  ave_u[2] = get_ave_uz( lattice, subs);
  max_rho = get_max_rho( lattice, subs);
  min_rho = get_min_rho( lattice, subs);
  ave_rho = get_ave_rho( lattice, subs);

  flux[0] = get_flux_u( lattice, subs);
  flux[1] = get_flux_x( lattice, subs);
  flux[2] = get_flux_y( lattice, subs);
  flux[3] = get_flux_z( lattice, subs);


  if( 1|| is_on_root_proc( lattice))
  {
    rho_ratio = ( ave_rho  != 0.) ? ( max_rho /ave_rho ):( 1.);
    u_x_ratio = ( ave_u[0] != 0.) ? ( max_u[0]/ave_u[0]):( 1.);
    u_y_ratio = ( ave_u[1] != 0.) ? ( max_u[1]/ave_u[1]):( 1.);
    u_z_ratio = ( ave_u[2] != 0.) ? ( max_u[2]/ave_u[2]):( 1.);

    fprintf( o,
      "%12d "
      "%12.7f %12.7f %12.7f %12.7f"
      "%12.7f %12.7f %12.7f "
      "%12.7f %12.7f %12.7f "
      "%12.7f %12.7f %12.7f "
      "%12.7f %12.7f %12.7f "
      "%12.7f %12.7f %12.7f "
      "%12.7f\n",
      lattice->time,
      flux[0], flux[1], flux[2], flux[3],
      min_u[0], min_u[1], min_u[2],
      max_u[0], max_u[1], max_u[2],
      ave_u[0], ave_u[1], ave_u[2],
      (fabs(u_x_ratio)<=9999.)?(u_x_ratio):(1./0.),
      (fabs(u_y_ratio)<=9999.)?(u_y_ratio):(1./0.),
      (fabs(u_z_ratio)<=9999.)?(u_z_ratio):(1./0.),
      min_rho,
      max_rho,
      ave_rho,
      (fabs(rho_ratio)<=9999.)?(rho_ratio):(1./0.) );

  fclose(o);

#if VERBOSITY_LEVEL > 0
    printf("%s %d %04d >> dump_frame_info() -- "
      "Wrote/appended to file \"%s\"\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);
#endif /* VERBOSITY_LEVEL > 0 */

#if VERBOSITY_LEVEL > 0
    printf("%s %d %04d >> dump_frame_info() -- "
      "frame = %d/%d = %d\n",
      __FILE__, __LINE__, get_proc_id(lattice),
      lattice->time,
      lattice->param.FrameRate,
      (int)((real)lattice->time/(real)lattice->param.FrameRate));
#endif /* VERBOSITY_LEVEL > 0 */
  }
 }

} /* void dump_frame_summary( struct lattice_struct *lattice) */

// void dump_macro_vars( struct lattice_struct *lattice)
//##############################################################################
//
// D U M P   M A C R O S C O P I C
//
//  - Output the macro_vars variables to files.
//
void dump_macro_vars( struct lattice_struct *lattice, int time)
{
  char   filename[1024],filename2[1024], fn3[1024];
  FILE   *o, *o_u, *o_rho, *o_ux, *o_uy, *o_ueq, *o_ueq_x, *o_ueq_y, *sat;
  int    *node_ptr;
  int    i, j, n, count, count1, count2, basek;
  real *macro_vars_ptr;
  real *ueq;
  int    frame;
  real min_u[3], max_u[3],  ave_u[3];
  real min_rho, max_rho,   ave_rho;
  real rho_ratio, u_x_ratio, u_y_ratio;
  int    subs;
  int    ni, nj, nk;
  int    j_slice, k_slice;

  frame = get_frame(lattice);

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
    max_rho = get_max_rho( lattice, subs);
    min_rho = get_min_rho( lattice, subs);
    ave_rho = get_ave_rho( lattice, subs);

   if( is_on_root_proc( lattice))
    {
    printf("%s %d %04d >> "
      "min_rho = %20.17f, max_rho = %20.17f, ave_rho = %20.17f\n",
      __FILE__,__LINE__, get_proc_id(lattice),
      min_rho, max_rho, ave_rho);
   } /* if( is_on_root_proc( lattice)) */

    write_rho_image( lattice, subs);
    write_u_image( lattice, subs);
    write_ux_image( lattice, subs);
    write_uy_image( lattice, subs);
    write_uz_image( lattice, subs);

    write_rho_dat( lattice, subs);
    write_ux_dat( lattice, subs);
    write_uy_dat( lattice, subs);
    write_uz_dat( lattice, subs);


    max_u[0] = get_max_ux( lattice, subs);
    min_u[0] = get_min_ux( lattice, subs);
    ave_u[0] = get_ave_ux( lattice, subs);
    max_u[1] = get_max_uy( lattice, subs);
    min_u[1] = get_min_uy( lattice, subs);
    ave_u[1] = get_ave_uy( lattice, subs);
    max_u[2] = get_max_uz( lattice, subs);
    min_u[2] = get_min_uz( lattice, subs);
    ave_u[2] = get_ave_uz( lattice, subs);

    if( is_on_root_proc( lattice))
    {
      printf("%s %d %04d >> "
        "min_u = [ %f %f %f], max_u = [ %f %f %f], ave_u = [ %f %f %f]\n",
        __FILE__,__LINE__, get_proc_id(lattice),
      min_u[0], min_u[1], min_u[2],
      max_u[0], max_u[1], max_u[2],
      ave_u[0], ave_u[1], ave_u[2] );
    } /* if( is_on_root_proc( lattice)) */

//------------------------------------------------------------------------------------------
//
//  Write Vmag files for post processing large domains that have memory issues
//
//  edit: peter
//
#if WRITE_VMAG

int eye;
real vx, vy, vz, vxsq, vysq, vzsq;
real mag=0;
real vmagmax=0, vmagmin=0;
real vmag[ni*nj*nk];
real vmagscale[ni*nj*nk];

  for(eye=0; eye<(ni*nj*nk);eye++)
  {
    vx=lattice->macro_vars[subs][eye].u[0];
    vxsq = vx*vx;
    vy=lattice->macro_vars[subs][eye].u[1];
    vysq = vy*vy;
    vz=lattice->macro_vars[subs][eye].u[2];
    vzsq = vz*vz;
    mag = sqrt(vxsq+vysq+vzsq);
    vmag[eye] = mag;
    if(vmag[eye]>vmagmax)
    {
      vmagmax = vmag[eye];
    }
    if(vmag[eye]<vmagmin)
    {
      vmagmin = vmag[eye];
    }
  }

char *ffilename[1024];
FILE *out;
  sprintf( ffilename,"./out/vmag_proc%04d.txt", get_proc_id( lattice));
  out = fopen( ffilename, "w");
    for(eye=0; eye<(ni*nj*nk);eye++)
    {
      fprintf(out,"%f\n",vmag[eye]);
    //  vmagscale[eye] = (255.*vmag[eye] / vmagmax);
    //  printf("vmag [%d] = %1.10f,   vmagscale[%d] = %3.5f \n", eye, vmag[eye], eye, vmagscale[eye]);
    }
  fclose(out);
//printf("vmagmax = %f, vmagmin = %f \n", vmagmax, vmagmin);
//printf("velocity [%d] = %1.10f, Squared = %1.10f \n", eye, lattice->macro_vars[subs][eye].u[0], (lattice->macro_vars[subs][eye].u[0]*lattice->macro_vars[subs][eye].u[0]));
#endif /*WRITE_VMAG*/
//-------------------------------------------------------------------------------------------

  } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

#if WRITE_PLOT_FILE
    if(NUM_FLUID_COMPONENTS == 2)
    {
    sprintf( filename, "./out/plt_%dx%dx%d_frame%04d_proc%04d.dat",
        ni, nj, nk, frame, get_proc_id( lattice));
    sprintf( filename2, "./out/R2_%dx%dx%d_frame%04d_proc%04d.dat",
        ni, nj, nk, frame, get_proc_id( lattice));
    printf("Huang tecplot 2 component\n");
    write_plt(
      lattice,
      &(get_rho(lattice,0,0)),
      &(get_rho(lattice,1,0)),
      filename, filename2);

    //get saturation
   basek = 0;
#if PARALLEL
  basek = get_g_SZ( lattice);
#endif
        sprintf( fn3, "./out/saturation_proc%04d.dat",
        get_proc_id( lattice));
  if(frame ==0)
  {
  if(  !( sat = fopen(fn3, "w+")))
   {
    printf(" %04d>> ERROR: "
      "fopen(\"%s\",\"r\") = NULL.  Bye, bye!\n",
      __FILE__, __LINE__,
      get_proc_id( lattice),
      fn3);
    process_exit(1);
   }
  }
  if(frame !=0)
  {
  if(  !( sat = fopen(fn3, "a+")))
   {
    printf(" %04d>> ERROR: "
      "fopen(\"%s\",\"r\") = NULL.  Bye, bye!\n",
      __FILE__, __LINE__,
      get_proc_id( lattice),
      fn3);
    process_exit(1);
   }
  }
  count =0;
  count1 =0;
  count2 =0;

 for( n=0; n<ni*nj*nk; n++)
  {

     if( (N2Z(n,ni,nj,nk)+basek >lattice->param.z1 &&
          N2Z(n,ni,nj,nk)+basek <lattice->param.z2    ) )
     {
       count++;
      if ( get_rho(lattice,0,n)> (lattice->param.rho_A[0]+lattice->param.rho_A[1])/2.) count1 ++;
      if ( lattice->solids[0][n].is_solid) count2++;
     }
  }
  fprintf( sat, " %7d, %7d, %7d, %7d\n",   count, count2, count-count2, count1  );


    sprintf( filename, "./out/plt_%dx%dx%d_frame%04d_uvw.dat",
        ni, nj, nk, frame, get_proc_id( lattice));
    printf("Huang tecplot 2 UVZ component\n");
    write_plt_uvw(
      lattice,
      (get_rho_ptr(lattice,0,0)),
      (get_rho_ptr(lattice,1,0)),
      filename);
 fclose(sat);
    }

    else
    {
    sprintf( filename, "./out/plt_%dx%dx%d_frame%04d_proc%04d.dat",
        ni, nj, nk, frame, get_proc_id( lattice));
    printf("Huang tecplot 1 component\n");
    write_plt_single(
      lattice,
      &(get_rho(lattice,0,0)),
      filename);
    }
#endif

} /* void dump_macro_vars( struct lattice_struct *lattice, int time) */

void read_solids( lattice_ptr lattice, char *filename)
{
  int size, size_read;
  unsigned char *raw;
  int n;
  int g_n;
  int subs;
  FILE *fd;
#if 1
  int i, j, k, p;
#endif

  for( n=0; n<lattice->NumNodes; n++)
  {
    is_solid( lattice, n);
  }

#if 1
  fd = fopen( filename, "r+");
  if( !fd)
  {
    printf("%s %d %04d >> ERROR: Can't open file \"%s\". (Exiting!)\n",
      __FILE__,__LINE__, get_proc_id(lattice), filename);
    process_exit(1);
  }

  size = get_g_NumNodes( lattice)*sizeof(unsigned char);
  if( !( raw = (unsigned char *)malloc(size)))
  {
    printf("%s %d %04d >> read_solids() -- "
        "ERROR: Can't malloc image buffer. (Exiting!)\n", __FILE__, __LINE__, get_proc_id(lattice));
    process_exit(1);
  }

  printf("%s %d %04d >> Reading %d bytes from file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), size, filename);

  //size_read = read( fd, raw, size);
  size_read = fread( raw, 1, size, fd);

  if( size_read != size)
  {
    printf("%s %d %04d >> read_solids() -- "
        "ERROR: Can't read image data: read = %d. (Exiting!)\n",
          __FILE__, __LINE__, get_proc_id(lattice), size_read);
    process_exit(1);
  }

  fclose( fd);

  g_n = get_g_StartNode( lattice); // Global node index.
  for( n=0; n<lattice->NumNodes; n++)
  {
    set_is_solid( lattice, n, raw[g_n]);
    g_n++;
  }

  free(raw);
#else
  int *a;
  a = (int*)malloc( (lattice->NumNodes)*sizeof(int));
  if( !a)
  {
    printf("%s %d %04d >> read_solids() -- "
           "ERROR: Can't malloc array of length %d. (Exiting!)\n",
           __FILE__, __LINE__, get_proc_id(lattice), get_NumNodes(lattice));
    process_exit(1);
  }
  read_raw( lattice, a, lattice->NumNodes, /*stride*/1, filename);
  for( n=0; n<get_NumNodes(lattice); n++)
  {
    set_is_solid( lattice, n, a[n]);
  }
  free(a);
#endif

#if 1
  if( /* Domain not too big. (Tweak to suit.) */
      ( get_LX(lattice) <= 12
      &&
        get_LY(lattice) <= 12
      &&
        get_LZ(lattice) <= 12
      )
    )
  {
#if PARALLEL
  for( p=0; p<get_num_procs( lattice); p++)
  {
    MPI_Barrier( MPI_COMM_WORLD);
    if( p == get_proc_id( lattice))
    {
      printf("%s %d %04d >> Solids:\n", __FILE__, __LINE__, p);
#endif
      for( j=0; j<get_LY( lattice); j++)
      {
        for( k=0; k<get_LZ( lattice); k++)
        {
          for( i=0; i<get_LX( lattice); i++)
          {
            if( is_solid( lattice
                        , XYZ2N( i, j, k
                               , get_LX( lattice)
                               , get_LY( lattice)) ) )
            {
              printf("#");
            }
            else
            {
              printf("-");
            }
          }
          printf(" ");
        }
        printf("\n");
      }
#if PARALLEL
    }
  }
  MPI_Barrier( MPI_COMM_WORLD);
#endif
  }
#endif

} /* read_solids( lattice_ptr lattice, char *filename) */

void read_raw(
       lattice_ptr lattice,
       int    *a,
       int     stride,
       char   *filename )
{
  int size, size_read;
  unsigned char *raw;
  int n;
  int g_n;
  int subs;
  FILE *fd;

  fd = fopen( filename, "r+");
  if( !fd)
  {
    printf("%s %d %04d >> ERROR: Can't open file \"%s\". (Exiting!)\n",
      __FILE__,__LINE__, get_proc_id(lattice), filename);
    process_exit(1);
  }

  size = get_g_NumNodes( lattice)*sizeof(unsigned char);
  if( !( raw = (unsigned char *)malloc(size)))
  {
    printf("%s %d %04d >> read_solids() -- "
           "ERROR: Can't malloc image buffer. (Exiting!)\n",
           __FILE__, __LINE__, get_proc_id(lattice));
    process_exit(1);
  }

  printf("%s %d %04d >> Reading %d bytes from file \"%s\".\n",
         __FILE__, __LINE__, get_proc_id(lattice), size, filename);

  size_read = fread( raw, 1, size, fd);

  if( size_read != size)
  {
    printf("%s %d %04d >> read_solids() -- "
        "ERROR: Can't read image data: read = %d. (Exiting!)\n",
          __FILE__, __LINE__, get_proc_id(lattice), size_read);
    process_exit(1);
  }

  fclose( fd);

  g_n = get_g_StartNode( lattice); // Global node index.
  for( n=0; n<lattice->NumNodes; n++)
  {
    //printf("%s %d >> raw[%d] = %d\n",__FILE__,__LINE__,g_n,raw[g_n]);
    for( subs=0; subs<(NUM_FLUID_COMPONENTS); subs++)
    {
      a[n] = raw[g_n];
    }
    //printf("%s %d %04d >> solids[%d] = %d.\n",
    //    __FILE__, __LINE__, get_proc_id(lattice), n, (int)lattice->solids[0][n].is_solid);
    g_n++;
  }

  free(raw);

} /* read_solids( lattice_ptr lattice, char *filename) */

#if POROUS_MEDIA
void read_ns( lattice_ptr lattice, char *filename)
{
  int size, size_read;
  unsigned char *raw;
  int n;
  int g_n;
  int subs;
  FILE *fd;
#if 1
  int i, j, k, p;
#endif

  for( n=0; n<lattice->NumNodes; n++)
  {
    lattice->ns[n].ns = 0;
  }

  fd = fopen( filename, "r+");
  if( !fd)
  {
    printf("%s %d %04d >> ERROR: Can't open file \"%s\". (Exiting!)\n",
      __FILE__,__LINE__, get_proc_id(lattice), filename);
    process_exit(1);
  }

  size = get_g_NumNodes( lattice)*sizeof(unsigned char);
  if( !( raw = (unsigned char *)malloc(size)))
  {
    printf("%s %d %04d >> read_solids() -- "
        "ERROR: Can't malloc image buffer. (Exiting!)\n", __FILE__, __LINE__, get_proc_id(lattice));
    process_exit(1);
  }

  printf("%s %d %04d >> Reading %d bytes from file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), size, filename);

  //size_read = read( fd, raw, size);
  size_read = fread( raw, 1, size, fd);

  if( size_read != size)
  {
    printf("%s %d %04d >> read_solids() -- "
        "ERROR: Can't read image data: read = %d. (Exiting!)\n",
          __FILE__, __LINE__, get_proc_id(lattice), size_read);
    process_exit(1);
  }

  fclose( fd);

  g_n = get_g_StartNode( lattice); // Global node index.
  for( n=0; n<lattice->NumNodes; n++)
  {
    //printf("%s %d >> raw[%d] = %d\n",__FILE__,__LINE__,g_n,raw[g_n]);
    lattice->ns[n].ns = raw[g_n]/255.0;
    //printf("%s %d %04d >> ns[%d] = %f.\n",
    //    __FILE__, __LINE__, get_proc_id(lattice), n,
    //    lattice->ns[n].ns);
    g_n++;
  }

  free(raw);

#if 0
  if( /* Domain not too big. (Tweak to suit.) */
      ( get_LX(lattice) <= 12
      &&
        get_LY(lattice) <= 12
      &&
        get_LZ(lattice) <= 12
      )
    )
  {
#if PARALLEL
  for( p=0; p<get_num_procs( lattice); p++)
  {
    MPI_Barrier( MPI_COMM_WORLD);
    if( p == get_proc_id( lattice))
    {
      printf("%s %d %04d >> Solids:\n", __FILE__, __LINE__, p);
#endif
      int ni = get_LX(lattice);
      int nj = get_LY(lattice);
      for( j=0; j<get_LY( lattice); j++)
      {
        for( k=0; k<get_LZ( lattice); k++)
        {
          for( i=0; i<get_LX( lattice); i++)
          {
            if(    lattice->ns[XYZ2N(i,j,k,ni,nj)].ns > 32
                &&
                   lattice->ns[XYZ2N(i,j,k,ni,nj)].ns < 127)
            {
              printf("%c", (char)(lattice->ns[XYZ2N(i,j,k,ni,nj)].ns));
            }
            else
            {
              if(    lattice->ns[XYZ2N(i,j,k,ni,nj)].ns-128 > 32
                  &&
                     lattice->ns[XYZ2N(i,j,k,ni,nj)].ns-128 != 127)
              {
                printf("%c", (char)(lattice->ns[XYZ2N(i,j,k,ni,nj)].ns-128));
              }
              else
              {
                printf(" ");
              }
            }
          }
          printf(" ");
        }
        printf("\n");
      }
#if PARALLEL
    }
  }
  MPI_Barrier( MPI_COMM_WORLD);
#endif
  }
#endif

} /* read_solids( lattice_ptr lattice, char *filename) */
#endif

/*void read_solids_from_plt( lattice_ptr lattice, char *filename)
{
  int size, size_read;
  unsigned char *raw;
  int n;
  int g_n;
  int subs;
  FILE *in;
#if 1
  int i, j, k, p;
#endif


  for( n=0; n<lattice->NumNodes; n++)
  {
    lattice->solids[0][n].is_solid = 0;
  }


   if( !( in = fopen( filename, "r")))
  {
    printf("%s %d %04d>> ERROR: "
      "fopen(\"%s\",\"r\") = NULL.  Bye, bye!\n",
      __FILE__, __LINE__,
      get_proc_id( lattice),
      infile);
    process_exit(1);
  }

  size = get_g_NumNodes( lattice)*sizeof(unsigned char);
  if( !( raw = (unsigned char *)malloc(size)))
  {
    printf("%s %d %04d >> read_solids() -- "
        "ERROR: Can't malloc image buffer. (Exiting!)\n", __FILE__, __LINE__, get_proc_id(lattice));
    process_exit(1);
  }
 for (i =1 ; i<=7; i++)
 {
   skip_label( in);
 }

  for(k =1; k<=lattice->param.LZ; k++)
  {
  for(j =1; j<=lattice->param.LY; j++)
  {
  for(i =1; i<=lattice->param.LX; i++)
  {
   fscanf( in, "%d", &( raw[n])         );
  }
  }
  }

  fclose( in);

  for( n=0; n<lattice->NumNodes; n++)
  {
    for( subs=0; subs<(NUM_FLUID_COMPONENTS); subs++)
    {
      lattice->solids[subs][n].is_solid = raw[g_n]*255;
    }
    //printf("%s %d %04d >> solids[%d] = %d.\n",
    //    __FILE__, __LINE__, get_proc_id(lattice), n, (int)lattice->solids[0][n].is_solid);
    g_n++;
  }

  free(raw);

#if 1
  if( // Domain not too big. (Tweak to suit.)
      ( get_LX(lattice) <= 12
      &&
        get_LY(lattice) <= 12
      &&
        get_LZ(lattice) <= 12
      )
    )
  {
#if PARALLEL
  for( p=0; p<get_num_procs( lattice); p++)
  {
    MPI_Barrier( MPI_COMM_WORLD);
    if( p == get_proc_id( lattice))
    {
      printf("%s %d %04d >> Solids:\n", __FILE__, __LINE__, p);
#endif
      for( j=0; j<get_LY( lattice); j++)
      {
        for( k=0; k<get_LZ( lattice); k++)
        {
          for( i=0; i<get_LX( lattice); i++)
          {
            if( is_solid( lattice,
                          XYZ2N( i, j, k,
                                 get_LX( lattice),
                                 get_LY( lattice)) ) )
            {
              printf("#");
            }
            else
            {
              printf(" ");
            }
          }
          printf(" ");
        }
        printf("\n");
      }
#if PARALLEL
    }
  }
  MPI_Barrier( MPI_COMM_WORLD);
#endif
  }
#endif
}
*/

void write_raw(
       lattice_ptr lattice,
       real *a,
       int     stride,
       real  a_max,
       real  a_min,
       char   *filename)
{
  int size;
  int n, ni, nj, nk;
  unsigned char *Xpix1;
  FILE *infile;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  size = ni*nj*nk*sizeof(unsigned char);

  if(!(Xpix1 = ( unsigned char*)malloc( size)))
  {
    printf("%s %d %04d >> ERROR: Can't malloc Xpix1. (Exiting!)\n",
      __FILE__,__LINE__, get_proc_id(lattice));
    process_exit(1);
  }
  for( n=0; n<ni*nj*nk; n++)
  {
    if( !is_solid( lattice, n))
    {
//      Xpix1[n] = (unsigned char)ROUND(255.*(a[stride*n]-a_min)/(a_max-a_min));
      Xpix1[n] = (unsigned char)ROUND(1.+254.*(a[stride*n]-a_min)/(a_max-a_min));
    }
    else
    {
      Xpix1[n] = (unsigned char)0;
    }
    //printf("%s %d %04d >> write_raw() -- Xpix1[%d] = 255*ROUND(%f) = %d.\n",
    //    __FILE__, __LINE__, get_proc_id(lattice),
    //  n, (a[stride*n]-a_min)/(a_max-a_min), (int)Xpix1[n]);
  }

  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fwrite( (char*)Xpix1, 1, size, infile);

  fclose(infile);

  printf("%s %d %04d >> write_raw() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);
} /* void write_raw( lattice_ptr lattice, real *a,  ... */

void write_raw_u(
       lattice_ptr lattice,
       real *ux,
       real *uy,
       real *uz,
       int     stride,
       real  ux_max,
       real  ux_min,
       real  uy_max,
       real  uy_min,
       real  uz_max,
       real  uz_min,
       char   *filename)
{
  int size;
  int n, ni, nj, nk;
  unsigned char *Xpix1;
  FILE *infile;
  real a_max = sqrt( ux_max*ux_max + uy_max*uy_max + uz_max*uz_max);
  real a_min = sqrt( ux_min*ux_min + uy_min*uy_min + uz_min*uz_min);
  real a;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  size = ni*nj*nk*sizeof(unsigned char);

  if(!(Xpix1 = ( unsigned char*)malloc( size)))
  {
    printf("%s %d %04d >> ERROR: Can't malloc Xpix1. (Exiting!)\n",
      __FILE__,__LINE__, get_proc_id(lattice));
    process_exit(1);
  }
  for( n=0; n<ni*nj*nk; n++)
  {
// printf("n=%d\n",n);
    if( !is_solid( lattice, n))
    {
      // Xpix1[n] =
      //   (unsigned char)ROUND(255.*(a[stride*n]-a_min)/(a_max-a_min));

// printf("&ux = %x\n",ux);
// printf("  ux = %f\n",ux[stride*n]);
// printf("  uy = %f\n",uy[stride*n]);
// printf("  uz = %f\n",
//   ((get_NumDims(lattice)==2)?(0):(uz[stride*n]*uz[stride*n])));

      a = sqrt( ux[stride*n]*ux[stride*n]
              + uy[stride*n]*uy[stride*n]
              + ((get_NumDims(lattice)==2)?(0):(uz[stride*n]*uz[stride*n])));
      Xpix1[n] = (unsigned char)ROUND(1.+254.*(a-a_min)/(a_max-a_min));
    }
    else
    {
      Xpix1[n] = (unsigned char)0;
    }
    //printf("%s %d %04d >> write_raw() -- Xpix1[%d] = 255*ROUND(%f) = %d.\n",
    //    __FILE__, __LINE__, get_proc_id(lattice),
    //  n, (a[stride*n]-a_min)/(a_max-a_min), (int)Xpix1[n]);
  }

  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fwrite( (char*)Xpix1, 1, size, infile);

  fclose(infile);

  printf("%s %d %04d >> write_raw() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);
} /* void write_raw( lattice_ptr lattice, real *a,  ... */


#ifdef __CUDACC__
// The definitions of rho2bmp and u2bmp need to be relocated to lattice.c for
// nvcc so that they will be defined when write_rho_image and write_u_image
// are defined. This is due to the ad hoc way we are compiling the multiple
// source code files and the way that nvcc works.
#else
// R H O 2 B M P  {{{
//##############################################################################
// void rho2bmp( char *filename, int time)
//
void rho2bmp( lattice_ptr lattice, int time)
{
  FILE   *in,
         *o;
  int    i, j,
         n, m;
  int    pad,
         bytes_per_row;
  int    frame;
  char   k;
  char   b;
  struct bitmap_file_header bmfh;
  struct bitmap_info_header bmih;
  struct rgb_quad rgb;
  int    *int_ptr;
  short  int *short_int_ptr;
  int    *width_ptr;
  int    *height_ptr;
  short  int *bitcount_ptr;
  char   filename[1024];
  char   red_val,
         green_val,
         blue_val,
         val;
  real fval;
  real min_rho, max_rho;
  int    subs;
  int    num_colors;

#if SAY_HI
  printf("rho2bmp() -- Hi!\n");
#endif /* SAY_HI */

 for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
 {
   frame = time/lattice->param.FrameRate;

   bmfh.bfType[0] = 'B';
   bmfh.bfType[1] = 'M';
   *((int*)bmfh.bfSize)= get_LY(lattice)*(
       (int)ceil( ( ((real)get_LX(lattice))*( /*depth*/24.))/8.) // bytes per row
       + ( 4 - (int)ceil( ( ((real)get_LX(lattice))*( /*depth*/24.))/8.) % 4) % 4 // pad
       );
   *((short int*)bmfh.bfReserved1) = 0;
   *((short int*)bmfh.bfReserved2) = 0;
   *((int*)bmfh.bfOffBits) = 54; // 14 byte file header and 40 byte info header

   *((int*)bmih.biSize) = 40;
   *((int*)bmih.biWidth) = get_LX(lattice);
   *((int*)bmih.biHeight) = get_LY(lattice);
   *((short int*)bmih.biPlanes) = 1;
   *((short int*)bmih.biBitCount) = 24;
   *((int*)bmih.biCompression) = 0;
   *((int*)bmih.biSizeImage) = 0;
   *((int*)bmih.biXPelsPerMeter) = 0;
   *((int*)bmih.biYPelsPerMeter) = 0;
   *((int*)bmih.biClrUsed) = 0;
   *((int*)bmih.biClrImportant) = 0;

   width_ptr = (int*)bmih.biWidth;
   height_ptr = (int*)bmih.biHeight;
   bitcount_ptr = (short int*)bmih.biBitCount;

   // Bytes per row of the bitmap.
   bytes_per_row =
     ((int)ceil(( ( ((real)(ENDIAN4(*width_ptr)))
                  * ((real)(ENDIAN2(*bitcount_ptr))) )/8.)));

   // Bitmaps pad rows to preserve 4-byte boundaries.
   // The length of a row in the file will be bytes_per_row + pad .
   pad = ((4) - bytes_per_row%4)%4;

   max_rho = get_max_rho( lattice, subs);
   min_rho = get_min_rho( lattice, subs);

   gen_filename( lattice, filename, "rho", get_frame(lattice), subs, ".bmp");
   if( !( o = fopen( filename, "w+")))
   {
     printf("ERROR: fopen( \"%s\", \"w+\") = NULL.  Bye, bye!\n", filename);
     process_exit(1);
   }

   fwrite( &bmfh, sizeof(struct bitmap_file_header), 1, o );
   fwrite( &bmih, sizeof(struct bitmap_info_header), 1, o );

   for( j=0; j<get_LY(lattice); j++)
   {
     n = j*get_LX(lattice);

     for( i=0; i<get_LX(lattice); i++, n++)
     {
       if( is_not_solid(lattice,n))
       {
         if( subs==0)
         {
           if( lattice->param.plot_scale_dynamic)
           {
             if( max_rho!=min_rho)
             {
               fval = ROUND( 255.*( get_rho(lattice,subs,n)
                     - min_rho)
                   /( max_rho-min_rho));
             }
             else
             {
               fval = 255.;
             }
           }
           else
           {
             fval = ROUND( 255.*(get_rho(lattice,subs,n)
                   /( (lattice->param.rho_A[subs]>lattice->param.rho_B[subs])
                     ?(lattice->param.rho_A[subs])
                     :(lattice->param.rho_B[subs]) )
                   ));
           }
           if( fval >= 0.)
           {
             if( fval <= 255.)
             {
               red_val   = (char)((int)(255. - fval)%256);
               green_val = (char)((int)(255. - fval)%256);
               blue_val  = (char)255;
             }
             else
             {
               red_val   = (char)0;
               green_val = (char)0;
               blue_val  = (char)255;
             }
           }
           else
           {
             red_val   = (char)((int)(255. + fval)%256);
             green_val = (char)((int)(255. + fval)%256);
             blue_val  = (char)((int)(255. + fval)%256);
             // TODO: Issue warning or something? Potential instability?
           }
         } /* if( subs==0) */

         else // subs == 1
         {
           if( lattice->param.plot_scale_dynamic)
           {
             if( max_rho!=min_rho)
             {
               fval = ROUND( 255.*( get_rho(lattice,subs,n)
                     - min_rho)
                   /( max_rho-min_rho));
             }
             else
             {
               fval = 0.;
             }
           }
           else
           {
             //printf("%s (%d) >> fval = %f -> ", __FILE__, __LINE__, fval);
#if INAMURO_SIGMA_COMPONENT
             fval = ROUND( 255.*(get_rho(lattice,subs,n))
                 /(lattice->param.rho_sigma));
#else /* !( INAMURO_SIGMA_COMPONENT) */
             fval = ROUND( 255.*(get_rho(lattice,subs,n)
                   /( (lattice->param.rho_A[subs]>lattice->param.rho_B[subs])
                     ?(lattice->param.rho_A[subs])
                     :(lattice->param.rho_B[subs]) )
                   ));
#endif /* INAMURO_SIGMA_COMPONENT */
             //printf("%f\n", fval);
           }
           if( fval >= 0.)
           {
             if( fval <= 255.)
             {
               red_val   = (char)255;
               green_val = (char)((int)(255. - fval)%256);
               blue_val  = (char)((int)(255. - fval)%256);
             }
             else
             {
               red_val   = (char)255;//((int)(255. - (fval - 255.))%256);
               green_val = (char)  0;//((int)(255. - (fval - 255.))%256);
               blue_val  = (char)  0;//((int)(255. - (fval - 255.))%256);
             }
           }
           else
           {
             red_val   = (char)((int)(255. + fval)%256);
             green_val = (char)((int)(255. + fval)%256);
             blue_val  = (char)((int)(255. + fval)%256);
             // TODO: Issue a warning or something?  Potential instability?
           }

         } /* if( subs==0) else */

       } /* if( lattice->bc[subs][ n].bc_type == 0) */

       else // lattice->bc[subs][ n].bc_type != 0
       {
#if SOLID_COLOR_IS_CHECKERBOARD
         // Checkerboard pattern over the solids and boundary conditions.
         if( (i+j)%2)
         {
           red_val   = (char)200;
           green_val = (char)200;
           blue_val  = (char)200;
           val       = (char)200;
         }
         else
         {
           red_val   = (char)184;
           green_val = (char)184;
           blue_val  = (char)184;
           val       = (char)184;
         }
#else /* !( SOLID_COLOR_IS_CHECKERBOARD) */
#if SOLID_COLOR_IS_BLACK
         red_val   = (char)0;
         green_val = (char)0;
         blue_val  = (char)0;
         val       = (char)0;
#else /* !( SOLID_COLOR_IS_BLACK) */
         red_val   = (char)255;
         green_val = (char)255;
         blue_val  = (char)255;
         val       = (char)255;
#endif /* SOLID_COLOR_IS_BLACK */
#endif /* SOLID_COLOR_IS_CHECKERBOARD */

       } /* if( lattice->bc[subs][ n].bc_type == 0) else */

       //printf("blue_val( %d, %d) = %d\n", i, j, (int)blue_val);

       if( fwrite( &blue_val, 1, 1, o) != 1) { printf("BOOM!\n"); process_exit(1);}
       //printf("BING %d %d\n", i, j);
       if( fwrite( &green_val, 1, 1, o) != 1) { printf("BOOM!\n"); process_exit(1);}
       //printf("BING %d %d\n", i, j);
       if( fwrite( &red_val, 1, 1, o) != 1) { printf("BOOM!\n"); process_exit(1);}
       //printf("BING %d %d\n", i, j);

     } /* for( i=0; i<get_LX(lattice); i++) */

     // Pad for 4-byte boundaries.
     val = (char)0;
     for( i=0; i<pad; i++)
     {
       if( fwrite( &val, 1, 1, o) != 1) { printf("BOOM!\n"); process_exit(1);}
     }

   } /* for( j=0; j<get_LY(lattice); j++) */

   fclose(o);

#if VERBOSITY_LEVEL > 0
   printf("rho2bmp() -- Wrote file \"%s\".\n", filename);
#endif /* VERBOSITY_LEVEL > 0 */

 } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

#if SAY_HI
  printf("rho2bmp() -- Bye!\n");
  printf("\n");
#endif /* SAY_HI */

} /* rho2bmp( lattice_ptr lattice, int time) */
// U 2 B M P  {{{
//##############################################################################
// void u2bmp( char *filename, int time)
//
void u2bmp( lattice_ptr lattice, int time)
{
  FILE   *in,
         *o_u,
         *o_ux,
         *o_uy;
  int    i, j,
         n, m;
  int    pad,
         bytes_per_row;
  int    frame;
  char   k;
  char   b;
  struct bitmap_file_header bmfh;
  struct bitmap_info_header bmih;
  struct rgb_quad rgb;
  int    *int_ptr;
  short  int *short_int_ptr;
  int    *width_ptr;
  int    *height_ptr;
  short  int *bitcount_ptr;
  char   filename[1024];
  char   red_val,
         green_val,
         blue_val,
         val;
  real max_u[2], maxu;
  real u_x, u_y, u;
  int    subs;

#if SAY_HI
  printf("u2bmp() -- Hi!\n");
#endif /* SAY_HI */

 for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
 {

  frame = time/lattice->param.FrameRate;

#if 0
  sprintf( filename, "./in/%dx%d_proc%04d.bmp",
      get_LX(lattice), get_LY(lattice), get_proc_id(lattice));
  if( !( in = fopen( filename, "r")))
  {
    printf("%s %d >> u2bmp() -- Error opening file \"%s\".\n",
      __FILE__,__LINE__,filename);
    process_exit(1);
  }

  // n = fread( void *BUF, size_t SIZE, size_t COUNT, FILE *FP);

  n = fread( &bmfh, sizeof(struct bitmap_file_header), 1, in );
  if( strncmp(bmfh.bfType,"BM",2))
  {
    printf("ERROR: Can't process this file type.  Exiting!\n");
    printf("\n");
    process_exit(1);
  }
  n = fread( &bmih, sizeof(struct bitmap_info_header), 1, in );
  int_ptr = (int*)bmih.biCompression;
  if( *int_ptr != 0)
  {
    printf("ERROR: Can't handle compression.  Exiting!\n");
    printf("\n");
    process_exit(1);
  }

#if 0
  *((int*)(bmih.biWidth)) = ENDIAN4(((int)(*((int*)(bmih.biWidth)))));
  *((int*)(bmih.biHeight)) = ENDIAN4(((int)(*((int*)(bmih.biHeight)))));
  *((short int*)(bmih.biBitCount)) = ENDIAN2(((short int)(*((short int*)(bmih.biBitCount)))));
#endif

  width_ptr    = (int*)bmih.biWidth;
  height_ptr   = (int*)bmih.biHeight;
  bitcount_ptr = (short int*)bmih.biBitCount;

printf("%s %d >> width    = %d\n",__FILE__,__LINE__, ENDIAN4(*width_ptr)   );
printf("%s %d >> height   = %d\n",__FILE__,__LINE__, ENDIAN4(*height_ptr)  );
printf("%s %d >> bitcount = %d\n",__FILE__,__LINE__, ENDIAN2(*bitcount_ptr));

  // Read palette entries, if applicable.
  if( ENDIAN2(*bitcount_ptr) < 24)
  {
    n = (int)pow(2.,(real)ENDIAN2(*bitcount_ptr)); // Num palette entries.
    for( i=0; i<n; i++)
    {
      k = fread( &rgb, sizeof(struct rgb_quad), 1, in );
      if( k!=1)
      {
        printf("Error reading palette entry %d.  Exiting!\n", i);
        process_exit(1);
      }
    }
  }

  fclose(in);

#else
  bmfh.bfType[0] = 'B';
  bmfh.bfType[1] = 'M';
  *((int*)bmfh.bfSize)= get_LY(lattice)*(
    (int)ceil( ( ((real)get_LX(lattice))*( /*depth*/24.))/8.) // bytes per row
  + ( 4 - (int)ceil( ( ((real)get_LX(lattice))*( /*depth*/24.))/8.) % 4) % 4 // pad
  );
  *((short int*)bmfh.bfReserved1) = 0;
  *((short int*)bmfh.bfReserved2) = 0;
  *((int*)bmfh.bfOffBits) = 54; // 14 byte file header and 40 byte info header

  *((int*)bmih.biSize) = 40;
  *((int*)bmih.biWidth) = get_LX(lattice);
  *((int*)bmih.biHeight) = get_LY(lattice);
  *((short int*)bmih.biPlanes) = 1;
  *((short int*)bmih.biBitCount) = 24;
  *((int*)bmih.biCompression) = 0;
  *((int*)bmih.biSizeImage) = 0;
  *((int*)bmih.biXPelsPerMeter) = 0;
  *((int*)bmih.biYPelsPerMeter) = 0;
  *((int*)bmih.biClrUsed) = 0;
  *((int*)bmih.biClrImportant) = 0;

  width_ptr = (int*)bmih.biWidth;
  height_ptr = (int*)bmih.biHeight;
  bitcount_ptr = (short int*)bmih.biBitCount;

#endif

  // Bytes per row of the bitmap.
  bytes_per_row =
    ((int)ceil(( (((real)(ENDIAN4(*width_ptr)))*((real)(ENDIAN2(*bitcount_ptr))))/8.)));

  // Bitmaps pad rows to preserve 4-byte boundaries.
  // The length of a row in the file will be bytes_per_row + pad .
  pad = ((4) - bytes_per_row%4)%4;

  max_u[0] = get_max_ux(lattice,subs);
  max_u[1] = get_max_uy(lattice,subs);
  max_u[2] = get_max_uz(lattice,subs);

  gen_filename( lattice, filename, "u", get_frame(lattice), subs, ".bmp");
  if( !( o_u = fopen( filename, "w+")))
  {
    printf("ERROR: fopen( \"%s\", \"w+\") = NULL.  Bye, bye!\n", filename);
    process_exit(1);
  }

  gen_filename( lattice, filename, "u_x", get_frame(lattice), subs, ".bmp");
  if( !( o_ux = fopen( filename, "w+")))
  {
    printf("ERROR: fopen( \"%s\", \"w+\") = NULL.  Bye, bye!\n", filename);
    process_exit(1);
  }

  gen_filename( lattice, filename, "u_y", get_frame(lattice), subs, ".bmp");
  if( !( o_uy = fopen( filename, "w+")))
  {
    printf("ERROR: fopen( \"%s\", \"w+\") = NULL.  Bye, bye!\n", filename);
    process_exit(1);
  }

  fwrite( &bmfh, sizeof(struct bitmap_file_header), 1, o_u );
  fwrite( &bmih, sizeof(struct bitmap_info_header), 1, o_u );

  fwrite( &bmfh, sizeof(struct bitmap_file_header), 1, o_ux );
  fwrite( &bmih, sizeof(struct bitmap_info_header), 1, o_ux );

  fwrite( &bmfh, sizeof(struct bitmap_file_header), 1, o_uy);
  fwrite( &bmih, sizeof(struct bitmap_info_header), 1, o_uy);

  //for( j=get_LY(lattice)-1; j>=0; j--)
  for( j=0; j<get_LY(lattice); j++)
  {
    n = j*get_LX(lattice);

    for( i=0; i<get_LX(lattice); i++, n++)
    {
      if( is_not_solid(lattice,n))
      {
#if 1
        blue_val  = (char)0;
        green_val = (char)0;
        red_val   = (char)0;

        u_x = (get_ux(lattice,subs,n));
        u_y = (get_uy(lattice,subs,n));

        u = sqrt(u_x*u_x + u_y*u_y);
        maxu = sqrt( max_u[0]*max_u[0] + max_u[1]*max_u[1]);

        if( is_solid(lattice,n))
        {
          blue_val  = (char)ROUND( 128.*fabs(u_x)/max_u[0]);
          green_val = (char)ROUND( 128.*fabs(u_y)/max_u[1]);
        }
        else
        {
#if 0
          blue_val  = (char)ROUND(  255.*fabs(u_x)/max_u[0]);
          green_val = (char)ROUND(  255.*fabs(u_y)/max_u[1]);
          red_val   = 0.;//(char)ROUND( 128.*fabs(u)/maxu);
#else
          blue_val  = (char)ROUND(  255.*((fabs(u_x)!=0.)
                    ? (fabs(u_x)/max_u[0])
                    : (0.)));
          green_val = (char)ROUND(  255.*((fabs(u_y)!=0.)
                    ? (fabs(u_y)/max_u[1])
                    : (0.)));
          red_val   = 0.;
          //red_val   = (char)ROUND( 128.*((fabs(u  )!=0.)
          //          ? (fabs(u  )/maxu)
          //          : (0.)));
#endif

        }
#else
        blue_val  = (char)255;
        green_val = (char)255;
        red_val   = (char)255;

        u = sqrt(u_x*u_x + u_y*u_y);
        maxu = sqrt( max_u[0]*max_u[0] + max_u[1]*max_u[1]);

        //if( fabs(u) > .1*maxu)
        //{
          green_val  = (char)ROUND( 255.-255.*fabs(u)/maxu);
          red_val    = (char)ROUND( 255.-255.*fabs(u)/maxu);
          blue_val    = (char)ROUND( 255.-255.*fabs(u)/maxu);
        //}
        //else
        //{
        //  green_val  = (char)0;
        //  red_val    = (char)0;
        //  blue_val    = (char)0;
        //}
#endif

        val = (char)0;

      } /* if( is_not_solid(lattice,n)) */

      else // is_solid(lattice,n)
      {
#if SOLID_COLOR_IS_CHECKERBOARD
        // Checkerboard pattern over the solids and boundary conditions.
        if( (i+j)%2)
        {
          red_val   = (char)200;
          green_val = (char)200;
          blue_val  = (char)200;
          val       = (char)200;
        }
        else
        {
          red_val   = (char)184;
          green_val = (char)184;
          blue_val  = (char)184;
          val       = (char)184;
        }
#else /* !( SOLID_COLOR_IS_CHECKERBOARD) */
#if SOLID_COLOR_IS_BLACK
        red_val   = (char)0;
        green_val = (char)0;
        blue_val  = (char)0;
        val       = (char)0;
#else /* !( SOLID_COLOR_IS_BLACK) */
        red_val   = (char)255;
        green_val = (char)255;
        blue_val  = (char)255;
        val       = (char)255;
#endif /* SOLID_COLOR_IS_BLACK */
#endif /* SOLID_COLOR_IS_CHECKERBOARD */

      } /* if( is_not_solid(lattice,n)) else */

#if MARK_ORIGIN_FOR_REFERENCE
 // Mark the origin for reference.
 if( ( i == 0 && j == 0))
 {
   red_val   = (char)255;
   green_val = (char)255;
   blue_val  = (char)255;
   val       = (char)255;
 }
#endif /* MARK_ORIGIN_FOR_REFERENCE */

      if( fwrite( &blue_val,  1, 1, o_u ) != 1)
      { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &green_val, 1, 1, o_u ) != 1)
      { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &red_val,   1, 1, o_u ) != 1)
      { printf("BOOM!\n"); process_exit(1);}

      if( is_not_solid(lattice,n))
      {
        blue_val  = (char)0;
        green_val = (char)0;
        red_val   = (char)0;

        u_x = (get_ux(lattice,subs,n));
        u_y = (get_uy(lattice,subs,n));

        val = (char)0;

        if( !is_solid( lattice, n))
        {
          if( u_x > 0)
          {
            red_val = val;
            blue_val = (char)ROUND( 128.*fabs(u_x)/max_u[0]);
          }
          else
          {
            red_val = (char)ROUND( 128.*fabs(u_x)/max_u[0]);
            blue_val = val;
          }

        }

        else
        {
          if( u_x > 0)
          {
            red_val = val;
            blue_val = (char)ROUND( 255.*((fabs(u_x)!=0.)?(fabs(u_x)/max_u[0]):(0.)));
          }
          else
          {
            red_val = (char)ROUND( 255.*((fabs(u_x)!=0.)?(fabs(u_x)/max_u[0]):(0.)));
            blue_val = val;
          }

        }

      } /* if( lattice->bc[subs][ n].bc_type == 0) */

      else // lattice->bc[subs][ n].bc_type != 0
      {
#if SOLID_COLOR_IS_CHECKERBOARD
        // Checkerboard pattern over the solids and boundary conditions.
        if( (i+j)%2)
        {
          red_val   = (char)200;
          green_val = (char)200;
          blue_val  = (char)200;
          val       = (char)200;
        }
        else
        {
          red_val   = (char)184;
          green_val = (char)184;
          blue_val  = (char)184;
          val       = (char)184;
        }
#else /* !( SOLID_COLOR_IS_CHECKERBOARD) */
#if SOLID_COLOR_IS_BLACK
        red_val   = (char)0;
        green_val = (char)0;
        blue_val  = (char)0;
        val       = (char)0;
#else /* !( SOLID_COLOR_IS_BLACK) */
        red_val   = (char)255;
        green_val = (char)255;
        blue_val  = (char)255;
        val       = (char)255;
#endif /* SOLID_COLOR_IS_BLACK */
#endif /* SOLID_COLOR_IS_CHECKERBOARD */

      } /* if( lattice->bc[subs][ n].bc_type == 0) else */

#if MARK_ORIGIN_FOR_REFERENCE
 // Mark the origin for reference.
 if( ( i == 0 && j == 0))
 {
   red_val   = (char)255;
   green_val = (char)255;
   blue_val  = (char)255;
   val       = (char)255;
 }
#endif /* MARK_ORIGIN_FOR_REFERENCE */

      if( fwrite( &blue_val,  1, 1, o_ux) != 1)
      { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &val,       1, 1, o_ux) != 1)
      { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &red_val,   1, 1, o_ux) != 1)
      { printf("BOOM!\n"); process_exit(1);}

      if( is_not_solid(lattice,n))
      {
        blue_val  = (char)0;
        green_val = (char)0;
        red_val   = (char)0;

        u_x = (get_ux(lattice,subs,n));
        u_y = (get_uy(lattice,subs,n));

        val = (char)0;

        if( !is_solid( lattice, n))
        {
          if( u_y > 0)
          {
            red_val = val;
            green_val = (char)ROUND( 128.*((fabs(u_y)!=0.)?(fabs(u_y)/max_u[1]):(0.)));
          }
          else
          {
            red_val = (char)ROUND( 128.*((fabs(u_y)!=0.)?(fabs(u_y)/max_u[1]):(0.)));
            green_val = val;
          }

        }

        else
        {
          if( u_y > 0)
          {
            blue_val = (char)ROUND( 255.*((fabs(u_y)!=0.)?(fabs(u_y)/max_u[1]):(0.)));
            green_val = val;
          }
          else
          {
            blue_val = val;
            green_val = (char)ROUND( 255.*((fabs(u_y)!=0.)?(fabs(u_y)/max_u[1]):(0.)));
          }

        }

      } /* if( lattice->bc[subs][ n].bc_type == 0) */

      else // lattice->bc[subs][ n].bc_type != 0
      {
#if SOLID_COLOR_IS_CHECKERBOARD
        // Checkerboard pattern over the solids and boundary conditions.
        if( (i+j)%2)
        {
          red_val   = (char)200;
          green_val = (char)200;
          blue_val  = (char)200;
          val       = (char)200;
        }
        else
        {
          red_val   = (char)184;
          green_val = (char)184;
          blue_val  = (char)184;
          val       = (char)184;
        }
#else /* !( SOLID_COLOR_IS_CHECKERBOARD) */
#if SOLID_COLOR_IS_BLACK
        red_val   = (char)0;
        green_val = (char)0;
        blue_val  = (char)0;
        val       = (char)0;
#else /* !( SOLID_COLOR_IS_BLACK) */
        red_val   = (char)255;
        green_val = (char)255;
        blue_val  = (char)255;
        val       = (char)255;
#endif /* SOLID_COLOR_IS_BLACK */
#endif /* SOLID_COLOR_IS_CHECKERBOARD */

      } /* if( lattice->bc[subs][ n].bc_type == 0) else */

#if MARK_ORIGIN_FOR_REFERENCE
 // Mark the origin for reference.
 if( ( i == 0 && j == 0))
 {
   red_val   = (char)255;
   green_val = (char)255;
   blue_val  = (char)255;
   val       = (char)255;
 }
#endif /* MARK_ORIGIN_FOR_REFERENCE */

      if( fwrite( &blue_val,  1, 1, o_uy) != 1)
      { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &green_val, 1, 1, o_uy) != 1)
      { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &red_val,   1, 1, o_uy) != 1)
      { printf("BOOM!\n"); process_exit(1);}

    } /* for( i=0; i<get_LY(lattice); i++) */

    // Pad for 4-byte boundaries.
    val = (char)0;
    for( i=0; i<pad; i++)
    {
      if( fwrite( &val, 1, 1, o_u ) != 1) { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &val, 1, 1, o_ux) != 1) { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &val, 1, 1, o_uy) != 1) { printf("BOOM!\n"); process_exit(1);}
    }

  } /* for( j=0; j<get_LY(lattice); j++) */

  fclose(o_u );
  fclose(o_ux);
  fclose(o_uy);

#if VERBOSITY_LEVEL > 0
  sprintf( filename, "%s/u%dx%d_frame%04d_subs%02d_proc%04d.bmp", get_out_path(lattice),
      get_LX(lattice),
      get_LY(lattice),
      frame, subs, get_proc_id(lattice));
  printf("u2bmp()   -- Wrote file \"%s\".\n", filename);

  sprintf( filename, "%s/u_x%dx%d_frame%04d_subs%02d_proc%04d.bmp", get_out_path(lattice),
      get_LX(lattice),
      get_LY(lattice),
      frame, subs, get_proc_id(lattice));
  printf("u2bmp()   -- Wrote file \"%s\".\n", filename);

  sprintf( filename, "%s/u_y%dx%d_frame%04d_subs%02d_proc%04d.bmp", get_out_path(lattice),
      get_LX(lattice),
      get_LY(lattice),
      frame, subs, get_proc_id(lattice));
  printf("u2bmp()   -- Wrote file \"%s\".\n", filename);
#endif /* VERBOSITY_LEVEL > 0 */

 } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

#if STORE_U_COMPOSITE

  frame = time/lattice->param.FrameRate;

#if !(PARALLEL)
  sprintf( filename, "./in/%dx%d.bmp",
      get_LX(lattice), get_LY(lattice));
#else /*#if PARALLEL*/
  sprintf( filename, "./in/%dx%d_proc%04d.bmp",
      get_LX(lattice), get_LY(lattice), get_proc_id(lattice));
#endif
  if( !( in = fopen( filename, "r")))
  {
    printf("%s %d >> u2bmp() -- Error opening file \"%s\".\n",
      __FILE__,__LINE__,filename);
    process_exit(1);
  }

  // n = fread( void *BUF, size_t SIZE, size_t COUNT, FILE *FP);

  n = fread( &bmfh, sizeof(struct bitmap_file_header), 1, in );
  if( strncmp(bmfh.bfType,"BM",2))
  {
    printf("ERROR: Can't process this file type.  Exiting!\n");
    printf("\n");
    process_exit(1);
  }
  n = fread( &bmih, sizeof(struct bitmap_info_header), 1, in );
  int_ptr = (int*)bmih.biCompression;
  if( *int_ptr != 0)
  {
    printf("ERROR: Can't handle compression.  Exiting!\n");
    printf("\n");
    process_exit(1);
  }

  width_ptr = (int*)bmih.biWidth;
  height_ptr = (int*)bmih.biHeight;
  bitcount_ptr = (short int*)bmih.biBitCount;

  // Read palette entries, if applicable.
  if( ENDIAN2(*bitcount_ptr) < 24)
  {
    n = (int)pow(2.,(real)ENDIAN2(*bitcount_ptr)); // Num palette entries.
    for( i=0; i<n; i++)
    {
      k = fread( &rgb, sizeof(struct rgb_quad), 1, in );
      if( k!=1)
      {
        printf("Error reading palette entry %d.  Exiting!\n", i);
        process_exit(1);
      }
    }
  }

  fclose(in);

  // Bytes per row of the bitmap.
  bytes_per_row =
    ((int)ceil(( (((real)(ENDIAN4(*width_ptr)))*((real)(ENDIAN2(*bitcount_ptr))))/8.)));

  // Bitmaps pad rows to preserve 4-byte boundaries.
  // The length of a row in the file will be bytes_per_row + pad .
  pad = ((4) - bytes_per_row%4)%4;

  compute_max_upr( lattice, max_u);

  sprintf( filename, "%s/upr%dx%d_frame%04d_proc%04d.bmp", get_out_path(lattice),
      get_LX(lattice),
      get_LY(lattice),
      frame, get_proc_id(lattice));
  if( !( o_u = fopen( filename, "w+")))
  {
    printf("ERROR: fopen( \"%s\", \"w+\") = NULL.  Bye, bye!\n", filename);
    process_exit(1);
  }

  sprintf( filename, "%s/upr_x%dx%d_frame%04d_proc%04d.bmp", get_out_path(lattice),
      get_LX(lattice),
      get_LY(lattice),
      frame, get_proc_id(lattice));
  if( !( o_ux = fopen( filename, "w+")))
  {
    printf("ERROR: fopen( \"%s\", \"w+\") = NULL.  Bye, bye!\n", filename);
    process_exit(1);
  }

  sprintf( filename, "%s/upr_y%dx%d_frame%04d_proc%04d.bmp", get_out_path(lattice),
      get_LX(lattice),
      get_LY(lattice),
      frame, get_proc_id(lattice));
  if( !( o_uy = fopen( filename, "w+")))
  {
    printf("ERROR: fopen( \"%s\", \"w+\") = NULL.  Bye, bye!\n", filename);
    process_exit(1);
  }

  fwrite( &bmfh, sizeof(struct bitmap_file_header), 1, o_u );
  fwrite( &bmih, sizeof(struct bitmap_info_header), 1, o_u );

  fwrite( &bmfh, sizeof(struct bitmap_file_header), 1, o_ux );
  fwrite( &bmih, sizeof(struct bitmap_info_header), 1, o_ux );

  fwrite( &bmfh, sizeof(struct bitmap_file_header), 1, o_uy);
  fwrite( &bmih, sizeof(struct bitmap_info_header), 1, o_uy);

  //for( j=get_LY(lattice)-1; j>=0; j--)
  for( j=0; j<get_LY(lattice); j++)
  {
    n = j*get_LX(lattice);

    for( i=0; i<get_LX(lattice); i++, n++)
    {
      if( lattice->bc[0][ n].bc_type == /*FLUID_NODE*/0)
      {
#if 1
        blue_val  = (char)0;
        green_val = (char)0;
        red_val   = (char)0;

        u_x = (lattice->upr[ n].u[0]);
        u_y = (lattice->upr[ n].u[1]);

        u = sqrt(u_x*u_x + u_y*u_y);
        maxu = sqrt( max_u[0]*max_u[0] + max_u[1]*max_u[1]);

        if( is_solid(lattice,n))
        {
          blue_val  = (char)ROUND( 128.*fabs(u_x)/max_u[0]);
          green_val = (char)ROUND( 128.*fabs(u_y)/max_u[1]);

        }

        else
        {
          blue_val  = (char)ROUND(  255.*fabs(u_x)/max_u[0]);
          green_val = (char)ROUND(  255.*fabs(u_y)/max_u[1]);
          //red_val    = (char)ROUND( 128.*fabs(u)/maxu);

        }
#else
        blue_val  = (char)255;
        green_val = (char)255;
        red_val   = (char)255;

        u = sqrt(u_x*u_x + u_y*u_y);
        maxu = sqrt( max_u[0]*max_u[0] + max_u[1]*max_u[1]);

        //if( fabs(u) > .1*maxu)
        //{
          green_val  = (char)ROUND( 255.-255.*fabs(u)/maxu);
          red_val    = (char)ROUND( 255.-255.*fabs(u)/maxu);
          blue_val    = (char)ROUND( 255.-255.*fabs(u)/maxu);
        //}
        //else
        //{
        //  green_val  = (char)0;
        //  red_val    = (char)0;
        //  blue_val    = (char)0;
        //}
#endif

        val = (char)0;

      } /* if( lattice->bc[subs][ n].bc_type == 0) */

      else // lattice->bc[subs][ n].bc_type != 0
      {
#if SOLID_COLOR_IS_CHECKERBOARD
        // Checkerboard pattern over the solids and boundary conditions.
        if( (i+j)%2)
        {
          red_val   = (char)200;
          green_val = (char)200;
          blue_val  = (char)200;
          val       = (char)200;
        }
        else
        {
          red_val   = (char)184;
          green_val = (char)184;
          blue_val  = (char)184;
          val       = (char)184;
        }
#else /* !( SOLID_COLOR_IS_CHECKERBOARD) */
#if SOLID_COLOR_IS_BLACK
        red_val   = (char)0;
        green_val = (char)0;
        blue_val  = (char)0;
        val       = (char)0;
#else /* !( SOLID_COLOR_IS_BLACK) */
        red_val   = (char)255;
        green_val = (char)255;
        blue_val  = (char)255;
        val       = (char)255;
#endif /* SOLID_COLOR_IS_BLACK */
#endif /* SOLID_COLOR_IS_CHECKERBOARD */

      } /* if( lattice->bc[subs][ n].bc_type == 0) else */

#if MARK_ORIGIN_FOR_REFERENCE
 // Mark the origin for reference.
 if( ( i == 0 && j == 0))
 {
   red_val   = (char)255;
   green_val = (char)255;
   blue_val  = (char)255;
   val       = (char)255;
 }
#endif /* MARK_ORIGIN_FOR_REFERENCE */

      if( fwrite( &blue_val,  1, 1, o_u ) != 1) { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &green_val, 1, 1, o_u ) != 1) { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &red_val,   1, 1, o_u ) != 1) { printf("BOOM!\n"); process_exit(1);}

      if( lattice->bc[0][ n].bc_type == /*FLUID_NODE*/0)
      {
        blue_val  = (char)0;
        green_val = (char)0;
        red_val   = (char)0;

        u_x = (lattice->upr[ n].u[0]);
        u_y = (lattice->upr[ n].u[1]);

        val = (char)0;

        if( is_solid(lattice,n))
        {
          if( u_x > 0)
          {
            red_val = val;
            blue_val = (char)ROUND( 128.*fabs(u_x)/max_u[0]);
          }
          else
          {
            red_val = (char)ROUND( 128.*fabs(u_x)/max_u[0]);
            blue_val = val;
          }

        }

        else
        {
          if( u_x > 0)
          {
            red_val = val;
            blue_val = (char)ROUND( 255.*fabs(u_x)/max_u[0]);
          }
          else
          {
            red_val = (char)ROUND( 255.*fabs(u_x)/max_u[0]);
            blue_val = val;
          }

        }

      } /* if( lattice->bc[subs][ n].bc_type == 0) */

      else // lattice->bc[subs][ n].bc_type != 0
      {
#if SOLID_COLOR_IS_CHECKERBOARD
        // Checkerboard pattern over the solids and boundary conditions.
        if( (i+j)%2)
        {
          red_val   = (char)200;
          green_val = (char)200;
          blue_val  = (char)200;
          val       = (char)200;
        }
        else
        {
          red_val   = (char)184;
          green_val = (char)184;
          blue_val  = (char)184;
          val       = (char)184;
        }
#else /* !( SOLID_COLOR_IS_CHECKERBOARD) */
#if SOLID_COLOR_IS_BLACK
        red_val   = (char)0;
        green_val = (char)0;
        blue_val  = (char)0;
        val       = (char)0;
#else /* !( SOLID_COLOR_IS_BLACK) */
        red_val   = (char)255;
        green_val = (char)255;
        blue_val  = (char)255;
        val       = (char)255;
#endif /* SOLID_COLOR_IS_BLACK */
#endif /* SOLID_COLOR_IS_CHECKERBOARD */

      } /* if( lattice->bc[subs][ n].bc_type == 0) else */

#if MARK_ORIGIN_FOR_REFERENCE
 // Mark the origin for reference.
 if( ( i == 0 && j == 0))
 {
   red_val   = (char)255;
   green_val = (char)255;
   blue_val  = (char)255;
   val       = (char)255;
 }
#endif /* MARK_ORIGIN_FOR_REFERENCE */

      if( fwrite( &blue_val,  1, 1, o_ux) != 1) { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &val,       1, 1, o_ux) != 1) { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &red_val,   1, 1, o_ux) != 1) { printf("BOOM!\n"); process_exit(1);}

      if( lattice->bc[0][ n].bc_type == /*FLUID_NODE*/0)
      {
        blue_val  = (char)0;
        green_val = (char)0;
        red_val   = (char)0;

        u_x = (lattice->upr[ n].u[0]);
        u_y = (lattice->upr[ n].u[1]);

        val = (char)0;

        if( is_solid(lattice,n))
        {
          if( u_y > 0)
          {
            red_val = val;
            green_val = (char)ROUND( 128.*fabs(u_y)/max_u[1]);
          }
          else
          {
            red_val = (char)ROUND( 128.*fabs(u_y)/max_u[1]);
            green_val = val;
          }

        }

        else
        {
          if( u_y > 0)
          {
            blue_val = (char)ROUND( 255.*fabs(u_y)/max_u[1]);
            green_val = val;
          }
          else
          {
            blue_val = val;
            green_val = (char)ROUND( 255.*fabs(u_y)/max_u[1]);
          }

        }

      } /* if( lattice->bc[subs][ n].bc_type == 0) */

      else // lattice->bc[subs][ n].bc_type != 0
      {
#if SOLID_COLOR_IS_CHECKERBOARD
        // Checkerboard pattern over the solids and boundary conditions.
        if( (i+j)%2)
        {
          red_val   = (char)200;
          green_val = (char)200;
          blue_val  = (char)200;
          val       = (char)200;
        }
        else
        {
          red_val   = (char)184;
          green_val = (char)184;
          blue_val  = (char)184;
          val       = (char)184;
        }
#else /* !( SOLID_COLOR_IS_CHECKERBOARD) */
#if SOLID_COLOR_IS_BLACK
        red_val   = (char)0;
        green_val = (char)0;
        blue_val  = (char)0;
        val       = (char)0;
#else /* !( SOLID_COLOR_IS_BLACK) */
        red_val   = (char)255;
        green_val = (char)255;
        blue_val  = (char)255;
        val       = (char)255;
#endif /* SOLID_COLOR_IS_BLACK */
#endif /* SOLID_COLOR_IS_CHECKERBOARD */

      } /* if( lattice->bc[subs][ n].bc_type == 0) else */

#if MARK_ORIGIN_FOR_REFERENCE
 // Mark the origin for reference.
 if( ( i == 0 && j == 0))
 {
   red_val   = (char)255;
   green_val = (char)255;
   blue_val  = (char)255;
   val       = (char)255;
 }
#endif /* MARK_ORIGIN_FOR_REFERENCE */

      if( fwrite( &blue_val,  1, 1, o_uy) != 1) { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &green_val, 1, 1, o_uy) != 1) { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &red_val,   1, 1, o_uy) != 1) { printf("BOOM!\n"); process_exit(1);}

    } /* for( i=0; i<get_LY(lattice); i++) */

    // Pad for 4-byte boundaries.
    val = (char)0;
    for( i=0; i<pad; i++)
    {
      if( fwrite( &val, 1, 1, o_u ) != 1) { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &val, 1, 1, o_ux) != 1) { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &val, 1, 1, o_uy) != 1) { printf("BOOM!\n"); process_exit(1);}
    }

  } /* for( j=0; j<get_LY(lattice); j++) */

  fclose(o_u );
  fclose(o_ux);
  fclose(o_uy);

#if VERBOSITY_LEVEL > 0
  sprintf( filename, "%s/upr%dx%d_frame%04d_proc%04d.bmp", get_out_path(lattice),
      get_LX(lattice), get_LY(lattice), frame, get_proc_id(lattice));
  printf("u2bmp()   -- Wrote file \"%s\".\n", filename);

  sprintf( filename, "%s/upr_x%dx%d_frame%04d_proc%04d.bmp", get_out_path(lattice),
      get_LX(lattice), get_LY(lattice), frame, get_proc_id(lattice));
  printf("u2bmp()   -- Wrote file \"%s\".\n", filename);

  sprintf( filename, "%s/upr_y%dx%d_frame%04d_proc%04d.bmp", get_out_path(lattice),
      get_LX(lattice), get_LY(lattice), frame, get_proc_id(lattice));
  printf("u2bmp()   -- Wrote file \"%s\".\n", filename);
#endif /* VERBOSITY_LEVEL > 0 */
#endif /* STORE_U_COMPOSITE */


#if SAY_HI
  printf("u2bmp() -- Bye!\n");
  printf("\n");
#endif /* SAY_HI */

} /* u2bmp( lattice_ptr lattice, int time) */
#endif



/*void write_raw1(
       lattice_ptr lattice,
       real *a,
       int     stride,
       real  a_max,
       real  a_min,
       char   *filename)
{
  int size;
  int n, ni, nj, nk;
  unsigned char *Xpix1;
  real half;
  FILE *infile;
//  half = (a_max+ a_min)/2.;
    half = lattice->param.rho_A[0]-0.5;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;
//       R      G       B
//14  88  0  0
//15  114  0  1
//16  140  0  1
//17  166  0  1        decane
//18  192  0  1
//19  218  0  1
//20  244  0  1
//
//33  255  0  2
//34  255  0  78
//35  255  0  171
//36  127  0  254
//37  0  0  254
//38  0  0  254      water
//
//64  0  0  254
//65  0  0  178
//66  88  36  101
//67  213  198  24
//68  254  254  0       solid

  size = ni*nj*nk*sizeof(unsigned char);

  if(!(Xpix1 = ( unsigned char*)malloc( size)))
  {
    printf("%s %d %04d >> ERROR: Can't malloc Xpix1. (Exiting!)\n",
      __FILE__,__LINE__, get_proc_id(lattice));
    process_exit(1);
  }
  for( n=0; n<ni*nj*nk; n++)
  {
    if( !lattice->solids[0][n].is_solid && a[stride*n]<half)
    {
      Xpix1[n] = (unsigned char)ROUND(64.-29.*(a[stride*n]-a_min)/(half-a_min));
//According to Utah color Table, please refer to Utah Color Table!  35--64
    }
    else if(!lattice->solids[0][n].is_solid && a[stride*n]>half)
    {
      Xpix1[n] = (unsigned char)ROUND(35.-15.*(a[stride*n]-half)/(a_max-half));
//According to Utah color Table, please refer to Utah Color Table!  20---35
    }

    else
    {
      Xpix1[n] = (unsigned char)255.;
    }
    //printf("%s %d %04d >> write_raw() -- Xpix1[%d] = 255*ROUND(%f) = %d.\n",
    //    __FILE__, __LINE__, get_proc_id(lattice),
    //  n, (a[stride*n]-a_min)/(a_max-a_min), (int)Xpix1[n]);
  }

  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fwrite( (char*)Xpix1, 1, size, infile);

  fclose(infile);

  printf("%s %d %04d >> write_raw() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);
} // void write_raw( lattice_ptr lattice, real *a,  ... */
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void write_dat(
       lattice_ptr lattice,
       real *a,
       int     stride,
       real  a_max,
       real  a_min,
       char   *filename)
{
  int size;
  int n, ni, nj, nk;
  real*Xpix1;
  FILE *infile ;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  size = ni*nj*nk*sizeof(real);

  if(!(Xpix1 = ( real*)malloc( size)))
  {
    printf("%s %d %04d >> ERROR: Can't malloc Xpix1. (Exiting!)\n",
      __FILE__,__LINE__, get_proc_id(lattice));
    process_exit(1);
  }

  for( n=0; n<ni*nj*nk; n++)
  {
    if( !is_solid( lattice, n))
    {
      Xpix1[n] = a[stride*n];
    }
    else
    {
      Xpix1[n] = 0.;
    }
  }

  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }

  fwrite( (real*)Xpix1, 1, size, infile);

  fclose(infile);

  printf("%s %d %04d >> write_dat() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);

} /* void write_dat( lattice_ptr lattice, real *a,  ... */

#if 0
// TODO: Revise this section for the new data structures or delete it.
#if (NUM_FLUID_COMPONENTS==2)
void write_plt(
       lattice_ptr lattice,
       real *a, real *b,
       char   *filename, char *filename2)
{
  int size;
  int n, ni, nj, nk, count, basek;
  FILE *infile, *infile2, *FL1, *FL2, *FL3, *FL4;
  real v_in, v_out, v_inupr, v_outupr;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  basek = 0;
#if PARALLEL
  basek = get_g_SZ( lattice);
#endif

  size = ni*nj*nk*sizeof(real);


  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  if (!(infile2 = fopen(filename2,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }


  //  fprintf( infile, "variables = x, y, z, rho1, rho2, solid\n"  );
  fprintf( infile, "variables =rho1\n"  );
  fprintf( infile, "zone i=%d, j=%d, k=%d,  f=point\n", ni,  nj, nk );
  fprintf( infile2, "variables =rho2\n"  );
  fprintf( infile2, "zone i=%d, j=%d, k=%d,  f=point\n", ni,  nj, nk );

  count = 0;

  if(get_proc_id(lattice) == get_num_procs(lattice)-1)
  {

  if (!(FL1 = fopen("./out/in_vels.dat","w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  if (!(FL3 = fopen("./out/zupper.plt","w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  //ptinf out the 2D plt file of upper inlet surface
  fprintf(FL3,"variables = i, j, rho1, u, v, w\n");
  fprintf(FL3,"zone i=%d, j=%d,  f=point\n", ni,  nj);
   v_in = 0.;
   v_inupr = 0.;
  }
  // endif get_proc_id(lattice) == get_num_procs(lattice)-1)
  if(get_proc_id(lattice) == 0)
  {
  if (!(FL2 = fopen("./out/out_vels.dat","w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }

  if (!(FL4 = fopen("./out/zbottom.plt","w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fprintf(FL4,"variables = i, j, rho1, u, v, w\n");
  fprintf(FL4,"zone i=%d, j=%d,  f=point\n", ni,  nj);

   v_out = 0.;
   v_outupr = 0.;
  }
  //endif (get_proc_id(lattice) == 0

  for( n=0; n<ni*nj*nk; n++)
  {
     if( is_solid( lattice, n))
     {
      a[4*n] = -1.;
      get_rho( lattice, 1, n);
     }

     if(  get_proc_id(lattice) == get_num_procs(lattice)-1 &&
N2Z(n,ni,nj,nk)+basek == (get_num_procs(lattice) *lattice->param.LZ-1) )
     {
       v_in= v_in + a[4*n+3];
       v_inupr= v_inupr +lattice->ueq[n].u[2];

    fprintf(FL3,"%4d %4d %7.4f  %2.9f  %2.9f %2.9f \n",N2X(n,ni,nj,nk),N2Y(n,ni,nj,nk), a[4*n], a[4*n+1], a[4*n+2], a[4*n+3] ); //MS changed format from 8.5 to 2.9
     }
     if(get_proc_id(lattice) == 0 && N2Z(n,ni,nj,nk)+basek == 1)
     {
       v_out= v_out + a[4*n+3];
             v_outupr= v_outupr +lattice->ueq[n].u[2];
    fprintf(FL4,"%4d %4d %7.4f  %2.9f  %2.9f %2.9f \n",N2X(n,ni,nj,nk),N2Y(n,ni,nj,nk), a[4*n], a[4*n+1], a[4*n+2], a[4*n+3] ); //MS changed format from 8.5 to 2.9

     }
     if( (N2Z(n,ni,nj,nk)+basek >lattice->param.z1 &&
    N2Z(n,ni,nj,nk)+basek <lattice->param.z2
           ) && (a[4*n] > (lattice->param.rho_A[0]+lattice->param.rho_A[1])/2.) )
     {count++;}
//  fprintf( infile, "%4d  %4d %4d  %7.4f  %7.4f  %2d\n", N2X(n,ni,nj,nk),N2Y(n,ni,nj,nk), N2Z(n,ni,nj,nk),
//  a[/*stride*/ 4*n], lattice->macro_vars[1][n].rho , lattice->solids[0][n].is_solid/255  );
  fprintf( infile, " %7.4f \n",   a[/*stride*/ 4*n]  );
  fprintf( infile2, " %7.4f \n",   b[/*stride*/ 4*n]  );
  }
  if(get_proc_id(lattice) == get_num_procs(lattice)-1)
  {
  fprintf(FL1,"position=%4d",(get_num_procs(lattice) *lattice->param.LZ-2 ));
  fprintf(FL1, "integration in =%12.7f\n", v_in );
  fprintf(FL1, "integration in upr=%12.7f\n", v_inupr );
  fclose(FL1);
  fclose(FL3);
  }
  if(get_proc_id(lattice) == 0 )
  {
  fprintf(FL2, "integration out=%12.7f\n", v_out );
  fprintf(FL2, "integration out upr =%12.7f\n", v_outupr );
  fclose(FL2);
  fclose(FL4);
  }

  fclose(infile);
  fclose(infile2);

  printf("%s %d %04d >> write_tecplot_file() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);
  printf("%s %d >> Non-Wetting fluid nodes: %d \n",
      __FILE__, __LINE__, count);

} /* void write_plt( lattice_ptr lattice, real *a,  ... */
#endif

void write_plt_uvw(
       lattice_ptr lattice,
       real *a, real *b,
       char   *filename)
{
  int size;
  int n, ni, nj, nk;
  FILE *infile;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  size = ni*nj*nk*sizeof(real);

  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fprintf( infile, "variables = x, y, z, u, v, w\n"  );
  fprintf( infile, "zone i=%d, j=%d, k=%d,  f=point\n", ni,  nj, nk );


 for( n=0; n<ni*nj*nk; n++)
  {
   if( lattice->solids[0][n].is_solid)
     {
      a[4*n] = -1.;
      lattice->macro_vars[1][n].rho = -1.;
     }

//   if(NUM_FLUID_COMPONENTS ==2)
   {
  fprintf( infile, "%4d  %4d %4d  %7.4f  %7.4f  %7.4f\n", N2X(n,ni,nj,nk),N2Y(n,ni,nj,nk), N2Z(n,ni,nj,nk),
    b[/*stride*/ 4*n+1], b[/*stride*/ 4*n+2], b[/*stride*/ 4*n+3]  );
   }

  }

  fclose(infile);

  printf("%s %d %04d >> write_tecplot_uvw_file() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);

} /* void write_plt_uvw( lattice_ptr lattice, real *a,  ... */


void write_plt_single(
       lattice_ptr lattice,
       real *a,
       char   *filename)
{
  int size;
  int n, ni, nj, nk, basek;
  real v_in, v_out;
  FILE *infile, *FL1, *FL2;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;
  basek = 0;
#if PARALLEL
  basek = get_g_SZ( lattice);
#endif

  size = ni*nj*nk*sizeof(real);

  if(lattice->time==0 && get_proc_id(lattice) == get_num_procs(lattice)-1)
  {
  if (!(FL1 = fopen("./out/in_vels.dat","w+")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fprintf(FL1, "time  integration in \n" );

  }//endif   if(get_proc_id(lattice) == get_num_procs(lattice)-1)
  if(lattice->time!=0 && get_proc_id(lattice) == get_num_procs(lattice)-1)
  {
  if (!(FL1 = fopen("./out/in_vels.dat","a+")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  }//endif   if(get_proc_id(lattice) == get_num_procs(lattice)-1)

  if (lattice->time==0 && get_proc_id(lattice) == 0)
  {
  if (!(FL2 = fopen("./out/out_vels.dat","w+")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fprintf(FL2, "time  integration in \n" );

  }
//if   if(get_proc_id(lattice) == 0
  if (lattice->time!=0 && get_proc_id(lattice) == 0)
  {
  if (!(FL2 = fopen("./out/out_vels.dat","a+")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  }
//if   if(get_proc_id(lattice) == 0

  if (!(infile = fopen(filename,"w")))
  {
    printf("Can't creat file\n");
    process_exit(-1);
  }
  fprintf( infile, "variables = x, y, z, rho1, u, v, w, solid\n"  );
//  fprintf( infile, "variables = rho1, w, solid\n"  );
  fprintf( infile, "zone i=%d, j=%d, k=%d,  f=point\n", ni,  nj, nk );

  v_in = 0.;
  v_out = 0.;
 for( n=0; n<ni*nj*nk; n++)
  {
   if( lattice->solids[0][n].is_solid)
     {
      a[4*n] = -1.;
     }
   else
     {
       if(  (get_proc_id(lattice) == get_num_procs(lattice)-1) &&
           N2Z(n,ni,nj,nk)+basek == get_num_procs(lattice) *lattice->param.LZ-1) v_in= v_in + a[4*n+3];
       if(  (get_proc_id(lattice) == 0) && N2Z(n,ni,nj,nk)+basek ==1) v_out= v_out + a[4*n+3];
     }
  fprintf( infile, "%4d  %4d %4d  %7.4f  %2.9f %2.9f %2.9f %2d\n", N2X(n,ni,nj,nk),N2Y(n,ni,nj,nk), N2Z(n,ni,nj,nk)+basek,
  a[/*stride*/ 4*n], a[/*stride*/ 4*n+1], a[/*stride*/ 4*n+2], a[/*stride*/ 4*n+3], lattice->solids[0][n].is_solid/255  ); //MS changed format from 8.5 to 2.9

//  fprintf( infile, " %7.4f  %7.4f  %2d\n",
//  a[/*stride*/ 4*n],  a[/*stride*/ 4*n+3], lattice->solids[0][n].is_solid/255  );
  }

if(get_proc_id(lattice) == get_num_procs(lattice)-1) //Works only if Parallel (MS)

  {
  fprintf(FL1, "%8d  %12.7f\n", lattice->time, v_in );
  fclose(FL1);
  }

  if(get_proc_id(lattice) == 0)
  {
  fprintf(FL2, "%8d  %12.7f\n", lattice->time, v_out );
  fclose(FL2);
  }

  fclose(infile);

  printf("%s %d %04d >> write_tecplot_file() -- Wrote to file \"%s\".\n",
      __FILE__, __LINE__, get_proc_id(lattice), filename);

} /* void write_plt( lattice_ptr lattice, real *a,  ... */
#endif

void write_rho_txt( lattice_ptr lattice)
{
  int i, j, k;
  int ni = get_ni(lattice);
  int nj = get_nj(lattice);
  int nk = get_nk(lattice);
  int n;
  int subs;
  char filename[1024];

  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
#if 0
    if( get_NumDims(lattice)==2)
    {
      sprintf( filename
             , "%s/rho%dx%d_frame%04d_subs%02d_proc%04d.txt"
             , get_out_path(lattice)
             , ni, nj
             , get_frame(lattice)
             , subs
             , get_proc_id(lattice) );
    }
    else
    {
      sprintf( filename
             , "%s/rho%dx%dx%d_frame%04d_subs%02d_proc%04d.txt"
             , get_out_path(lattice)
             , ni, nj, nk
             , get_frame(lattice)
             , subs
             , get_proc_id(lattice) );
    }
#else
    gen_filename( lattice, filename, "rho", get_frame(lattice), subs, ".txt");
#endif

    FILE* fout;
    fout = fopen(filename,"w+");

    n = 0;
    for( k=0; k<nk; k++)
    {
      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          fprintf(fout," %18.16f",get_rho(lattice,subs,n));
          n++;
        }
        fprintf(fout,"\n");
      }
      fprintf(fout,"\n");
    }

    fclose(fout);
  }
}

// S P Y   B M P  {{{
//##############################################################################
// void spy_bmp( char *filename, int ***spy)
//
//  - Returns matrix 'spy' of ones and zeros.
//
//  - Zeros for white pixels.
//
//  - Ones for non-white pixels.
//
void spy_bmp( lattice_ptr lattice, char *filename, int **matrix)
{
  FILE *in, *o;
  int i, j, n, m;
  int g_i, g_j;
  int pad, bytes_per_row;
  char k;
  char b, g, r;
  struct bitmap_file_header bmfh;
  struct bitmap_info_header bmih;
  struct rgb_quad rgb;
  int *int_ptr;
  short int *short_int_ptr;
  int *width_ptr;
  int *height_ptr;
  short int *bitcount_ptr;
  char ctemp;
  int  itemp;
  int **spy;

  printf("spy_bmp() -- Hi!\n");

  spy = (int**)malloc( get_g_LY(lattice)*sizeof(int*));
  for( j=0; j<get_g_LY(lattice); j++)
  {
    spy[j] = (int*)malloc( get_g_LX(lattice)*sizeof(int));
  }

  // Clear the spy array.
  for( j=0; j<get_g_LY(lattice); j++)
  {
    for( i=0; i<get_g_LX(lattice); i++)
    {
      spy[j][i] = 0;
    }
  }

  // Clear the matrix array.
  for( j=0; j<get_LY(lattice); j++)
  {
    for( i=0; i<get_LX(lattice); i++)
    {
      matrix[j][i] = 0;
    }
  }

  if(/*ignore_solids*/0) return;

  if( !( in = fopen( filename, "r")))
  {
#if 1
    printf("%s %d >> spy_bmp() -- Error opening file \"%s\".\n",
      __FILE__, __LINE__, filename);
    process_exit(1);
#else
    printf(" %s::spy_bmp() %d >> File \"%s\" cannot be opened for reading.\n",
        __FILE__, __LINE__, filename);
    if( !( o = fopen( filename, "w+")))
    {
      // TODO: Write blank bmp file.
    }
    printf(" %s::spy_bmp() %d >> Wrote a blank \"%s\" file.\n",
        __FILE__, __LINE__, filename);
    printf(" %s::spy_bmp() %d >> Returning all zeros!\n", __FILE__, __LINE__);
    fclose( o);
    return;
#endif
  }

  // n = fread( void *BUF, size_t SIZE, size_t COUNT, FILE *FP);

  n = fread( &bmfh, sizeof(struct bitmap_file_header), 1, in );
  if( strncmp(bmfh.bfType,"BM",2))
  {
    printf("ERROR: Can't process this file type.  Exiting!\n");
    printf("\n");
    process_exit(1);
  }
  n = fread( &bmih, sizeof(struct bitmap_info_header), 1, in );

  *((int*)(bmih.biWidth)) = ENDIAN4(((int)(*((int*)(bmih.biWidth)))));
  *((int*)(bmih.biHeight)) = ENDIAN4(((int)(*((int*)(bmih.biHeight)))));
  *((short int*)(bmih.biBitCount)) =
    ENDIAN2(((short int)(*((short int*)(bmih.biBitCount)))));

  int_ptr = (int*)bmih.biCompression;
  if( *int_ptr != 0)
  {
    printf("%s %d >> ERROR: Can't handle compression.  Exiting!\n"
          ,__FILE__,__LINE__);
    printf("\n");
    process_exit(1);
  }

  width_ptr = (int*)bmih.biWidth;
  height_ptr = (int*)bmih.biHeight;
  bitcount_ptr = (short int*)bmih.biBitCount;

  if( *width_ptr != get_g_LX(lattice))
  {
    printf("%s %d >> ERROR: LX %d does not match the "
        "width %d of the BMP file. Exiting!\n"
        "Note that, if the width stated here seems absurd, you\n"
        "might need to recompile with the SWAP_BYTE_ORDER flag.\n"
        "This can be done by \"make swap\".\n",
        __FILE__, __LINE__, get_g_LX(lattice), *width_ptr);
    process_exit(1);
  }
  printf("%s %d >> biWidth = %d \n", __FILE__, __LINE__, (int)*bmih.biWidth);
  printf("%s %d >> width_ptr = %d \n", __FILE__, __LINE__, (int)*width_ptr);

  if( *height_ptr != get_g_LY(lattice))
  {
    printf("%s %d >> ERROR: LY %d does not match the "
        "height %d of the BMP file. Exiting!\n",
        __FILE__, __LINE__, get_g_LY(lattice), *height_ptr);
    process_exit(1);
  }

  if( (*bitcount_ptr) < 24)
  {
    n = (int)pow(2.,(double)(*bitcount_ptr)); // Num palette entries.
    for( i=0; i<n; i++)
    {
      k = fread( &rgb, sizeof(struct rgb_quad), 1, in );
      if( k!=1)
      {
        printf("%s %d >> Error reading palette entry %d.  Exiting!\n"
              , __FILE__, __LINE__, i);
        process_exit(1);
      }
    }
  }

  // Bytes per row of the bitmap.
  bytes_per_row =
    ((int)ceil(( (((double)(*width_ptr))*((double)((*bitcount_ptr))))/8.)));

  // Bitmaps pad rows to preserve 4-byte boundaries.
  // The length of a row in the file will be bytes_per_row + pad .
  pad = ((4) - bytes_per_row%4)%4;

  n = 0;
  m = 0;
  n+=( k = fread( &b, 1, 1, in ));
  i = 0;
  j = 0;
  while( !feof(in))
  {
    switch((*bitcount_ptr))
    {
      case 1: // Monochrome.
        printf("%s %d >> spy_bmp() -- "
            "Support for Monochrome BMPs is pending.  "
            "Exiting!\n", __FILE__, __LINE__);
        process_exit(1);

        if( i < *width_ptr) { (spy)[j][i] = ( (b & 0x80) == 0); }
        i++;
        if( i < *width_ptr) { (spy)[j][i] = ( (b & 0x40) == 0); }
        i++;
        if( i < *width_ptr) { (spy)[j][i] = ( (b & 0x20) == 0); }
        i++;
        if( i < *width_ptr) { (spy)[j][i] = ( (b & 0x10) == 0); }
        i++;
        if( i < *width_ptr) { (spy)[j][i] = ( (b & 0x08) == 0); }
        i++;
        if( i < *width_ptr) { (spy)[j][i] = ( (b & 0x04) == 0); }
        i++;
        if( i < *width_ptr) { (spy)[j][i] = ( (b & 0x02) == 0); }
        i++;
        if( i < *width_ptr) { (spy)[j][i] = ( (b & 0x01) == 0); }
        i++;
        break;

      case 4: // 16 colors.
        printf("%s %d >> spy_bmp() -- "
            "Support for 16 color BMPs is pending.  "
            "Exiting!\n", __FILE__, __LINE__);
        process_exit(1);

        if( i < *width_ptr) { (spy)[j][i] = ( (b&0xf0)>>4 != 15); }
        i++;
        if( i < *width_ptr) { (spy)[j][i] = ( (b&0x0f) != 15); }
        i++;
        break;

      case 8: // 256 colors.
        printf("%s %d >> spy_bmp() -- "
            "Support for 256 color BMPs is pending.  "
            "Exiting!\n", __FILE__, __LINE__);
        process_exit(1);

        if( i < *width_ptr) { (spy)[j][i] = ( (b&0xff) != 255); }
        i++;
        break;

      case 24: // 24-bit colors.
        if( i < 3*(*width_ptr))
        {
          i++; n+=( k = fread( &g, 1, 1, in ));
          i++; n+=( k = fread( &r, 1, 1, in ));

          if( ( (b&0xff) == 0) &&( (g&0xff) == 0) &&( (r&0xff) == 0) )
          {
            (spy)[j][(int)floor((double)i/3.)] = 1;
          }

#if 0
          if( ( (b&0xff) == 0) &&( (g&0xff) == 0) &&( (r&0xff) == 255) )
          {
            // Red ==> Inflow, Pressure boundaries.
            if(    (int)floor((double)i/3.) == 0
                || (int)floor((double)i/3.) == get_g_LX(lattice)-1 )
            {
              if( !( j==0 || j == get_g_LY(lattice)-1))
              {
                lattice->periodic_x[subs] = 0;
              }
            }
            if(    j == 0
                || j == get_g_LY(lattice)-1 )
            {
              if( !(   (int)floor((double)i/3.) == 0
                    || (int)floor((double)i/3.) == get_g_LX(lattice)-1))
              {
                lattice->periodic_y[subs] = 0;
              }
            }
          }

          if( ( (b&0xff) == 0) &&( (g&0xff) == 255) &&( (r&0xff) == 0) )
          {
            // Green ==> Outflow, Pressure boundaries.
            if(    (int)floor((double)i/3.) == 0
                || (int)floor((double)i/3.) == get_g_LX(lattice)-1 )
            {
              if( !( j==0 || j == get_g_LY(lattice)-1))
              {
                lattice->periodic_x[subs] = 0;
              }
            }
            if(    j == 0
                || j == get_g_LY(lattice)-1 )
            {
              if( !(   (int)floor((double)i/3.) == 0
                    || (int)floor((double)i/3.) == get_g_LX(lattice)-1))
              {
                lattice->periodic_y[subs] = 0;
              }
            }
          }
#endif
        }
        i++;
        break;

      default: // 32-bit colors?
        printf("%s %d >> ERROR: Unhandled color depth, "
            "BitCount = %d. Exiting!\n", __FILE__, __LINE__, *bitcount_ptr);
        process_exit(1);
        break;

    } /* switch(*(bmih.biBitCount)) */

    if( !(n%(bytes_per_row+pad))) { m++; i=0; j++;}
    n+=( k = fread( &b, 1, 1, in ));

  } /* while( !feof(in)) */

  if( (bytes_per_row+pad)*m!=n)
  {
    printf("WARNING: Num bytes read = %d versus num bytes predicted = %d .\n",
        n, (bytes_per_row+pad)*m);
  }

  if( m != *height_ptr)
  {
    printf("WARNING: m (%d) != bmih.biHeight (%d).\n", m, *height_ptr);
  }

  fclose(in);

  for( j=0, g_j=get_g_SY(lattice); j<get_LY(lattice); j++, g_j++)
  {
    for( i=0, g_i=get_g_SX(lattice); i<get_LX(lattice); i++, g_i++)
    {
      matrix[j][i] = spy[g_j][g_i];
    }
  }

  for( j=0; j<get_g_LY(lattice); j++)
  {
    free(spy[j]);
  }
  free(spy);

  printf("spy_bmp() -- Bye!\n");
  printf("\n");

} /* spy_bmp( char *filename, int **spy) */
                                                                        // }}}


