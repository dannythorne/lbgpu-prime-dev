//##############################################################################
//
// lattice.c
//

int get_LX( lattice_ptr lattice) { return lattice->param.LX;}
int get_LY( lattice_ptr lattice) { return lattice->param.LY;}
int get_LZ( lattice_ptr lattice) { return lattice->param.LZ;}
int get_ni( lattice_ptr lattice) { return lattice->param.LX;}
int get_nj( lattice_ptr lattice) { return lattice->param.LY;}
int get_nk( lattice_ptr lattice) { return lattice->param.LZ;}
int* get_ni_ptr( lattice_ptr lattice) { return &(lattice->param.LX);}
int* get_nj_ptr( lattice_ptr lattice) { return &(lattice->param.LY);}
int* get_nk_ptr( lattice_ptr lattice) { return &(lattice->param.LZ);}

void set_LX( lattice_ptr lattice, const int arg_LX)
{
  lattice->param.LX = arg_LX;
}
void set_LY( lattice_ptr lattice, const int arg_LY)
{
  lattice->param.LY = arg_LY;
}
void set_LZ( lattice_ptr lattice, const int arg_LZ)
{
  lattice->param.LZ = arg_LZ;
}

int get_FrameRate( lattice_ptr lattice)
{
  return lattice->param.FrameRate;
}

int get_frame( lattice_ptr lattice)
{
  return lattice->frame;
}

int get_NumNodes( lattice_ptr lattice) { return lattice->NumNodes;}
int* get_NumNodes_ptr( lattice_ptr lattice) { return &(lattice->NumNodes);}
void set_NumNodes( lattice_ptr lattice)
{
  lattice->NumNodes = get_LX( lattice)*get_LY( lattice)*get_LZ( lattice);
}

int set_NumDims( lattice_ptr lattice, int value)
{
  lattice->NumDims = value;
}

int get_NumDims( lattice_ptr lattice)
{
  return lattice->NumDims;
}
int* get_NumDims_ptr( lattice_ptr lattice)
{
  return &(lattice->NumDims);
}

int set_NumVelDirs( lattice_ptr lattice, int value)
{
  lattice->NumVelDirs = value;
}

int get_NumVelDirs( lattice_ptr lattice)
{
  return lattice->NumVelDirs;
}
int* get_NumVelDirs_ptr( lattice_ptr lattice)
{
  return &(lattice->NumVelDirs);
}

int set_NumSubs( lattice_ptr lattice, int value)
{
  lattice->NumSubs = value;
}

int get_NumSubs( lattice_ptr lattice)
{
  return lattice->NumSubs;
}

int* get_NumSubs_ptr( lattice_ptr lattice)
{
  return &(lattice->NumSubs);
}


int is_solid( lattice_ptr lattice, const int n)
{
  if( lattice->solids_memblock[n])
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

void set_is_solid(
  lattice_ptr lattice
, const int n
, const unsigned char val)
{
  lattice->solids_memblock[n] = val;
}

int is_not_solid( lattice_ptr lattice, const int n)
{
  return !is_solid( lattice, n);
}

real get_rho( lattice_ptr lattice, const int subs, const int n)
{
  return lattice->vars[subs].macrovars1d[0][n];
}

real* get_rho_ptr( lattice_ptr lattice, const int subs, const int n)
{
  return lattice->vars[subs].macrovars1d[0] + n;
}

real get_ux( lattice_ptr lattice, const int subs, const int n)
{
  return lattice->vars[subs].macrovars1d[1][n];
}

real get_uy( lattice_ptr lattice, const int subs, const int n)
{
  return lattice->vars[subs].macrovars1d[2][n];
}

real get_uz( lattice_ptr lattice, const int subs, const int n)
{
  if( get_NumDims(lattice)==2)
  {
    return 0.;
  }
  else
  {
    return lattice->vars[subs].macrovars1d[3][n];
  }
}

real* get_ux_ptr( lattice_ptr lattice, const int subs, const int n)
{
  //return &(lattice->vars[subs].macrovars1d[1][n]);
  return lattice->vars[subs].macrovars1d[1] + n;
}

real* get_uy_ptr( lattice_ptr lattice, const int subs, const int n)
{
  //return &(lattice->vars[subs].macrovars1d[2][n]);
  return lattice->vars[subs].macrovars1d[2] + n;
}

real* get_uz_ptr( lattice_ptr lattice, const int subs, const int n)
{
  //return (get_NumDims(lattice)==2)?(NULL):(&(lattice->vars[subs].macrovars1d[3][n]));
  return
    (get_NumDims(lattice)==2)
    ?(NULL)
    :((lattice->vars[subs].macrovars1d[3] + n));
}

real get_max_ux( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.max_macrovars[subs][1];
}

real get_max_uy( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.max_macrovars[subs][2];
}

real get_max_uz( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.max_macrovars[subs][3];
}

real get_min_ux( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.min_macrovars[subs][1];
}

real get_min_uy( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.min_macrovars[subs][2];
}

real get_min_uz( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.min_macrovars[subs][3];
}

real get_ave_ux( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.ave_macrovars[subs][1];
}

real get_ave_uy( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.ave_macrovars[subs][2];
}

real get_ave_uz( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.ave_macrovars[subs][3];
}

real get_max_rho( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.max_macrovars[subs][0];
}

real get_min_rho( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.min_macrovars[subs][0];
}

real get_ave_rho( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.ave_macrovars[subs][0];
}

real get_flux_u( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.flux[subs][0];
}

real get_flux_x( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.flux[subs][1];
}

real get_flux_y( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.flux[subs][2];
}

real get_flux_z( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.flux[subs][3];
}

const char* get_out_path( lattice_ptr lattice)
{
  return "./out";
}

int do_post_streaming_bcs( lattice_ptr lattice)
{
  return !lattice->param.GZL
      && !lattice->param.AllBoundaryPeriodic;
}

int do_post_collision_bcs( lattice_ptr lattice)
{
  return  lattice->param.GZL
      && !lattice->param.AllBoundaryPeriodic;
}

void set_tic( lattice_ptr lattice, real t)
{
  lattice->tic =  t;
}

void set_toc( lattice_ptr lattice, real t)
{
  lattice->toc =  t;
}

real display_etime( lattice_ptr lattice)
{
  if( is_on_root_proc(lattice))
  {
    real t = lattice->toc - lattice->tic;
    printf("%s %d %04d >> elapsed time = "
      "%f seconds (%f minutes, %f hours, %f days)\n",
      __FILE__,__LINE__,get_proc_id(lattice)
    , t
    , t / 60.0
    , t / 60.0 / 60.0
    , t / 60.0 / 60.0 / 24.0);
  }
}

// THORNE 20111103
// {
void gen_filename(
  lattice_ptr lattice
, char* filename
, const char* prefix
, int frame
, int subs
, const char* suffix )
{
  char size_str[32];
  char frame_str[32];
  char subs_str[32];

  if( get_NumDims(lattice)==2)
  {
    sprintf(size_str,"%dx%d",get_g_LX(lattice),get_g_LY(lattice));
  }
  else
  {
    sprintf(
      size_str
    , "%dx%dx%d"
    , get_g_LX(lattice)
    , get_g_LY(lattice)
    , get_g_LY(lattice) );
  }

  if( frame >= 0)
  {
    sprintf( frame_str, "_frame%04d", frame);
  }
  else
  {
    frame_str[0] = '\0';
  }

  if( get_NumSubs(lattice) == 1)
  {
    subs_str[0] = '\0';
  }
  else
  {
    sprintf( subs_str, "_subs%02d", subs);
  }

  if( get_num_procs( lattice) > 1)
  {
    sprintf(
      filename
    ,"%s/%s%s%s%s_proc%04d%s"
    , get_out_path( lattice)
    , prefix
    , size_str
    , frame_str
    , subs_str
    , get_proc_id( lattice)
    , suffix
    );
  }
  else
  {
    sprintf(
      filename
    ,"%s/%s%s%s%s%s"
    , get_out_path( lattice)
    , prefix
    , size_str
    , frame_str
    , subs_str
    , suffix
    );
  }
}
// }



#ifdef __CUDACC__
// The following definitions of rho2bmp and u2bmp were relocated from lbio.c
// for nvcc so that they will be defined when write_rho_image and write_u_image
// are defined. This is due to the ad hoc way we are compiling the multiple
// source code files and the way that nvcc works.

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


   sprintf( filename
          , "%s/rho%dx%d_frame%04d_subs%02d_proc%04d.bmp"
          , get_out_path(lattice)
          , get_LX(lattice)
          , get_LY(lattice)
          , frame
          , subs
          , get_proc_id(lattice) );
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

  sprintf( filename, "%s/u%dx%d_frame%04d_subs%02d_proc%04d.bmp"
         , get_out_path(lattice)
         , get_LX(lattice)
         , get_LY(lattice)
         , frame, subs
         , get_proc_id(lattice));
  if( !( o_u = fopen( filename, "w+")))
  {
    printf("ERROR: fopen( \"%s\", \"w+\") = NULL.  Bye, bye!\n", filename);
    process_exit(1);
  }

  sprintf( filename
         , "%s/u_x%dx%d_frame%04d_subs%02d_proc%04d.bmp"
         , get_out_path(lattice)
         , get_LX(lattice)
         , get_LY(lattice)
         , frame, subs, get_proc_id(lattice));
  if( !( o_ux = fopen( filename, "w+")))
  {
    printf("ERROR: fopen( \"%s\", \"w+\") = NULL.  Bye, bye!\n", filename);
    process_exit(1);
  }

  sprintf( filename
         , "%s/u_y%dx%d_frame%04d_subs%02d_proc%04d.bmp"
         , get_out_path(lattice)
         , get_LX(lattice)
         , get_LY(lattice)
         , frame, subs, get_proc_id(lattice));
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

        if( is_solid( lattice, n))
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

        if( is_solid( lattice, n))
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

void write_rho_image( lattice_ptr lattice, int subs)
{
  if( get_NumDims(lattice)==2)
  {
#if WRITE_MACRO_VAR_BMP_FILES
    rho2bmp( lattice, lattice->time);
#endif
  }
  else
  {
#if WRITE_MACRO_VAR_RAW_FILES
    char filename[1024];
    gen_filename( lattice, filename, "rho", get_frame(lattice), subs, ".raw");
    write_raw(
      lattice,
      get_rho_ptr(lattice,subs,0),
      /*stride*/ 1,
      /*max_u[0]*/get_max_rho(lattice,subs),
      /*min_u[0]*/get_min_rho(lattice,subs),
      filename);
#endif
  }
}
void write_u_image( lattice_ptr lattice, int subs)
{
  if( get_NumDims(lattice)==2)
  {
#if WRITE_MACRO_VAR_BMP_FILES
    u2bmp( lattice, lattice->time);
#endif
  }
  else
  {
#if WRITE_MACRO_VAR_RAW_FILES
    char filename[1024];
    gen_filename( lattice, filename, "u", get_frame(lattice), subs, ".raw");
    write_raw_u(
      lattice,
      get_ux_ptr(lattice,subs,0),
      get_uy_ptr(lattice,subs,0),
      get_uz_ptr(lattice,subs,0),
      /*stride*/ 1,
      /*max_u[0]*/get_max_ux(lattice,subs),
      /*TODO: Want to hardcode min=0 for some reason???*/
      /*min_u[0]*/get_min_ux(lattice,subs),
      /*max_u[1]*/get_max_uy(lattice,subs),
      /*min_u[1]*/get_min_uy(lattice,subs),
      /*max_u[2]*/get_max_uz(lattice,subs),
      /*min_u[2]*/get_min_uz(lattice,subs),
      filename);
#endif
  }
}

void write_ux_image( lattice_ptr lattice, int subs)
{
  if( get_NumDims(lattice)==3)
  {
#if WRITE_MACRO_VAR_RAW_FILES
    char filename[1024];
    gen_filename( lattice, filename, "u_x", get_frame(lattice), subs, ".raw");
    write_raw(
      lattice,
      get_ux_ptr(lattice,subs,0),
      /*stride*/ 1,
      /*max_u[0]*/get_max_ux(lattice,subs),
      /*TODO: Want to hardcode min=0 for some reason???*/
      /*min_u[0]*/get_min_ux(lattice,subs),
      filename);
#endif
  }
}
void write_uy_image( lattice_ptr lattice, int subs)
{
  if( get_NumDims(lattice)==3)
  {
#if WRITE_MACRO_VAR_RAW_FILES
    char filename[1024];
    gen_filename( lattice, filename, "u_y", get_frame(lattice), subs, ".raw");
    write_raw(
      lattice,
      get_uy_ptr(lattice,subs,0),
      /*stride*/ 1,
      /*max_u[0]*/get_max_uy(lattice,subs),
      /*TODO: Want to hardcode min=0 for some reason???*/
      /*min_u[0]*/get_min_uy(lattice,subs),
      filename);
#endif
  }
}
void write_uz_image( lattice_ptr lattice, int subs)
{
  if( get_NumDims(lattice)==3)
  {
#if WRITE_MACRO_VAR_RAW_FILES
    char filename[1024];
    gen_filename( lattice, filename, "u_z", get_frame(lattice), subs, ".raw");
    write_raw(
      lattice,
      get_uz_ptr(lattice,subs,0),
      /*stride*/ 1,
      /*max_u[0]*/get_max_uz(lattice,subs),
      /*TODO: Want to hardcode min=0 for some reason???*/
      /*min_u[0]*/get_min_uz(lattice,subs),
      filename);
#endif
  }
}

void write_rho_dat( lattice_ptr lattice, int subs)
{
#if WRITE_MACRO_VAR_DAT_FILES
  char filename[1024];
  gen_filename( lattice, filename, "rho", get_frame(lattice), subs, ".dat");
  write_dat(
    lattice
  , get_rho_ptr(lattice,subs,0)
  , /*stride*/ 1
  , /*max_u[1]*/0.
  , /*min_u[1]*/0.
  , filename );
#endif
}
void write_ux_dat( lattice_ptr lattice, int subs)
{
#if WRITE_MACRO_VAR_DAT_FILES
  char filename[1024];
  gen_filename( lattice, filename, "u_x", get_frame(lattice), subs, ".dat");
  write_dat(
    lattice,
    get_ux_ptr(lattice,subs,0),
    /*stride*/ 1,
    /*max_u[1]*/0.,
    /*min_u[1]*/0.,
    filename);
#endif
}
void write_uy_dat( lattice_ptr lattice, int subs)
{
#if WRITE_MACRO_VAR_DAT_FILES
  char filename[1024];
  gen_filename( lattice, filename, "u_y", get_frame(lattice), subs, ".dat");
  write_dat(
    lattice,
    get_uy_ptr(lattice,subs,0),
    /*stride*/ 1,
    /*max_u[1]*/0.,
    /*min_u[1]*/0.,
    filename);
#endif
}
void write_uz_dat( lattice_ptr lattice, int subs)
{
#if WRITE_MACRO_VAR_DAT_FILES
  char filename[1024];
  if( get_NumDims(lattice)==2) { return;}
  gen_filename( lattice, filename, "u_z", get_frame(lattice), subs, ".dat");
  write_dat(
    lattice,
    get_uz_ptr(lattice,subs,0),
    /*stride*/ 1,
    /*max_u[2]*/0.,
    /*min_u[2]*/0.,
    filename);
#endif
}

//##############################################################################
// DT 20111025
// {

real get_rho_A( lattice_ptr lattice, int subs)
{
  return lattice->param.rho_A[subs];
}

real get_rho_B( lattice_ptr lattice, int subs)
{
  return lattice->param.rho_B[subs];
}

void set_rho(
       lattice_ptr lattice
     , const int subs
     , const int n
     , const real rho)
{
  if( n>get_NumNodes(lattice))
  {
    printf("ERROR: n=%d > NumNodes=%d", n, get_NumNodes(lattice));
    process_exit(1);
  }
  lattice->vars[subs].macrovars1d[0][n] = rho;
}

void set_ux( lattice_ptr lattice, const int subs, const int n, const real ux)
{
  lattice->vars[subs].macrovars1d[1][n] = ux;
}

void set_uy( lattice_ptr lattice, const int subs, const int n, const real uy)
{
  lattice->vars[subs].macrovars1d[2][n] = uy;
}

void set_uz( lattice_ptr lattice, const int subs, const int n, const real uz)
{
  if( get_NumDims(lattice) == 3)
  {
    lattice->vars[subs].macrovars1d[3][n] = uz;
  }
}

void set_f(
       lattice_ptr lattice
     , const int subs
     , const int n
     , const int a
     , const real f)
{
  lattice->vars[subs].f1d[a][n] = f;
}

int do_diagnostic_init_of_macrovars( lattice_ptr lattice)
{
  return 0; // TODO: params.in or flags.in
}

real* get_max_ux_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.max_macrovars[subs] + 1;
}

real* get_max_uy_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.max_macrovars[subs] + 2;
}

real* get_max_uz_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.max_macrovars[subs] + 3;
}

real* get_min_ux_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.min_macrovars[subs] + 1;
}

real* get_min_uy_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.min_macrovars[subs] + 2;
}

real* get_min_uz_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.min_macrovars[subs] + 3;
}

real* get_ave_ux_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.ave_macrovars[subs] + 1;
}

real* get_ave_uy_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.ave_macrovars[subs] + 2;
}

real* get_ave_uz_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.ave_macrovars[subs] + 3;
}

real* get_max_rho_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.max_macrovars[subs] + 0;
}

real* get_min_rho_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.min_macrovars[subs] + 0;
}

real* get_ave_rho_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.ave_macrovars[subs] + 0;
}

real* get_flux_u_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.flux[subs] + 0;
}

real* get_flux_x_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.flux[subs] + 1;
}

real* get_flux_y_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.flux[subs] + 2;
}

real* get_flux_z_ptr( lattice_ptr lattice, const int subs)
{
  return lattice->attrib.flux[subs] + 3;
}

void set_time( lattice_ptr lattice, int time)
{
  lattice->time = time;
}

void set_frame( lattice_ptr lattice, int frame)
{
  lattice->frame = frame;
}

int get_NumTimeSteps( lattice_ptr lattice)
{
  return lattice->NumTimeSteps;
}

int is_end_of_frame( lattice_ptr lattice, int time)
{
  return !(time%get_FrameRate(lattice));
}

real* get_fptr( lattice_ptr lattice, int subs, int i, int j, int k, int a)
{
  int ni = get_ni(lattice);
  int nj = get_nj(lattice);
  int nk = get_nk(lattice);

  if( i<0) { i+=ni;}
  if( j<0) { j+=nj;}
  if( k<0) { k+=nk;}

  if( i>=ni) { i%=ni;}
  if( j>=nj) { j%=nj;}
  if( k>=nk) { k%=nk;}

  return lattice->vars[subs].f1d[a] + i + ni*j + ni*nj*k;
}

#ifdef __CUDACC__

__constant__ int vx_c[19];
__constant__ int vy_c[19];
__constant__ int vz_c[19];
__constant__ real wt_c[19];

__constant__ real tau_c[2];    // hardcoded for up to 2 fluid components
__constant__ real gaccel_c[6]; // hardcoded for up to 2 3D fluid components

__constant__ int numsubs_c;
__constant__ int numdims_c;
__constant__ int numdirs_c;
__constant__ int ni_c;
__constant__ int nj_c;
__constant__ int nk_c;
__constant__ int numnodes_c;


__device__ real get_f1d_d( real* f_mem_d, int subs, int i, int j, int k, int a)
{
  if( i<0) { i+=ni_c;}
  if( j<0) { j+=nj_c;}
  if( k<0) { k+=nk_c;}

  if( i>=ni_c) { i%=ni_c;}
  if( j>=nj_c) { j%=nj_c;}
  if( k>=nk_c) { k%=nk_c;}

  return f_mem_d[ subs * numnodes_c * numdirs_c
                + a * numnodes_c
                + i + ni_c*j + ni_c*nj_c*k];

//  return lattice_d->vars[subs].f1d[a][i + ni*j + ni*nj*k];
}
__device__ void set_mv_d( real* mv_mem_d, int subs,
                          int i, int j, int k, int a, real value)
{
  mv_mem_d[ subs*numnodes_c*(1 + numdims_c)
          + a*numnodes_c
          + i + ni_c*j + ni_c*nj_c*k] = value;

}

__device__ void set_f1d_d( real* f_mem_d, int subs,
                            int i, int j, int k, int a, real value)
{
  if( i<0) { i+=ni_c;}
  if( j<0) { j+=nj_c;}
  if( k<0) { k+=nk_c;}

  if( i>=ni_c) { i%=ni_c;}
  if( j>=nj_c) { j%=nj_c;}
  if( k>=nk_c) { k%=nk_c;}

  f_mem_d[ subs*numnodes_c*numdirs_c
         + a*numnodes_c
         + i + ni_c*j + ni_c*nj_c*k] = value;
}

//TODO? put __constant__ real vx[19], vy[19], vz[19] here?
//can't be dynamically allocated, but this may be OK - we have 64kB of const. mem.
//maybe calc and put a vectors of weights in const. mem. as well?


__device__ void calc_f_tilde_d(
                  real* f_mem_d
                , int subs
                , int dir
                , int thread
                , int block_size
                , real* f_temp
                , real usq)
{
  f_temp[thread + dir*block_size] *= (1. - 1. / tau_c[subs]);

  real vdotu = ((real) vx_c[dir])*f_temp[thread + (numdirs_c+1)*block_size]
             + ((real) vy_c[dir])*f_temp[thread + (numdirs_c+2)*block_size];
  if( numdims_c==3)
  {
    vdotu += ((real) vz_c[dir])*f_temp[thread + (numdirs_c+3)*block_size];
  }

  f_temp[thread + dir*block_size] += wt_c[dir]
                                   * f_temp[ thread + numdirs_c*block_size]
                                   * ( 1. + 3.*vdotu
                                          + 4.5*vdotu*vdotu
                                          + 1.5*usq
                                     ) / tau_c[subs];
}
//maybe later an equivalent function for the forcing term in the LBE, as per Guo.
__device__ void apply_accel_mv(
                  int subs
                , int cmpnt   //1, 2 or 3
                , int thread
                , int block_size
                , real* f_temp)
{
  f_temp[thread + (numdirs_c+cmpnt)*block_size]
    += gaccel_c[ subs*numdims_c + cmpnt-1];
}
#endif

#if 0
real* get_fptr( lattice_ptr lattice, int subs, int n, int a)
{
  if( n<0)
  {
    if( n==-1)
    {
      n+=get_ni(lattice);
    }
    else if( n<=-get_ni(lattice) && n>=-2*get_ni(lattice)+1)
    {
      n+=get_ni(lattice)*(get_nj(lattice)+0);
    }
    else
    {
      n+=get_NumNodes(lattice);
    }
  }
  if( n>=get_NumNodes(lattice))
  {
    if( n==get_NumNodes(lattice))
    {
      n-=get_ni(lattice);
    }
    else if( n>=get_NumNodes(lattice)+get_ni(lattice)-1
          && n<=get_NumNodes(lattice)+2*(get_ni(lattice)-1))
    {
      n-=get_ni(lattice)*get_nj(lattice);
    }
    else
    {
      n-=get_NumNodes(lattice);
    }
  }
  if( n<0 || n>=get_NumNodes(lattice))
  {
    printf("%s %d ERROR: linear index n=%d out of bounds.\n"
          ,__FILE__,__LINE__,n);
    process_exit(lattice);
  }
  return lattice->vars[subs].f1d[a] + n;
}
#endif

real get_tau( lattice_ptr lattice, int subs)
{
  return lattice->param.tau[subs];
}
real* get_tau_ptr( lattice_ptr lattice)
{
  return lattice->param.tau;
}

real get_gaccel( lattice_ptr lattice, int subs, int cmpt)
{
  return lattice->param.gforce[subs][cmpt];
}
real* get_gaccel_ptr( lattice_ptr lattice, int subs)
{
  return lattice->param.gforce[subs];
}

int write_debug_txt_files( lattice_ptr lattice)
{
  return 1; // TODO: params.in or flags.in
}

int get_time( lattice_ptr lattice)
{
  return lattice->time;
}

void set_nk( lattice_ptr lattice, int nk)
{
  lattice->param.LZ = nk;
}

// }
//##############################################################################
