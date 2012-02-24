//##############################################################################
//
// lattice.c
//

int get_LX( lattice_ptr lattice) { return lattice->param.LX;}
int get_LY( lattice_ptr lattice) { return lattice->param.LY;}
int get_LZ( lattice_ptr lattice) { return lattice->param.LZ;}
int get_BX( lattice_ptr lattice) { return lattice->param.BX;}
int get_BY( lattice_ptr lattice) { return lattice->param.BY;}
int get_BZ( lattice_ptr lattice) { return lattice->param.BZ;}
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

int get_NumFrames( lattice_ptr lattice)
{
  return lattice->param.NumFrames;
}

int get_frame( lattice_ptr lattice)
{
  return lattice->frame;
}

void set_NumNodes( lattice_ptr lattice)
{
  lattice->NumNodes = get_LX( lattice)*get_LY( lattice)*get_LZ( lattice);
}

int get_NumNodes( lattice_ptr lattice)
{ 
  return lattice->NumNodes;
}

int* get_NumNodes_ptr( lattice_ptr lattice) 
{ 
  return &(lattice->NumNodes);
}


void set_NumDims( lattice_ptr lattice, int value)
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

void set_NumVelDirs( lattice_ptr lattice, int value)
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

void set_NumBoundDirs( lattice_ptr lattice, int value)
{
  lattice->NumBoundDirs = value;
}

int get_NumBoundDirs( lattice_ptr lattice)
{
  return lattice->NumBoundDirs;  
  //#if NUM_DIMENSIONS == 2
  //  return 6;
  //#else
  //  return 10;
  //#endif
}

int* get_NumBoundDirs_ptr( lattice_ptr lattice)
{
  return &(lattice->NumBoundDirs);  
  //#if NUM_DIMENSIONS == 2
  //  return 6;
  //#else
  //  return 10;
  //#endif
}

void set_NumUnboundDirs( lattice_ptr lattice, int value)
{
  lattice->NumUnboundDirs = value;
}

int get_NumUnboundDirs( lattice_ptr lattice)
{
  return lattice->NumUnboundDirs;
}

int* get_NumUnboundDirs_ptr( lattice_ptr lattice)
{
  return &(lattice->NumUnboundDirs);
}

void set_EndBoundSize( lattice_ptr lattice, int value)
{
  lattice->EndBoundSize = value;
}

int get_EndBoundSize( lattice_ptr lattice)
{
  return lattice->EndBoundSize;
}

int* get_EndBoundSize_ptr( lattice_ptr lattice)
{
  return &(lattice->EndBoundSize);
}

void set_NumSubs( lattice_ptr lattice, int value)
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

void set_ns(
    lattice_ptr lattice
    , const int n
    , const real val)
{
  lattice->ns_memblock[n] = val;
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

void display_etime( lattice_ptr lattice)
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
        ,"%s/%s%s%s%s%s_proc%04d%s"
        , get_out_path( lattice)
#ifdef __CUDACC__
        , "cuda_"
#else
        , ""
#endif
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
        ,"%s/%s%s%s%s%s%s"
        , get_out_path( lattice)
#ifdef __CUDACC__
        , "cuda_"
#else
        , ""
#endif
        , prefix
        , size_str
        , frame_str
        , subs_str
        , suffix
        );
  }
}

#if 1 // __CUDACC__
// The definitions of rho2bmp and u2bmp need to be relocated here from lbio.c
// until the way this project is compiled is fixed to allow separate
// compilation and linking.

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
  FILE   *o;
  int    i, j,
         n;
  int    pad,
         bytes_per_row;
  struct bitmap_file_header bmfh;
  struct bitmap_info_header bmih;
  int    *width_ptr;
  short  int *bitcount_ptr;
  char   filename[1024];
  char   red_val,
         green_val,
         blue_val,
         val;
  real fval;
  real min_rho, max_rho;
  int    subs;

#if SAY_HI
  printf("rho2bmp() -- Hi!\n");
#endif /* SAY_HI */

  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
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
  FILE   *o_u,
         *o_ux,
         *o_uy;
  int    i, j,
         n;
  int    pad,
         bytes_per_row;
  int    frame;
  struct bitmap_file_header bmfh;
  struct bitmap_info_header bmih;
  int    *width_ptr;
  short  int *bitcount_ptr;
  char   filename[1024];
  char   red_val,
         green_val,
         blue_val,
         val;
  real max_u[2]; //, maxu, u;
  real u_x, u_y;
  int    subs;

#if SAY_HI
  printf("u2bmp() -- Hi!\n");
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
    bitcount_ptr = (short int*)bmih.biBitCount;

    // Bytes per row of the bitmap.
    bytes_per_row =
      ((int)ceil(( (((real)(ENDIAN4(*width_ptr)))*((real)(ENDIAN2(*bitcount_ptr))))/8.)));

    // Bitmaps pad rows to preserve 4-byte boundaries.
    // The length of a row in the file will be bytes_per_row + pad .
    pad = ((4) - bytes_per_row%4)%4;

    max_u[0] = get_max_ux(lattice,subs);
    max_u[1] = get_max_uy(lattice,subs);

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

          //u = sqrt(u_x*u_x + u_y*u_y);
          //maxu = sqrt( max_u[0]*max_u[0] + max_u[1]*max_u[1]);

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
  FILE *in;
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
// W R I T E _ E M P T Y _ B M P  {{{
//##############################################################################
// void write_empty_bmp( char *filename, int time)
//
void write_empty_bmp( lattice_ptr lattice)
{
  FILE   *o;
  int    i, j;
  int    pad,
         bytes_per_row;
  struct bitmap_file_header bmfh;
  struct bitmap_info_header bmih;
  int    *width_ptr;
  short  int *bitcount_ptr;
  char   filename[1024];
  char   red_val,
         green_val,
         blue_val,
         val;

#if SAY_HI
  printf("write_empty_bmp() -- Hi!\n");
#endif /* SAY_HI */

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
  bitcount_ptr = (short int*)bmih.biBitCount;

  // Bytes per row of the bitmap.
  bytes_per_row =
    ((int)ceil(( ( ((real)(ENDIAN4(*width_ptr)))
          * ((real)(ENDIAN2(*bitcount_ptr))) )/8.)));

  // Bitmaps pad rows to preserve 4-byte boundaries.
  // The length of a row in the file will be bytes_per_row + pad .
  pad = ((4) - bytes_per_row%4)%4;

  sprintf(filename,"./in/%dx%d.bmp",get_ni(lattice),get_nj(lattice));

  if( !( o = fopen( filename, "w+")))
  {
    printf("ERROR: fopen( \"%s\", \"w+\") = NULL.  Bye, bye!\n", filename);
    process_exit(1);
  }

  fwrite( &bmfh, sizeof(struct bitmap_file_header), 1, o );
  fwrite( &bmih, sizeof(struct bitmap_info_header), 1, o );

  for( j=0; j<get_LY(lattice); j++)
  {
    for( i=0; i<get_LX(lattice); i++)
    {
      printf("writing pixel (%d,%d)\n",i,j);
      red_val   = (char)255;
      green_val = (char)255;
      blue_val  = (char)255;
      val       = (char)255;

      if( fwrite( &blue_val, 1, 1, o) != 1) { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &green_val, 1, 1, o) != 1) { printf("BOOM!\n"); process_exit(1);}
      if( fwrite( &red_val, 1, 1, o) != 1) { printf("BOOM!\n"); process_exit(1);}

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
  printf("write_empty_bmp() -- Wrote file \"%s\".\n", filename);
#endif /* VERBOSITY_LEVEL > 0 */

#if SAY_HI
  printf("write_empty_bmp() -- Bye!\n");
  printf("\n");
#endif /* SAY_HI */

} /* write_empty_bmp( lattice_ptr lattice, int time) */
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
#if INAMURO_SIGMA_COMPONENT
real get_rho_sigma( lattice_ptr lattice)
{
  return lattice->param.rho_sigma;
}
#endif
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

real get_f(
    lattice_ptr lattice
    , const int subs
    , const int n
    , const int a )
{
  return lattice->vars[subs].f1d[a][n];
}

int do_diagnostic_init_of_rho( lattice_ptr lattice)
{
  return 0; // TODO: params.in or flags.in
}

int do_diagnostic_init_of_u( lattice_ptr lattice)
{
  return 0; // TODO: params.in or flags.in
}

int do_diagnostic_init_of_f( lattice_ptr lattice)
{
  return 0; // TODO: params.in or flags.in
}

int skip_collision_step( lattice_ptr lattice)
{
  return 0; // TODO: params.in or flags.in
}

int skip_body_force_term( lattice_ptr lattice)
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

real get_gaccel_ux( lattice_ptr lattice, int subs)
{
  return lattice->param.gforce[subs][0];
}
real get_gaccel_uy( lattice_ptr lattice, int subs)
{
  return lattice->param.gforce[subs][1];
}
real get_gaccel_uz( lattice_ptr lattice, int subs)
{
  if( get_NumDims(lattice)==3)
  {
    return lattice->param.gforce[subs][2];
  }
  else
  {
    return 0.;
  }
}

#if 0
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
#else
real* get_fptr(
    lattice_ptr lattice
    , int subs
    , int i0
    , int j0
    , int k0
    , int di
    , int dj
    , int dk
    , int a )
{
  int ni = get_ni(lattice);
  int nj = get_nj(lattice);
  int nk = get_nk(lattice);

  if( di!=0 || dj!=0 || dk!=0)
  {
    // Getting f from neighboring node (i+di,j+dj,k+dk). This is for the
    // stream_collide_stream step.
    int i = i0+di;
    int j = j0+dj;
    int k = k0+dk;

    if( i<0) { i+=ni;}
    if( j<0) { j+=nj;}
    if( k<0) { k+=nk;}

    if( i>=ni) { i%=ni;}
    if( j>=nj) { j%=nj;}
    if( k>=nk) { k%=nk;}

    int n = i + ni*j + ni*nj*k;

    if( is_not_solid( lattice, n))
    {
      return lattice->vars[subs].f1d[a] + n;
    }
    else
    {
      // If neighboring node is a solid, return the f at node (i0,j0,k0) that
      // would be streamed out for halfway bounceback.
      return lattice->vars[subs].f1d[ a + ((!(a%2))?(-1):(1)) ]
        + i0 + ni*j0 + ni*nj*k0;
    }
  }
  else
  {
    // Getting f from node (i0,j0,k0). This is for the even collide step that
    // wraps up what the stream_collide_stream step started.
    if( i0<0) { i0+=ni;}
    if( j0<0) { j0+=nj;}
    if( k0<0) { k0+=nk;}

    if( i0>=ni) { i0%=ni;}
    if( j0>=nj) { j0%=nj;}
    if( k0>=nk) { k0%=nk;}

    return lattice->vars[subs].f1d[a] + i0 + ni*j0 + ni*nj*k0;
  }
}
#endif

unsigned char* get_solids_ptr( lattice_ptr lattice, const int n)
{
  return lattice->solids_memblock + n;
}

real* get_ns_ptr( lattice_ptr lattice, const int n)
{
  return lattice->ns_memblock + n;
}


#ifdef __CUDACC__

__constant__ int vx_c[19];     //
__constant__ int vy_c[19];     // Enough space for D3Q19; first 9
__constant__ int vz_c[19];     // components constitute D2Q9 model
__constant__ real wt_c[19];    // 
__constant__ int cumul_stride_c[20]; // For variable stride created by boundary regions

__constant__ int is_end_of_frame_mem_c;

__constant__ real tau_c[2];    // Hardcoded for up to 2 fluid components
__constant__ real gaccel_c[6]; // Hardcoded for up to 2 3D fluid components

__constant__ int numsubs_c;
__constant__ int numdims_c;
__constant__ int numdirs_c;
__constant__ int numbounddirs_c;
__constant__ int numunbounddirs_c;
__constant__ int end_bound_c;
__constant__ int proc_id_c;
__constant__ int num_procs_c;
__constant__ int ni_c;
__constant__ int nj_c;
__constant__ int nk_c;
__constant__ int kloop_c;
__constant__ int bixbj_c;
__constant__ int bjxbk_c;
__constant__ int bixbk_c;
__constant__ int nixnj_c;
__constant__ int blocksize_c;
__constant__ int numnodes_c;
__constant__ int subs_c;

__constant__ real fixed_bound_var_c;
#if SIMPLE_NS_CUSTOM
__constant__ real ns_c;
#endif

#if TEXTURE_FETCH
texture<unsigned char, 1, cudaReadModeElementType> tex_solid;
#endif


__device__ int d_is_solid( unsigned char* solids_mem_d, int n)
{
#if TEXTURE_FETCH
  if( tex1Dfetch(tex_solid, n) == 1)
  {
    return 1.;
  }
  else
  {
    return 0.;
  }
#else
  if( solids_mem_d[n] == 1)
  {
    return 1.;
  }
  else
  {
    return 0.;
  }
#endif
}

__device__ int d_is_not_solid( unsigned char* solids_mem_d, int n)
{ 
#if TEXTURE_FETCH
  if( tex1Dfetch(tex_solid, n) == 1)
  {
    return 0.;
  }
  else
  {
    return 1.;
  }
#else
  if( solids_mem_d[n] == 1)
  {
    return 0.;
  }
  else
  {
    return 1.;
  }
#endif
}

void pdf_boundary_swap( lattice_ptr lattice
    , real* f_mem_d
    , int* cstr
    , int subs
    , int dir
    , int time 
    )
{
  int in_boundary, in_main;
  if (dir&1) //dir%2 -- if dir is odd
  {
    in_boundary = subs * cstr[get_NumVelDirs( lattice)]
      + cstr[dir] - get_EndBoundSize( lattice);
    in_main = in_boundary + get_NumNodes( lattice);
  }
  else
  {
    in_main = subs * cstr[get_NumVelDirs( lattice)]
      + cstr[dir];
    in_boundary = in_boundary + get_NumNodes( lattice);
  }


  if( time&1) //time%2 -- if time is odd
  {
    cudaMemcpy( f_mem_d + in_boundary
        , f_mem_d + in_main
        , get_EndBoundSize( lattice)
        *sizeof(real)
        , cudaMemcpyDeviceToDevice);

    checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
  }
  else
  {
    cudaMemcpy( f_mem_d + in_main
        , f_mem_d + in_boundary
        , get_EndBoundSize( lattice)
        *sizeof(real)
        , cudaMemcpyDeviceToDevice);

    checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
  }
}

#if PARALLEL
void pdf_bound_DtH_odd_1( lattice_ptr lattice
    , real* f_mem_d
    , int* cstr
    , int subs
    , int dir
    )
{
  cudaMemcpy( lattice->process.pos_dir_pdf_to_send
      + subs * get_NumBoundDirs( lattice) * get_EndBoundSize( lattice) / 2
      + (dir - get_NumUnboundDirs( lattice)) * get_EndBoundSize( lattice) / 2
      , f_mem_d + subs * cstr[get_NumVelDirs( lattice)]
      + cstr[dir] - get_EndBoundSize( lattice)+ get_NumNodes( lattice)
      , get_EndBoundSize( lattice)
      *sizeof(real)
      , cudaMemcpyDeviceToHost);

  checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
}

void pdf_bound_DtH_even_1( lattice_ptr lattice
    , real* f_mem_d
    , int* cstr
    , int subs
    , int dir
    )
{
  cudaMemcpy( lattice->process.neg_dir_pdf_to_send
      + subs * get_NumBoundDirs( lattice) * get_EndBoundSize( lattice) / 2
      + (dir - get_NumUnboundDirs( lattice) - 1) * get_EndBoundSize( lattice) / 2
      , f_mem_d + subs * cstr[get_NumVelDirs( lattice)]
      + cstr[dir]
      , get_EndBoundSize( lattice)
      *sizeof(real)
      , cudaMemcpyDeviceToHost);

  checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
}

void pdf_bound_DtH_odd_2( lattice_ptr lattice
    , real* f_mem_d
    , int* cstr
    , int subs
    , int dir
    )
{
  cudaMemcpy( lattice->process.neg_dir_pdf_to_send
      + subs * get_NumBoundDirs( lattice) * get_EndBoundSize( lattice) / 2
      + (dir - get_NumUnboundDirs( lattice)) * get_EndBoundSize( lattice) / 2
      , f_mem_d + subs * cstr[get_NumVelDirs( lattice)]
      + cstr[dir] - get_EndBoundSize( lattice)
      , get_EndBoundSize( lattice)
      *sizeof(real)
      , cudaMemcpyDeviceToHost);

  checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
}

void pdf_bound_DtH_even_2( lattice_ptr lattice
    , real* f_mem_d
    , int* cstr
    , int subs
    , int dir
    )
{
  cudaMemcpy( lattice->process.pos_dir_pdf_to_send
      + subs * get_NumBoundDirs( lattice) * get_EndBoundSize( lattice) / 2
      + (dir - get_NumUnboundDirs( lattice) - 1) * get_EndBoundSize( lattice) / 2
      , f_mem_d + subs * cstr[get_NumVelDirs( lattice)]
      + cstr[dir] + get_NumNodes( lattice)
      , get_EndBoundSize( lattice)
      *sizeof(real)
      , cudaMemcpyDeviceToHost);

  checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
}

void pdf_bound_HtD_odd_1( lattice_ptr lattice
    , real* f_mem_d
    , int* cstr
    , int subs
    , int dir
    )
{
  cudaMemcpy( f_mem_d + subs * cstr[get_NumVelDirs( lattice)]
      + cstr[dir] - get_EndBoundSize( lattice)
      , lattice->process.pos_dir_pdf_to_recv
      + subs * get_NumBoundDirs( lattice) * get_EndBoundSize( lattice) / 2
      + (dir - get_NumUnboundDirs( lattice)) * get_EndBoundSize( lattice) / 2
      , get_EndBoundSize( lattice)
      *sizeof(real)
      , cudaMemcpyHostToDevice);

  checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
}
void pdf_bound_HtD_even_1( lattice_ptr lattice
    , real* f_mem_d
    , int* cstr
    , int subs
    , int dir
    )
{
  cudaMemcpy( f_mem_d + subs * cstr[get_NumVelDirs( lattice)]
      + cstr[dir] + get_NumNodes( lattice)
      , lattice->process.neg_dir_pdf_to_recv
      + subs * get_NumBoundDirs( lattice) * get_EndBoundSize( lattice) / 2
      + (dir - get_NumUnboundDirs( lattice) - 1) * get_EndBoundSize( lattice) / 2
      , get_EndBoundSize( lattice)
      *sizeof(real)
      , cudaMemcpyHostToDevice);

  checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
}

void pdf_bound_HtD_odd_2( lattice_ptr lattice
    , real* f_mem_d
    , int* cstr
    , int subs
    , int dir
    )
{
  cudaMemcpy( f_mem_d + subs * cstr[get_NumVelDirs( lattice)]
      + cstr[dir] - get_EndBoundSize( lattice) + get_NumNodes( lattice)
      , lattice->process.neg_dir_pdf_to_recv
      + subs * get_NumBoundDirs( lattice) * get_EndBoundSize( lattice) / 2
      + (dir - get_NumUnboundDirs( lattice)) * get_EndBoundSize( lattice) / 2
      , get_EndBoundSize( lattice)
      *sizeof(real)
      , cudaMemcpyHostToDevice);

  checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
}

void pdf_bound_HtD_even_2( lattice_ptr lattice
    , real* f_mem_d
    , int* cstr
    , int subs
    , int dir
    )
{
  cudaMemcpy( f_mem_d + subs * cstr[get_NumVelDirs( lattice)]
      + cstr[dir]
      , lattice->process.pos_dir_pdf_to_recv
      + subs * get_NumBoundDirs( lattice) * get_EndBoundSize( lattice) / 2
      + (dir - get_NumUnboundDirs( lattice) - 1) * get_EndBoundSize( lattice) / 2
      , get_EndBoundSize( lattice)
      *sizeof(real)
      , cudaMemcpyHostToDevice);

  checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
}


void pdf_boundary_parallel( lattice_ptr lattice
    , real* f_mem_d
    , int* cstr
    , int subs
    , int dir
    , int time 
    , int going_to_host
    )
{
  int in_boundary, in_main;

#if 1
  if ( dir&1) //dir%2 -- if dir is odd
  {
    in_boundary = subs * cstr[get_NumVelDirs( lattice)]
      + cstr[dir] - get_EndBoundSize( lattice);
    in_main = in_boundary + get_NumNodes( lattice);

    if( going_to_host)
    {
      if( time&1) //time%2 -- if time is odd
      {
        cudaMemcpy( lattice->process.pos_dir_pdf_to_send
            + (dir - get_NumUnboundDirs( lattice)) * get_EndBoundSize( lattice) / 2
            , f_mem_d + in_main
            , get_EndBoundSize( lattice)
            *sizeof(real)
            , cudaMemcpyDeviceToHost);

        checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
      }
      else
      {
        cudaMemcpy( lattice->process.neg_dir_pdf_to_send
            + (dir - get_NumUnboundDirs( lattice)) * get_EndBoundSize( lattice) / 2
            , f_mem_d + in_boundary
            , get_EndBoundSize( lattice)
            *sizeof(real)
            , cudaMemcpyDeviceToHost);

        checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
      }

    }
    else  // !(going_to_host)
    {
      if( time&1) //time%2 -- if time is odd
      {
        cudaMemcpy( f_mem_d + in_boundary
            , lattice->process.pos_dir_pdf_to_recv
            + (dir - get_NumUnboundDirs( lattice)) * get_EndBoundSize( lattice) / 2
            , get_EndBoundSize( lattice)
            *sizeof(real)
            , cudaMemcpyHostToDevice);

        checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
      }
      else
      {
        cudaMemcpy( f_mem_d + in_main
            , lattice->process.neg_dir_pdf_to_recv
            + (dir - get_NumUnboundDirs( lattice)) * get_EndBoundSize( lattice) / 2
            , get_EndBoundSize( lattice)
            *sizeof(real)
            , cudaMemcpyHostToDevice);

        checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
      }

    }

  }
  else
  {
    in_main = subs * cstr[get_NumVelDirs( lattice)]
      + cstr[dir];
    in_boundary = in_main + get_NumNodes( lattice);

    if( going_to_host)
    {
      if( time&1) //time%2 -- if time is odd
      {
        cudaMemcpy( lattice->process.neg_dir_pdf_to_send
            + (dir - get_NumUnboundDirs( lattice) - 1) * get_EndBoundSize( lattice) / 2
            , f_mem_d + in_main
            , get_EndBoundSize( lattice)
            *sizeof(real)
            , cudaMemcpyDeviceToHost);

        checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
      }
      else
      {
        cudaMemcpy( lattice->process.pos_dir_pdf_to_send
            + (dir - get_NumUnboundDirs( lattice) - 1) * get_EndBoundSize( lattice) / 2
            , f_mem_d + in_boundary
            , get_EndBoundSize( lattice)
            *sizeof(real)
            , cudaMemcpyDeviceToHost);

        checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
      }

    }
    else  // !(going_to_host)
    {
      if( time&1) //time%2 -- if time is odd
      {
        cudaMemcpy( f_mem_d + in_boundary
            , lattice->process.neg_dir_pdf_to_recv
            + (dir - get_NumUnboundDirs( lattice) - 1) * get_EndBoundSize( lattice) / 2
            , get_EndBoundSize( lattice)
            *sizeof(real)
            , cudaMemcpyHostToDevice);

        checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
      }
      else
      {
        cudaMemcpy( f_mem_d + in_main
            , lattice->process.pos_dir_pdf_to_recv
            + (dir - get_NumUnboundDirs( lattice) - 1) * get_EndBoundSize( lattice) / 2
            , get_EndBoundSize( lattice)
            *sizeof(real)
            , cudaMemcpyHostToDevice);

        checkCUDAError( __FILE__, __LINE__, "boundary swap"); 
      }

    }

  }

#endif
}
#endif

__device__ real get_f1d_d(
    real* f_mem_d
    , unsigned char* solids_mem_d
    , real* ns_mem_d
    , int subs
    , int i0
    , int j0
    , int k0
    , int n0
    , int di
    , int dj
    , int dk
    , int a
    , int da
#if SIMPLE_NS_CUSTOM
    , int alpha_switch 
#endif
    )
{
  // Getting f_a from node (i+di,j+dj,k+dk).

#if SIMPLE_NS_CUSTOM
  int i = i0+di;
  int j = j0+dj;
  int k = k0+dk;

#if PARALLEL
  if( i<0) { i+=ni_c;}
  if( i==ni_c) { i=0;}

  if( numdims_c == 3)
  {
    if( j<0) { j+=nj_c;}
    if( j==nj_c) { j=0;}
  }
#else
  if( i<0) { i+=ni_c;}
  if( i==ni_c) { i=0;}

  if( j<0) { j+=nj_c;}
  if( j==nj_c) { j=0;}

  if( k<0) { k+=nk_c;}
  if( k==nk_c) { k=0;}
#endif

  int n = i + ni_c*j + nixnj_c*k;

  real alpha = (1. - ns_mem_d[n + end_bound_c]) * (1. - ns_mem_d[n0 + end_bound_c]);
  //real alpha = (real) ((1 - solids_mem_d[n + end_bound_c]) * (1 - solids_mem_d[n0 + end_bound_c]));

  if( alpha_switch)     // k_stream_collide_stream
  {
    return alpha 
      * f_mem_d[ subs*cumul_stride_c[numdirs_c] 
      + cumul_stride_c[a] + n]
      + (1. - alpha)
      * f_mem_d[ subs*cumul_stride_c[numdirs_c] 
      + cumul_stride_c[a+da] + n0];
  }
  else  // k_collide
  {
    return (1. - alpha) 
      * f_mem_d[ subs*cumul_stride_c[numdirs_c] 
      + cumul_stride_c[a] + n]
      + alpha
      * f_mem_d[ subs*cumul_stride_c[numdirs_c] 
      + cumul_stride_c[a+da] + n0];
  }
#else   // !(SIMPLE_NS_CUSTOM)
#if COMPUTE_ON_SOLIDS

  if( d_is_not_solid( solids_mem_d, n0 + end_bound_c))
  {
    int i = i0+di;
    int j = j0+dj;
    int k = k0+dk;

#if PARALLEL
    if( i<0) { i+=ni_c;}
    if( i==ni_c) { i=0;}

    if( numdims_c == 3)
    {
      if( j<0) { j+=nj_c;}
      if( j==nj_c) { j=0;}
    }
#else
    if( i<0) { i+=ni_c;}
    if( i==ni_c) { i=0;}

    if( j<0) { j+=nj_c;}
    if( j==nj_c) { j=0;}

    if( k<0) { k+=nk_c;}
    if( k==nk_c) { k=0;}
#endif

    int n = i + ni_c*j + nixnj_c*k;

    if( d_is_not_solid( solids_mem_d, n + end_bound_c))
    {
      return f_mem_d[ subs*cumul_stride_c[numdirs_c] 
        + cumul_stride_c[a] + n];
    }
    else
    { 
      // If neighboring node is a solid, return the f at node (i0,j0,k0) that
      // would be streamed out for halfway bounceback.
      return f_mem_d[ subs*cumul_stride_c[numdirs_c]
        + cumul_stride_c[a+da] + n0];
    }
  }
  else
  {
    return 0.;
  }

#else   // !(COMPUTE_ON_SOLIDS)
  int i = i0+di;
  int j = j0+dj;
  int k = k0+dk;

#if PARALLEL
  if( i<0) { i+=ni_c;}
  if( i==ni_c) { i=0;}

  if( numdims_c == 3)
  {
    if( j<0) { j+=nj_c;}
    if( j==nj_c) { j=0;}
  }
#else
  if( i<0) { i+=ni_c;}
  if( i==ni_c) { i=0;}

  if( j<0) { j+=nj_c;}
  if( j==nj_c) { j=0;}

  if( k<0) { k+=nk_c;}
  if( k==nk_c) { k=0;}
#endif

  int n = i + ni_c*j + nixnj_c*k;

  if( d_is_not_solid( solids_mem_d, n + end_bound_c))
  {
    return f_mem_d[ subs*cumul_stride_c[numdirs_c] 
      + cumul_stride_c[a] + n];
  }
  else
  { 
    // If neighboring node is a solid, return the f at node (i0,j0,k0) that
    // would be streamed out for halfway bounceback.
    return f_mem_d[ subs*cumul_stride_c[numdirs_c]
      + cumul_stride_c[a+da] + n0];
  }
#endif  // COMPUTE_ON_SOLIDS
#endif  // SIMPLE_NS_CUSTOM

}



__device__ void set_mv_d( real* mv_mem_d, int subs,
    int n, int a, real value)
{
  mv_mem_d[ subs*numnodes_c*(1 + numdims_c)
    + a*numnodes_c + n ] = value;

}

__device__ real get_mv_d( real* mv_mem_d, int subs,
    int n, int a)
{
  return mv_mem_d[ subs*numnodes_c*(1 + numdims_c)
    + a*numnodes_c + n ];

}


__device__ void set_f1d_d(
    real* f_mem_d
    , unsigned char* solids_mem_d
    , real* ns_mem_d
    , int subs
    , int i0
    , int j0
    , int k0
    , int n0
    , int di
    , int dj
    , int dk
    , int a 
    , real value)
{
  // Setting f to node (i+di,j+dj,k+dk). The 'd_is_not_solid' conditional
  // statement is vital for the correct functioning of the bounceback
  // boundary conditions. If !(COMPUTE_ON_SOLIDS), this conditional statement
  // exists inside the kernels instead.

#if SIMPLE_NS_CUSTOM
  int i = i0+di;
  int j = j0+dj;
  int k = k0+dk;

#if PARALLEL
  if( i<0) { i+=ni_c;}
  if( i==ni_c) { i=0;}

  if( numdims_c == 3)
  {
    if( j<0) { j+=nj_c;}
    if( j==nj_c) { j=0;}
  }
#else
  if( i<0) { i+=ni_c;}
  if( i==ni_c) { i=0;}

  if( j<0) { j+=nj_c;}
  if( j==nj_c) { j=0;}

  if( k<0) { k+=nk_c;}
  if( k==nk_c) { k=0;}
#endif

  int n = i + ni_c*j + nixnj_c*k;

  f_mem_d[ subs*cumul_stride_c[numdirs_c] 
    + cumul_stride_c[a] + n] = value;

#else   // !(SIMPLE_NS_CUSTOM)
#if COMPUTE_ON_SOLIDS
  if( d_is_not_solid( solids_mem_d, n0 + end_bound_c))
  {
    int i = i0+di;
    int j = j0+dj;
    int k = k0+dk;

#if PARALLEL
    if( i<0) { i+=ni_c;}
    if( i==ni_c) { i=0;}

    if( numdims_c == 3)
    {
      if( j<0) { j+=nj_c;}
      if( j==nj_c) { j=0;}
    }
#else
    if( i<0) { i+=ni_c;}
    if( i==ni_c) { i=0;}

    if( j<0) { j+=nj_c;}
    if( j==nj_c) { j=0;}

    if( k<0) { k+=nk_c;}
    if( k==nk_c) { k=0;}
#endif

    int n = i + ni_c*j + nixnj_c*k;

    if( d_is_not_solid( solids_mem_d, n + end_bound_c))
    {
      f_mem_d[ subs*cumul_stride_c[numdirs_c] 
        + cumul_stride_c[a] + n] = value;
    }
    else
    {
      f_mem_d[ subs*cumul_stride_c[numdirs_c] 
        + cumul_stride_c[a] + n] = 0.;
    }
  }

#else   // !(COMPUTE_ON_SOLIDS)
  int i = i0+di;
  int j = j0+dj;
  int k = k0+dk;

#if PARALLEL
  if( i<0) { i+=ni_c;}
  if( i==ni_c) { i=0;}

  if( numdims_c == 3)
  {
    if( j<0) { j+=nj_c;}
    if( j==nj_c) { j=0;}
  }
#else
  if( i<0) { i+=ni_c;}
  if( i==ni_c) { i=0;}

  if( j<0) { j+=nj_c;}
  if( j==nj_c) { j=0;}

  if( k<0) { k+=nk_c;}
  if( k==nk_c) { k=0;}
#endif

  int n = i + ni_c*j + nixnj_c*k;

  if( d_is_not_solid( solids_mem_d, n + end_bound_c))
  {
    f_mem_d[ subs*cumul_stride_c[numdirs_c] 
      + cumul_stride_c[a] + n] = value;
  }

#endif  // COMPUTE_ON_SOLIDS
#endif  // SIMPLE_NS_CUSTOM


}
__device__ void calc_f_tilde_d(
    real* f_mem_d
    , int subs
    , int dir
    , int thread
    , int block_size
    , real* f_temp
    , real usq)
{
#if INAMURO_SIGMA_COMPONENT
  if( subs == 1)
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
      * ( 1. + 3.*vdotu) / tau_c[subs];

    return;
  }
#endif

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
        - 1.5*usq
      ) / tau_c[subs];
}
// Maybe later an equivalent function for the forcing term in the LBE, as per
// Guo.

__device__ void apply_accel_mv(
    int subs
    , int cmpnt   //1, 2 or 3
    , int thread
    , int block_size
    , real* f_temp)
{
#if 1
  f_temp[thread + (numdirs_c+cmpnt)*block_size]
    += gaccel_c[ subs*3 + cmpnt-1];
#else
  // DT: Testing
  if( subs==0 && cmpnt==1)
  {
    f_temp[thread + (numdirs_c+cmpnt)*block_size] += 0.00001;
  }
#endif
}

  __device__
int d_skip_collision_step()
{
  return 0; // TODO: params.in or flags.in
}

  __device__
int d_skip_body_force_term()
{
  return 0; // TODO: params.in or flags.in
}

  __device__
int d_skip_updating_macrovars()
{
  return 1; // TODO: params.in or flags.in
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

int do_write_rho_txt_file( lattice_ptr lattice)
{
  return 1; // TODO: params.in or flags.in
}

int do_write_u_txt_file( lattice_ptr lattice)
{
  return 1; // TODO: params.in or flags.in
}

int do_write_f_txt_file( lattice_ptr lattice)
{
  return 0; // TODO: params.in or flags.in
}

int get_time( lattice_ptr lattice)
{
  return lattice->time;
}

void set_nk( lattice_ptr lattice, int nk)
{
  lattice->param.LZ = nk;
}

real* pressure_n_in0( lattice_ptr lattice, int subs)
{
  return lattice->bcs_in[subs].pressure_n_in0;
}

real** pressure_n_in0_ptr( lattice_ptr lattice, int subs)
{
  return &( lattice->bcs_in[subs].pressure_n_in0);
}

int num_pressure_n_in0( lattice_ptr lattice, int subs)
{
  return lattice->bcs_in[subs].num_pressure_n_in0;
}

int* num_pressure_n_in0_ptr( lattice_ptr lattice, int subs)
{
  return &( lattice->bcs_in[subs].num_pressure_n_in0);
}

real* pressure_s_in0( lattice_ptr lattice, int subs)
{
  return lattice->bcs_in[subs].pressure_s_in0;
}

real** pressure_s_in0_ptr( lattice_ptr lattice, int subs)
{
  return &( lattice->bcs_in[subs].pressure_s_in0);
}

int num_pressure_s_in0( lattice_ptr lattice, int subs)
{
  return lattice->bcs_in[subs].num_pressure_s_in0;
}

int* num_pressure_s_in0_ptr( lattice_ptr lattice, int subs)
{
  return &( lattice->bcs_in[subs].num_pressure_s_in0);
}


//-------------------------------------------------------------//
// From Dr Dobbs "CUDA: Supercomputing for the masses, Part 3" //
// http://drdobbs.com/architecture-and-design/207200659        //
//-------------------------------------------------------------//
void checkCUDAError(const char *file, int line, const char *msg)
{
#if CUDA_ERROR_REPORTING
  // Don't forget to do cudaThreadSynchronize(); before
  // using this function on a kernel

  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err) 
  {
    fprintf(stderr, " CUDA Error in file:  %s\n"
        " Line number:         %d\n"
        " Point of failure:    %s\n"
        " Error message:       %s\n"
        , file, line, msg 
        , cudaGetErrorString( err) );

    process_exit(EXIT_FAILURE);
  }                   
#endif      
}

// }
//##############################################################################
