#include <stdio.h>

// Convert little endian to big endian.
// c.f. http://www.netrino.com/Publications/Glossary/Endianness.html
#define LE2BE(A) ( (((short int)(A) & 0xff00) >> 8) | \
                   (((short int)(A) & 0x00ff) << 8))

main()
{
  short int   ni, nj, nk, frame, subs;
  char  filename[1024];
  char  filename_raw[1024];
  char  filename_df3[1024];
  FILE *raw, *df3;
  int   size;
  int   i, j, k;
  unsigned char *val;
  int   n;

  // Big endians
  short    int  *be_int;

  ni    = 34;
  nj    = 34;
  nk    = 34;
  frame = 32;
  subs  = 0;

  be_int  = (short    int *)malloc( sizeof(short    int ));
  val = (unsigned char*)malloc(sizeof(unsigned char));

 for( frame=0; frame<=128; frame++)
 {
   printf("\n");

//printf("%s %d >> 8*sizeof(short int) = %d\n",
//      __FILE__,__LINE__, 8*sizeof(short int));
//printf("%s %d >> 8*sizeof(char) = %d\n",
//      __FILE__,__LINE__, 8*sizeof(char));

  sprintf( filename,
          "rho_%dx%dx%d_frame%04d_subs%02d",
           ni, nj, nk, frame, subs);

  sprintf( filename_raw, "%s.raw", filename);
  sprintf( filename_df3, "%s.df3", filename);

  raw = fopen( filename_raw, "r");
        if(!raw) { printf("%s %d >> ERROR: fopen() = %d\n",
                     __FILE__,__LINE__, raw); process_exit(1); }

  df3 = fopen( filename_df3, "w");
        if(!df3) { printf("%s %d >> ERROR: fopen() = %d\n",
                     __FILE__,__LINE__, df3); process_exit(1); }

  printf("%s %d >> Writing 6-byte header to file \"%s\"\n",
    __FILE__,__LINE__, filename_df3);

  *be_int = LE2BE(ni);
  size = fwrite( (short int *)be_int, 2, 1, df3);
         if(size!=1) { printf("%s %d >> ERROR: fwrite() = %d\n",
                        __FILE__,__LINE__, size); process_exit(1); }
  *be_int = LE2BE(nj);
  size = fwrite( (short int *)be_int, 2, 1, df3);
         if(size!=1) { printf("%s %d >> ERROR: fwrite() = %d\n",
                        __FILE__,__LINE__, size); process_exit(1); }
  *be_int = LE2BE(nk);
  size = fwrite( (short int *)be_int, 2, 1, df3);
         if(size!=1) { printf("%s %d >> ERROR: fwrite() = %d\n",
                        __FILE__,__LINE__, size); process_exit(1); }

  printf("%s %d >> "
      "Copying %d*%d*%d=%d 1-byte values from file \"%s\" to file \"%s\"\n",
    __FILE__,__LINE__, ni, nj, nk, ni*nj*nk, filename_raw, filename_df3);

  n = 0;
  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
//printf("%s %d >> (%2d,%2d,%2d)\n",
//  __FILE__,__LINE__, i, j, k);

        size = fread(  (char*)val, 1, 1, raw);
               if(size!=1) { printf("%s %d >> ERROR: fread() = %d\n",
                              __FILE__,__LINE__, size); process_exit(1); }
        size = fwrite( (char*)val, 1, 1, df3);
               if(size!=1) { printf("%s %d >> ERROR: fwrite() = %d\n",
                              __FILE__,__LINE__, size); process_exit(1); }
        n++;
      }
    }
  }

  fclose(raw);
  fclose(df3);
 }

  printf("\n");

  free( be_int);
  free( val);

  printf("%s %d >> Done! (n=%d)\n", __FILE__,__LINE__, n);

} /* main() */
