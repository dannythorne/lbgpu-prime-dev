//##############################################################################
//
// bcs.c
//
//  - Boundary conditions.
//
//  - void bcs( lattice_ptr lattice);
//
//  - void process_bcs( char *filename, int **bcs);
//

#define RHO0_TEST 0

#if SAVE_MEMO
//                                   B C S
//##############################################################################
//
// Apply boundary conditions.
//
void bcs( lattice_ptr lattice)
{
#if 0
  int    i, j, k, n, n1, a;
  int    ni, nj, nk;
  int    subs;
  real *ftemp, *ftemp1, *feq, *feq1, temp[Q];
  real  v, rho_in, rho_out, rho1;
  real  c0;
  real  D;
  real  c2;

  real u_x,  u_y, u_z, usq,udotx;
  real u_in[2][2],
       u_out[2][2],
       u,
       rho;
  real c;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;
  ftemp  = (real*)malloc(19*sizeof(real));
  ftemp1 = (real*)malloc(19*sizeof(real));
  feq    = (real*)malloc(19*sizeof(real));
  feq1   = (real*)malloc(19*sizeof(real));

  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
    // GZL P R E S S U R E // V E L O C I T Y   T O P   I N F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //  -- Velocity boundary on top side using inflow velocity condition.
    if( lattice->param.GZL && lattice->param.PressureBC)//lattice->param.velocity_t_in[subs])
    {
      k = nk-1;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          ftemp = lattice->pdf[subs][n].ftemp;

#if STORE_UEQ

#else /* !( STORE_UEQ) */
          //       lattice->macro_vars[subs][n].u[0] =  0.02 ;
          //       lattice->macro_vars[subs][n].u[1] =  0. ;
          //       lattice->macro_vars[subs][n].u[2] =  0. ;
#endif /* STORE_UEQ */

          n1 = IJK2N( i, j, k-1, ni, nj);
          ftemp1 = lattice->pdf[subs][n1].ftemp;

          switch(NUM_FLUID_COMPONENTS)
          {
            case 1:
              {
                rho =  /*lattice->macro_vars[subs][n1].rho;*/lattice->param.rho_in;
                u_x =  0.;//lattice->macro_vars[subs][n1].u[0] ;
                u_y =  0.;//lattice->macro_vars[subs][n1].u[1] ;
                u_z =  lattice->macro_vars[subs][n1].u[2] ;
                lattice->macro_vars[subs][n].u[2] = u_z;
                break;
              }
            case 2:
              {
#if STORE_UEQ
                if(subs==0)
                {
                  rho = lattice->param.rho_in;
                  //     rho = 8+ 0.001*(real)lattice->time;//lattice->param.rho_in ;    //lattice->param.uz_in; subs= 0 Non-wetting
                  //     if( lattice->time >1000)  {rho = 9.2;}
                }
                if(subs==1) rho = 0.001 ;    //lattice->param.uz_in; subs= 1 Wetting

                //          rho =  lattice->param.rho_A[subs];
                u_x =  0.;// lattice->ueq[n1].u[0] ;
                u_y =  0.;// lattice->ueq[n1].u[1] ;
                u_z =  lattice->ueq[n1].u[2] ;
                lattice->macro_vars[subs][n].u[2] = u_z;
                lattice->ueq[n].u[2] = u_z;

#endif
                break;
              }
          }
          lattice->macro_vars[subs][n].rho = rho;
          lattice->macro_vars[subs][n].u[0] = u_x;
          lattice->macro_vars[subs][n].u[1] = u_y;

          //rho1 = lattice->macro_vars[subs][n1].rho
          //u1   = lattice->macro_vars[subs][n1].u[0];
          //v1   = lattice->macro_vars[subs][n1].u[1];
          //w1   = lattice->macro_vars[subs][n1].u[2];

          if( !lattice->solids[subs][n].is_solid)
          {
            compute_a_feq( feq, rho, u_x, u_y, u_z );
            compute_single_feq(  lattice, n1, subs, feq1);
            for (a= 0; a<Q; a++)
            {
              temp[a] =  ftemp1[a] - feq1[a];
              ftemp[a]=  feq[a]+ temp[a];
            }


          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_GZL) */

    //GZL P R E S S U R E // V E L O C I T Y   B O T T O M   O U T F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // velocity bottom inflow
    //  -- Velocity boundary on bottom side using inflow velocity condition.
    if(  lattice->param.GZL && lattice->param.PressureBC)//lattice->param.velocity_b_in[subs])
    {

      k = 0;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          ftemp = lattice->pdf[subs][n].ftemp;

          n1 = IJK2N( i, j, k+1, ni, nj);
          ftemp1   = lattice->pdf[subs][n1].ftemp;

          switch(NUM_FLUID_COMPONENTS)
          {
            case 1:
              {
                rho =  /*lattice->macro_vars[subs][n1].rho;*/lattice->param.rho_out;
                u_x =  0.;//lattice->macro_vars[subs][n1].u[0] ;
                u_y =  0.;//lattice->macro_vars[subs][n1].u[1] ;
                u_z =  lattice->macro_vars[subs][n1].u[2] ;
                lattice->macro_vars[subs][n].u[2] = u_z;
                break;
              }
            case 2:
              {
#if STORE_UEQ
                if(subs==0) rho = 0.001;
                if(subs==1)
                {
                  rho = lattice->param.rho_out;
                  //rho = 8.0- 0.001*(real)lattice->time;//
                  //if( lattice->time >1000)   {rho = 6.8;}
                }
                //          rho =  lattice->param.rho_B[subs];
                u_x =  0.;//lattice->ueq[n1].u[0] ;
                u_y =  0.;//lattice->ueq[n1].u[1] ;
                u_z =  lattice->ueq[n1].u[2] ;
                lattice->macro_vars[subs][n].u[2] = u_z;
                lattice->ueq[n].u[2] = u_z;
#endif
                break;
              }
          }
          lattice->macro_vars[subs][n].rho = rho;
          lattice->macro_vars[subs][n].u[0] = u_x;
          lattice->macro_vars[subs][n].u[1] = u_y;

          //rho1 = lattice->macro_vars[subs][n1].rho;
          //u1   = lattice->macro_vars[subs][n1].u[0];
          //v1   = lattice->macro_vars[subs][n1].u[1];
          //w1   = lattice->macro_vars[subs][n1].u[2];

          if( !lattice->solids[subs][n].is_solid)
          {
            compute_a_feq( feq, rho, u_x, u_y, u_z );
            compute_single_feq(  lattice, n1, subs, feq1);

            for (a= 0; a<Q; a++)
            {
              temp[a] =  ftemp1[a] - feq1[a];
              ftemp[a]=  feq[a]+ temp[a];
            }

          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_GZL) */
  } /*for SUBS*/

#else
  // TODO
#endif
}
#else
//                                   BCS_1
//##############################################################################
//
// Apply boundary conditions before stream_collide_stream.
//
void bcs_1( lattice_ptr lattice, real *f_mem_d, real *mv_mem_d, unsigned char *solids_mem_d)
{
#ifdef __CUDACC__
  int subs;

  dim3 blockXboundDim(1, 1, 1);
  dim3 gridXboundDim(1, 1, 1);

  if( get_NumDims( lattice) == 2)
  { 
    blockXboundDim.y = get_BY( lattice); 
    gridXboundDim.y = get_LY( lattice) / get_BY( lattice); 
  }
  if( get_NumDims( lattice) == 3)
  {  
    blockXboundDim.y = get_BY( lattice); 
    gridXboundDim.y = get_LY( lattice) / get_BY( lattice); 
    blockXboundDim.z = get_BZ( lattice); 
    gridXboundDim.z = get_LZ( lattice) / get_BZ( lattice); 
  }

  dim3 blockYboundDim(1, 1, 1);
  dim3 gridYboundDim(1, 1, 1);

  if( get_NumDims( lattice) == 2)
  { 
    blockYboundDim.x = get_BX( lattice); 
    gridYboundDim.x = get_LX( lattice) / get_BX( lattice); 
  }
  if( get_NumDims( lattice) == 3)
  {  
    blockYboundDim.x = get_BX( lattice); 
    gridYboundDim.x = get_LX( lattice) / get_BX( lattice); 
    blockYboundDim.z = get_BZ( lattice); 
    gridYboundDim.z = get_LZ( lattice) / get_BZ( lattice); 
  }

  dim3 blockZboundDim(1, 1, 1);
  dim3 gridZboundDim(1, 1, 1);

  if( get_NumDims( lattice) == 3)
  { 
    blockZboundDim.x = get_BX( lattice); 
    gridZboundDim.x = get_LX( lattice) / get_BX( lattice); 
    blockZboundDim.y = get_BY( lattice); 
    gridZboundDim.y = get_LY( lattice) / get_BY( lattice); 
  }


  for( subs=0; subs<get_NumSubs( lattice); subs++)
  {
    cudaMemcpyToSymbol( subs_c, &subs, sizeof(int));

    checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
    // PRESSURE NORTH IN
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Constant pressure Zou-He boundary condition on north face

    if( get_proc_id( lattice) == get_num_procs( lattice)-1
        && lattice->param.pressure_n_in[subs] )
    {
      real rho;

      /* if( lattice->param.pressure_n_in[subs]==2)
         {
         rho = *( pressure_n_in0( lattice, subs)
         + get_time(lattice)%num_pressure_n_in0(lattice,subs));
         }
         else
         {*/
      rho = lattice->param.rho_in;
      //}

      cudaMemcpyToSymbol( fixed_bound_var_c
          , &rho, sizeof(real));
      checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");

      // Kernel for setting system boundary conditions
      k_sysbound_pressure_n_1
        <<<
        gridYboundDim
        , blockYboundDim
        , sizeof(real)
        * ( get_NumVelDirs( lattice) + 1 )
        * get_BX( lattice)
        * get_BZ( lattice)
        >>>( f_mem_d, mv_mem_d, solids_mem_d);
      cudaThreadSynchronize();
      checkCUDAError( __FILE__, __LINE__, "k_sysbound");

    } /* if( lattice->param.pressure_n_in[subs] ) */

    //********************************************************************


    // PRESSURE NORTH OUT
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Constant pressure Zou-He boundary condition on north face

    if( get_proc_id( lattice) == get_num_procs( lattice)-1
        && lattice->param.pressure_n_out[subs] )
    {
      cudaMemcpyToSymbol( fixed_bound_var_c
          , &(lattice->param.rho_out), sizeof(real));

      checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
      // Kernel for setting system boundary conditions
      k_sysbound_pressure_n_1
        <<<
        gridYboundDim
        , blockYboundDim
        , sizeof(real)
        * ( get_NumVelDirs( lattice) + 1 )
        * get_BX( lattice)
        * get_BZ( lattice)
        >>>( f_mem_d, mv_mem_d, solids_mem_d);
      cudaThreadSynchronize();
      checkCUDAError( __FILE__, __LINE__, "k_sysbound");

    } /* if( lattice->param.pressure_n_in[subs] ) */

    //********************************************************************


    // PRESSURE SOUTH IN
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Constant pressure Zou-He boundary condition on south face

    if( get_proc_id( lattice) == 0
        && lattice->param.pressure_s_in[subs] )
    {

      real rho;
      /*if( lattice->param.pressure_s_in[subs]==2)
        {
        rho = *( pressure_s_in0( lattice, subs)
        + get_time(lattice)%num_pressure_s_in0(lattice,subs));
        }
        else
        {*/
      rho = lattice->param.rho_in;
      //}

      cudaMemcpyToSymbol( fixed_bound_var_c
          , &rho, sizeof(real));
      checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");

      // Kernel for setting system boundary conditions
      k_sysbound_pressure_s_1
        <<<
        gridYboundDim
        , blockYboundDim
        , sizeof(real)
        * ( get_NumVelDirs( lattice) + 1 )
        * get_BX( lattice)
        * get_BZ( lattice)
        >>>( f_mem_d, solids_mem_d);
      cudaThreadSynchronize();
      checkCUDAError( __FILE__, __LINE__, "k_sysbound");

    } /* if( lattice->param.pressure_n_in[subs] ) */

    //********************************************************************


    // PRESSURE SOUTH OUT
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Constant pressure Zou-He boundary condition on south face

    if( get_proc_id( lattice) == 0
        && lattice->param.pressure_s_out[subs] )
    {
      cudaMemcpyToSymbol( fixed_bound_var_c
          , &(lattice->param.rho_out), sizeof(real));

      checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
      // Kernel for setting system boundary conditions
      k_sysbound_pressure_s_1
        <<<
        gridYboundDim
        , blockYboundDim
        , sizeof(real)
        * ( get_NumVelDirs( lattice) + 1 )
        * get_BX( lattice)
        * get_BZ( lattice)
        >>>( f_mem_d, solids_mem_d);
      cudaThreadSynchronize();
      checkCUDAError( __FILE__, __LINE__, "k_sysbound");

    } /* if( lattice->param.pressure_n_in[subs] ) */

    //********************************************************************

  }  // for( subs=0; subs<get_NumSubs( lattice); subs++)

#else
  int    i, j, k, n, n1, a, id;
  int    ni, nj, nk;
  int    subs;
  real *ftemp, *f, *f1, *feq, *feq1, *rhoo, temp[Q];
  real  v, rho_in, rho_out, rho1;
  real  c0;
  real  D;
  real  c2;

#if RHO0_TEST
  //------------------------------------------------------------------[ TEST ]----
  real *rho0;
  //------------------------------------------------------------------[ TEST ]----
#endif /* RHO0_TEST */
  real u_x,  u_y, u_z, usq,udotx;
  real u_in[2][2],
       u_out[2][2],
       u,
       rho,
       hua;
  real c;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  id = get_proc_id(lattice);

  // for( subs=0; subs<(NUM_FLUID_COMPONENTS)-(INAMURO_SIGMA_COMPONENT); subs++)
  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
    // P R E S S U R E   T O P   I N F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // pressure top inflow
    //  -- Pressure boundary on top side using inflow pressure condition.
    if(   (id ==(get_num_procs(lattice)-1))
        && !lattice->param.GZL
        && lattice->param.pressure_t_in[subs] )
    {
      k = nk-1;
      switch(NUM_FLUID_COMPONENTS)
      {
        case 1:
          { rho_in = lattice->param.rho_in; break;}
        case 2:
          {
            if(subs==0) { rho_in = lattice->param.rho_in;}
            if(subs==1) { rho_in = lattice->param.rho_A[subs];}
            break;
          }
      }

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          ftemp = lattice->pdf[subs][n].ftemp;

          // Top, Inflow
          if( lattice->param.incompressible)
          {
            u_z = -rho_in
              + ( ftemp[ C]
                  + ftemp[ N] + ftemp[ S] + ftemp[ E] + ftemp[ W]
                  + 2.*( ftemp[T ]
                    + ftemp[TW] + ftemp[TE] + ftemp[TS] + ftemp[TN]));
            c = u_z;
          }
          else // compressible
          {
            u_z = -1.
              + ( ftemp[C]
                  + ftemp[N ] + ftemp[S ] + ftemp[ E] + ftemp[ W]
                  + ftemp[NE] + ftemp[NW] + ftemp[SE] + ftemp[SW]
                  + 2.*( ftemp[T ]
                    + ftemp[TW] + ftemp[TE] + ftemp[TS] + ftemp[TN]))
              / rho_in;
            c = u_z*rho_in;
          }
          //rev_Huang
          ftemp[B ] = ftemp[T ] - (1./3.)*c;
          ftemp[BW] = ftemp[TE] - (1./6.)*c
            + 0.5*(-ftemp[ W] - ftemp[NW] - ftemp[SW]
                + ftemp[ E] + ftemp[NE] + ftemp[SE]);
          ftemp[BE] = ftemp[TW] - (1./6.)*c
            - 0.5*(-ftemp[ W] - ftemp[NW] - ftemp[SW]
                + ftemp[ E] + ftemp[NE] + ftemp[SE]);
          ftemp[BS] = ftemp[TN] - (1./6.)*c
            + 0.5*( ftemp[N ] + ftemp[NW] + ftemp[NE]
                - ftemp[S ] - ftemp[SW] - ftemp[SE]);
          ftemp[BN] = ftemp[TS] - (1./6.)*c
            - 0.5*( ftemp[N ] + ftemp[NW] + ftemp[NE]
                - ftemp[S ] - ftemp[SW] - ftemp[SE]);
#if 1 // Confirm density

          rho = ftemp[ C]
            + ftemp[ E] + ftemp[ W]
            + ftemp[ N] + ftemp[ S]
            + ftemp[ T] + ftemp[ B]
            + ftemp[NE] + ftemp[NW] + ftemp[SE] + ftemp[SW]
            + ftemp[TN] + ftemp[TS] + ftemp[BN] + ftemp[BS]
            + ftemp[TE] + ftemp[TW] + ftemp[BE] + ftemp[BW];
          if( rho_in - rho > 1e-6)
          {
            printf("%s %d >> ERROR: pressure_top_in FAIL! "
                "rho_in = %f, rho = %f\n",__FILE__,__LINE__, rho_in, rho);
          }
#endif

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.pressure_t_in[subs] ) */


    //********************************************************************************


    // P R E S S U R E   B O T T O M   O U T F L O W   B C  (peter's attempt?)
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // pressure bottom outflow
    //  -- Pressure boundary on bottom side using outflow pressure condition.

    if((id ==0) && !lattice->param.GZL && lattice->param.pressure_b_out[subs] )
    {
      k = 0;

      switch(NUM_FLUID_COMPONENTS)
      {
        case 1:
          {  rho_out = lattice->param.rho_out; break; }
        case 2:
          {
            if(subs==0) { rho_out = lattice->param.rho_out;}
            if(subs==1) { rho_out = lattice->param.rho_A[subs];}
            break;
          }
      }

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          ftemp = lattice->pdf[subs][n].ftemp;

          // Bottom, rho_out
          if( lattice->param.incompressible)
          {
            u_z = rho_out
              - ( ftemp[C]
                  + ftemp[N] + ftemp[S] + ftemp[E] + ftemp[W]
                  + 2.*( ftemp[B]
                    + ftemp[BW] + ftemp[BE] + ftemp[BS] + ftemp[BN]));
            c = u_z;
          }
          else // compressible
          {
            u_z =  1.
              - ( ftemp[C]
                  + ftemp[N ] + ftemp[S ] + ftemp[ E] + ftemp[ W]
                  + ftemp[NW] + ftemp[NE] + ftemp[SW] + ftemp[SE]
                  + 2.*( ftemp[B ]
                    + ftemp[BW] + ftemp[BE] + ftemp[BS] + ftemp[BN]))
              / rho_out;
            c = u_z*rho_out;
          }

          ftemp[T ] = ftemp[B ] + (1./3.)*c;
          ftemp[TW] = ftemp[BE] + (1./6.)*c
            + 0.5*(-ftemp[ W] - ftemp[NW] - ftemp[SW]
                + ftemp[ E] + ftemp[NE] + ftemp[SE]);
          ftemp[TE] = ftemp[BW] + (1./6.)*c
            - 0.5*(-ftemp[ W] - ftemp[NW] - ftemp[SW]
                + ftemp[ E] + ftemp[NE] + ftemp[SE]);
          ftemp[TS] = ftemp[BN] + (1./6.)*c
            + 0.5*( ftemp[N ] + ftemp[NW] + ftemp[NE]
                - ftemp[S ] - ftemp[SW] - ftemp[SE]);
          ftemp[TN] = ftemp[BS] + (1./6.)*c
            - 0.5*( ftemp[N ] + ftemp[NW] + ftemp[NE]
                - ftemp[S ] - ftemp[SW] - ftemp[SE]);

#if 1 // Confirm density
          rho = ftemp[ C]
            + ftemp[ E] + ftemp[ W]
            + ftemp[ N] + ftemp[ S]
            + ftemp[ T] + ftemp[ B]
            + ftemp[NE] + ftemp[NW] + ftemp[SE] + ftemp[SW]
            + ftemp[TN] + ftemp[TS] + ftemp[BN] + ftemp[BS]
            + ftemp[TE] + ftemp[TW] + ftemp[BE] + ftemp[BW];
          if( rho_out - rho > 1e-6)
          {
            printf("%s %d >> ERROR: pressure_bottom_out FAIL! "
                "rho_out = %f, rho = %f\n",__FILE__,__LINE__, rho_out, rho);
          }
#endif

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.pressure_b_out[subs] ) */


    //********************************************************************************
    // V E L O C I T Y   T O P   I N F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // velocity top inflow
    //  -- Velocity boundary on top side using inflow velocity condition.
    if( id == get_num_procs( lattice) - 1 && !lattice->param.GZL && 
        lattice->param.velocity_t_in[subs])
    {

      //printf("%s %d >> BOOM!",__FILE__,__LINE__);
      if(subs==0) u = lattice->param.uz_in;    //lattice->param.uz_in; Non-wetting
      if(subs==1) u = lattice->param.uz_in;    //lattice->param.uz_in;

      k = nk-1;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          ftemp = lattice->pdf[subs][n].ftemp;

          if( !lattice->solids[subs][n].is_solid)
          {
            rho = ( ftemp[C ]
                + ftemp[W ] + ftemp[E ] + ftemp[N ] + ftemp[S ]
                + ftemp[NW] + ftemp[NE] + ftemp[SW] + ftemp[SE]
                + 2.*( ftemp[T ]
                  + ftemp[TW] + ftemp[TE] + ftemp[TS] + ftemp[TN])
                )
              / ( 1. + u);
            c = rho*u;

            ftemp[B ] = ftemp[T ] - (1./3.)*c;
            ftemp[BW] = ftemp[TE] - (1./6.)*c + 0.5* (-ftemp[W ]- ftemp[NW]- ftemp[SW]
                + ftemp[E] +ftemp[NE] + ftemp[SE]);
            ftemp[BE] = ftemp[TW] - (1./6.)*c - 0.5* (-ftemp[W ]- ftemp[NW]- ftemp[SW]
                + ftemp[E] +ftemp[NE] + ftemp[SE]);
            ftemp[BS] = ftemp[TN] - (1./6.)*c + 0.5* (ftemp[N ] + ftemp[NW] + ftemp[NE]
                -ftemp[S ] -ftemp[SW] - ftemp[SE]);
            ftemp[BN] = ftemp[TS] - (1./6.)*c - 0.5* (ftemp[N ] + ftemp[NW] + ftemp[NE]
                -ftemp[S ] -ftemp[SW] - ftemp[SE]);

          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_t_in[subs]) */

    // V E L O C I T Y   B O T T O M   O U T F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // velocity bottom inflow
    //  -- Velocity boundary on bottom side using outflow velocity condition.
    if( id == 0 && !lattice->param.GZL && lattice->param.velocity_b_out[subs])
    {
      //printf("%s %d >> BOOM!",__FILE__,__LINE__);
      if(subs==0) u = lattice->param.uz_out;    //lattice->param.uz_out;  Non-wetting
      if(subs==1) u = lattice->param.uz_out;    //lattice->param.uz_out;

      k = 0;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          ftemp = lattice->pdf[subs][n].ftemp;

          if( !lattice->solids[subs][n].is_solid)
          {
            rho = ( ftemp[C ]
                + ftemp[W ] + ftemp[E ] + ftemp[N ] + ftemp[S ]
                + ftemp[NW] + ftemp[NE] + ftemp[SW] + ftemp[SE]
                + 2.*( ftemp[B ]
                  + ftemp[BW] + ftemp[BE] + ftemp[BS] + ftemp[BN])
                )
              / ( 1. - u);
            c = rho*u;


            ftemp[T ] = ftemp[B ] + (1./3.)*c;
            ftemp[TW] = ftemp[BE] + (1./6.)*c + 0.5* (-ftemp[W ]- ftemp[NW]- ftemp[SW]
                + ftemp[E] +ftemp[NE] + ftemp[SE]);
            ftemp[TE] = ftemp[BW] + (1./6.)*c - 0.5* (-ftemp[W ]- ftemp[NW]- ftemp[SW]
                + ftemp[E] +ftemp[NE] + ftemp[SE]);
            ftemp[TS] = ftemp[BN] + (1./6.)*c  + 0.5* (ftemp[N ] + ftemp[NW] + ftemp[NE]
                -ftemp[S ] -ftemp[SW] - ftemp[SE]);
            ftemp[TN] = ftemp[BS] + (1./6.)*c  - 0.5* (ftemp[N ] + ftemp[NW] + ftemp[NE]
                -ftemp[S ] -ftemp[SW] - ftemp[SE]);


          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_b_in[subs]) */

    //********************************************************************************


    // GZL P R E S S U R E // V E L O C I T Y   T O P   I N F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //  -- Velocity boundary on top side using inflow velocity condition.
#if PARALLEL
    id = get_proc_id(lattice);
#else
    id = get_num_procs(lattice)-1;
#endif

    if( (id ==get_num_procs(lattice)-1) && lattice->param.GZL && lattice->param.PressureBC)//lattice->param.velocity_t_in[subs])
    {
      printf("%s %d >> BOOM!",__FILE__,__LINE__);
      k = nk-1;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          f = lattice->pdf[subs][n].f;

#if STORE_UEQ
          //       u_x = &lattice->ueq[n].u[0] = 0.02;
          //       lattice->ueq[n].u[1] = 0.;
          //       lattice->ueq[n].u[2] = 0.;
#else /* !( STORE_UEQ) */
          //       lattice->macro_vars[subs][n].u[0] =  0.02 ;
          //       lattice->macro_vars[subs][n].u[1] =  0. ;
          //       lattice->macro_vars[subs][n].u[2] =  0. ;
#endif /* STORE_UEQ */

          n1 = IJK2N( i, j, k-1, ni, nj);
          f1 = lattice->pdf[subs][n1].f;
          feq1 = lattice->pdf[subs][n1].feq;

          //  rhoo = &( lattice->macro_vars[subs][n].rho);
          switch(NUM_FLUID_COMPONENTS)
          {case 1:
            {
              rho =  /*lattice->macro_vars[subs][n1].rho;*/lattice->param.rho_in;
              u_x =  0.;//lattice->macro_vars[subs][n1].u[0] ;
              u_y =  0.;//lattice->macro_vars[subs][n1].u[1] ;
              u_z =  lattice->macro_vars[subs][n1].u[2] ;
              lattice->macro_vars[subs][n].u[2] = u_z;
              break;
            }
            case 2:
            {
#if STORE_UEQ
              if(subs==0)
              {
                rho =lattice->param.rho_in;// 8+ 0.001*(real)lattice->time;// ;    subs= 0 Non-wetting
                if( lattice->time >1000)  {rho = lattice->param.rho_in;}
              }
              if(subs==1) rho = 0.000 ;    //lattice->param.uz_in; subs= 1 Wetting

              //          rho =  lattice->param.rho_A[subs];
              u_x =  0.;// lattice->ueq[n1].u[0] ;
              u_y =  0.;// lattice->ueq[n1].u[1] ;
              u_z =  lattice->ueq[n1].u[2] ;
              lattice->macro_vars[subs][n].u[2] = u_z;
              lattice->ueq[n].u[2] = u_z;

#endif
              break;
            }
          }
          lattice->macro_vars[subs][n].rho = rho;
          lattice->macro_vars[subs][n].u[0] = u_x;
          lattice->macro_vars[subs][n].u[1] = u_y;


          if( !lattice->solids[subs][n].is_solid)
          {
            compute_a_feq( feq, rho, u_x, u_y, u_z );

            for (a= 0; a<Q; a++)
            {
              temp[a] =  f1[a] - feq1[a];
              f[a]=feq[a]+ temp[a];
            }


          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_GZL) */

    //GZL P R E S S U R E // V E L O C I T Y   B O T T O M   O U T F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // velocity bottom inflow
    //  -- Velocity boundary on bottom side using inflow velocity condition.
#if PARALLEL
    id = get_proc_id(lattice);
#else
    id = 0;
#endif
    if( (id ==0) && lattice->param.GZL && lattice->param.PressureBC)//lattice->param.velocity_b_in[subs])
    {
      printf("%s %d >> BOOM!",__FILE__,__LINE__);
      k = 0;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          f = lattice->pdf[subs][n].f;

          n1 = IJK2N( i, j, k+1, ni, nj);
          f1   = lattice->pdf[subs][n1].f;
          feq1 = lattice->pdf[subs][n1].feq;

          switch(NUM_FLUID_COMPONENTS)
          {case 1:
            {
              rho =  /*lattice->macro_vars[subs][n1].rho;*/lattice->param.rho_out;
              u_x =  0.;//lattice->macro_vars[subs][n1].u[0] ;
              u_y =  0.;//lattice->macro_vars[subs][n1].u[1] ;
              u_z =  lattice->macro_vars[subs][n1].u[2] ;
              lattice->macro_vars[subs][n].u[2] = u_z;
              break;
            }
            case 2:
            {
#if STORE_UEQ
              if(subs==0) rho = 0.000;
              if(subs==1)
              {
                rho =lattice->param.rho_out; //8.0- 0.001*(real)lattice->time;//
                if( lattice->time >6000)   {rho = lattice->param.rho_out;}
              }
              //          rho =  lattice->param.rho_B[subs];
              u_x =  0.;//lattice->ueq[n1].u[0] ;
              u_y =  0.;//lattice->ueq[n1].u[1] ;
              u_z =  lattice->ueq[n1].u[2] ;
              lattice->macro_vars[subs][n].u[2] = u_z;
              lattice->ueq[n].u[2] = u_z;
#endif
              break;
            }
          }
          lattice->macro_vars[subs][n].rho = rho;
          lattice->macro_vars[subs][n].u[0] = u_x;
          lattice->macro_vars[subs][n].u[1] = u_y;



          if( !lattice->solids[subs][n].is_solid)
          {
            compute_a_feq( feq, rho, u_x, u_y, u_z );

            for (a= 0; a<Q; a++)
            {
              temp[a] =f1[a] - feq1[a];
              f[a] =  feq[a]+ temp[a];
            }

          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_GZL) */



    //PAN   V E L O C I T Y   T O P   I N F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //  -- Velocity boundary on top side using inflow velocity condition.
    if( lattice->param.GZL && (!lattice->param.PressureBC))//lattice->param.velocity_t_in[subs])
    {
      printf("%s %d >> BOOM!",__FILE__,__LINE__);
      k = nk-1;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          f = lattice->pdf[subs][n].f;
          //  feq =  lattice->pdf[subs][n].feq;


          n1 = IJK2N( i, j, k-1, ni, nj);
          f1 = lattice->pdf[subs][n1].f;
          feq1 = lattice->pdf[subs][n1].feq;

          rhoo = &( lattice->macro_vars[subs][n].rho);
          switch(NUM_FLUID_COMPONENTS)
          {case 1:
            {
              rho =  /*lattice->macro_vars[subs][n1].rho;*/lattice->param.rho_in;
              u_x =  0.;//lattice->macro_vars[subs][n1].u[0] ;
              u_y =  0.;//lattice->macro_vars[subs][n1].u[1] ;
              u_z =  -0.0;//lattice->macro_vars[subs][n1].u[2] ;
              break;
            }
            case 2:
            {
#if STORE_UEQ
              rho    =  lattice->param.rho_A[subs];
              u_x = 0.;// lattice->ueq[n1].u[0] ;
              u_y = 0.;// lattice->ueq[n1].u[1] ;
              u_z = 0.;// lattice->ueq[n1].u[2] ;
#endif
              break;
            }
          }
          lattice->macro_vars[subs][n].rho = rho;
          lattice->macro_vars[subs][n].u[0] = u_x;
          lattice->macro_vars[subs][n].u[1] = u_y;
          lattice->macro_vars[subs][n].u[2] = u_z;


          if( !lattice->solids[subs][n].is_solid)
          {
            compute_a_feq( feq, rho, u_x, u_y, u_z );

            for (a= 0; a<Q; a++)
            {f[a] = feq[a];}//+(f1[a] - feq1[a]);


          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_GZL) */

    //PAN V E L O C I T Y   B O T T O M   I N F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // velocity bottom inflow
    //  -- Velocity boundary on bottom side using inflow velocity condition.
    if(  lattice->param.GZL && (!lattice->param.PressureBC))//lattice->param.velocity_b_in[subs])
    {
      printf("%s %d >> BOOM!",__FILE__,__LINE__);
      k = 0;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          f = lattice->pdf[subs][n].f;
          //  feq =  lattice->pdf[subs][n].feq;

          n1 = IJK2N( i, j, k+1, ni, nj);
          f1 = lattice->pdf[subs][n1].f;
          feq1 = lattice->pdf[subs][n1].feq;


          switch(NUM_FLUID_COMPONENTS)
          {case 1:
            {
              rho =  /*lattice->macro_vars[subs][n1].rho;*/lattice->param.rho_out;
              u_x =  0.;//lattice->macro_vars[subs][n1].u[0] ;
              u_y =  0.;//lattice->macro_vars[subs][n1].u[1] ;
              u_z =  -0.0;//lattice->macro_vars[subs][n1].u[2] ;
              break;
            }
            case 2:
            {
#if STORE_UEQ
              rho =  lattice->param.rho_B[subs];
              u_x =  0.;//lattice->ueq[n1].u[0] ;
              u_y =  0.;//lattice->ueq[n1].u[1] ;
              u_z =  0.;//lattice->ueq[n1].u[2] ;
#endif
              break;
            }
          }

          lattice->macro_vars[subs][n].rho = rho;
          lattice->macro_vars[subs][n].u[0] = u_x;
          lattice->macro_vars[subs][n].u[1] = u_y;
          lattice->macro_vars[subs][n].u[2] = u_z;

          if( !lattice->solids[subs][n].is_solid)
          {
            compute_a_feq( feq, rho, u_x, u_y, u_z );
            if(subs==1)
              compute_a_feq( feq1, 2.0- 0.0003*(real)lattice->time, u_x, u_y, u_z );

            for (a= 0; a<Q; a++)
            {
              f[a] = feq[a];

            }// +(f1[a] - feq1[a]);

          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_GZL) */


  } /* for( subs=0; subs<(NUM_FLUID_COMPONENTS)-(INAMURO_SIGMA_COMPONENT); subs++) */


#endif
}

//                                   BCS_2
//##############################################################################
//
// Apply boundary conditions after stream_collide_stream.
//
void bcs_2( lattice_ptr lattice, real *f_mem_d, unsigned char *solids_mem_d)
{
#ifdef __CUDACC__
  int subs;

  dim3 blockXboundDim(1, 1, 1);
  dim3 gridXboundDim(1, 1, 1);

  if( get_NumDims( lattice) == 2)
  { 
    blockXboundDim.y = get_BY( lattice); 
    gridXboundDim.y = get_LY( lattice) / get_BY( lattice); 
  }
  if( get_NumDims( lattice) == 3)
  {  
    blockXboundDim.y = get_BY( lattice); 
    gridXboundDim.y = get_LY( lattice) / get_BY( lattice); 
    blockXboundDim.z = get_BZ( lattice); 
    gridXboundDim.z = get_LZ( lattice) / get_BZ( lattice); 
  }

  dim3 blockYboundDim(1, 1, 1);
  dim3 gridYboundDim(1, 1, 1);

  if( get_NumDims( lattice) == 2)
  { 
    blockYboundDim.x = get_BX( lattice); 
    gridYboundDim.x = get_LX( lattice) / get_BX( lattice); 
  }
  if( get_NumDims( lattice) == 3)
  {  
    blockYboundDim.x = get_BX( lattice); 
    gridYboundDim.x = get_LX( lattice) / get_BX( lattice); 
    blockYboundDim.z = get_BZ( lattice); 
    gridYboundDim.z = get_LZ( lattice) / get_BZ( lattice); 
  }

  dim3 blockZboundDim(1, 1, 1);
  dim3 gridZboundDim(1, 1, 1);

  if( get_NumDims( lattice) == 3)
  { 
    blockZboundDim.x = get_BX( lattice); 
    gridZboundDim.x = get_LX( lattice) / get_BX( lattice); 
    blockZboundDim.y = get_BY( lattice); 
    gridZboundDim.y = get_LY( lattice) / get_BY( lattice); 
  }


  for( subs=0; subs<get_NumSubs( lattice); subs++)
  {
    cudaMemcpyToSymbol( subs_c, &subs, sizeof(int));

    checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
    // PRESSURE NORTH IN
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Constant pressure Zou-He boundary condition on north face

    if( get_proc_id( lattice) == get_num_procs( lattice)-1
        && lattice->param.pressure_n_in[subs] )
    {
      real rho;

      /*if( lattice->param.pressure_n_in[subs]==2)
        {
        rho = *( pressure_n_in0( lattice, subs)
        + get_time(lattice)%num_pressure_n_in0(lattice,subs));
        }
        else
        {*/
      rho = lattice->param.rho_in;
      //}

      cudaMemcpyToSymbol( fixed_bound_var_c
          , &rho, sizeof(real));
      checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");

      // Kernel for setting system boundary conditions
      k_sysbound_pressure_n_2
        <<<
        gridYboundDim
        , blockYboundDim
        , sizeof(real)
        * ( get_NumVelDirs( lattice) + 1 )
        * get_BX( lattice)
        * get_BZ( lattice)
        >>>( f_mem_d, solids_mem_d);
      cudaThreadSynchronize();
      checkCUDAError( __FILE__, __LINE__, "k_sysbound");

    } /* if( lattice->param.pressure_n_in[subs] ) */

    //********************************************************************


    // PRESSURE NORTH OUT
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Constant pressure Zou-He boundary condition on north face

    if( get_proc_id( lattice) == get_num_procs( lattice)-1
        && lattice->param.pressure_n_out[subs] )
    {
      cudaMemcpyToSymbol( fixed_bound_var_c
          , &(lattice->param.rho_out), sizeof(real));

      checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
      // Kernel for setting system boundary conditions
      k_sysbound_pressure_n_2
        <<<
        gridYboundDim
        , blockYboundDim
        , sizeof(real)
        * ( get_NumVelDirs( lattice) + 1 )
        * get_BX( lattice)
        * get_BZ( lattice)
        >>>( f_mem_d, solids_mem_d);
      cudaThreadSynchronize();
      checkCUDAError( __FILE__, __LINE__, "k_sysbound");

    } /* if( lattice->param.pressure_n_in[subs] ) */

    //********************************************************************


    // PRESSURE SOUTH IN
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Constant pressure Zou-He boundary condition on south face

    if( get_proc_id( lattice) == 0
        && lattice->param.pressure_s_in[subs] )
    {
      real rho;
      /*if( lattice->param.pressure_s_in[subs]==2)
        {
        rho = *( pressure_s_in0( lattice, subs)
        + get_time(lattice)%num_pressure_s_in0(lattice,subs));
        }
        else
        {*/
      rho = lattice->param.rho_in;
      //}

      cudaMemcpyToSymbol( fixed_bound_var_c
          , &rho, sizeof(real));
      checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");

      // Kernel for setting system boundary conditions
      k_sysbound_pressure_s_2
        <<<
        gridYboundDim
        , blockYboundDim
        , sizeof(real)
        * ( get_NumVelDirs( lattice) + 1 )
        * get_BX( lattice)
        * get_BZ( lattice)
        >>>( f_mem_d, solids_mem_d);
      cudaThreadSynchronize();
      checkCUDAError( __FILE__, __LINE__, "k_sysbound");

    } /* if( lattice->param.pressure_n_in[subs] ) */

    //********************************************************************


    // PRESSURE SOUTH OUT
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Constant pressure Zou-He boundary condition on south face

    if( get_proc_id( lattice) == 0
        && lattice->param.pressure_s_out[subs] )
    {
      cudaMemcpyToSymbol( fixed_bound_var_c
          , &(lattice->param.rho_out), sizeof(real));

      checkCUDAError( __FILE__, __LINE__, "cudaMemcpyToSymbol");
      // Kernel for setting system boundary conditions
      k_sysbound_pressure_s_2
        <<<
        gridYboundDim
        , blockYboundDim
        , sizeof(real)
        * ( get_NumVelDirs( lattice) + 1 )
        * get_BX( lattice)
        * get_BZ( lattice)

        >>>( f_mem_d, solids_mem_d);
      cudaThreadSynchronize();
      checkCUDAError( __FILE__, __LINE__, "k_sysbound");

    } /* if( lattice->param.pressure_n_in[subs] ) */

    //********************************************************************


  }  // for( subs=0; subs<get_NumSubs( lattice); subs++)

#else
  int    i, j, k, n, n1, a, id;
  int    ni, nj, nk;
  int    subs;
  real *ftemp, *f, *f1, *feq, *feq1, *rhoo, temp[Q];
  real  v, rho_in, rho_out, rho1;
  real  c0;
  real  D;
  real  c2;

#if RHO0_TEST
  //------------------------------------------------------------------[ TEST ]----
  real *rho0;
  //------------------------------------------------------------------[ TEST ]----
#endif /* RHO0_TEST */
  real u_x,  u_y, u_z, usq,udotx;
  real u_in[2][2],
       u_out[2][2],
       u,
       rho,
       hua;
  real c;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  id = get_proc_id(lattice);

  // for( subs=0; subs<(NUM_FLUID_COMPONENTS)-(INAMURO_SIGMA_COMPONENT); subs++)
  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
    // P R E S S U R E   T O P   I N F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // pressure top inflow
    //  -- Pressure boundary on top side using inflow pressure condition.
    if(   (id ==(get_num_procs(lattice)-1))
        && !lattice->param.GZL
        && lattice->param.pressure_t_in[subs] )
    {
      k = nk-1;
      switch(NUM_FLUID_COMPONENTS)
      {
        case 1:
          { rho_in = lattice->param.rho_in; break;}
        case 2:
          {
            if(subs==0) { rho_in = lattice->param.rho_in;}
            if(subs==1) { rho_in = lattice->param.rho_A[subs];}
            break;
          }
      }

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          ftemp = lattice->pdf[subs][n].ftemp;

          // Top, Inflow
          if( lattice->param.incompressible)
          {
            u_z = -rho_in
              + ( ftemp[ C]
                  + ftemp[ N] + ftemp[ S] + ftemp[ E] + ftemp[ W]
                  + 2.*( ftemp[T ]
                    + ftemp[TW] + ftemp[TE] + ftemp[TS] + ftemp[TN]));
            c = u_z;
          }
          else // compressible
          {
            u_z = -1.
              + ( ftemp[C]
                  + ftemp[N ] + ftemp[S ] + ftemp[ E] + ftemp[ W]
                  + ftemp[NE] + ftemp[NW] + ftemp[SE] + ftemp[SW]
                  + 2.*( ftemp[T ]
                    + ftemp[TW] + ftemp[TE] + ftemp[TS] + ftemp[TN]))
              / rho_in;
            c = u_z*rho_in;
          }
          //rev_Huang
          ftemp[B ] = ftemp[T ] - (1./3.)*c;
          ftemp[BW] = ftemp[TE] - (1./6.)*c
            + 0.5*(-ftemp[ W] - ftemp[NW] - ftemp[SW]
                + ftemp[ E] + ftemp[NE] + ftemp[SE]);
          ftemp[BE] = ftemp[TW] - (1./6.)*c
            - 0.5*(-ftemp[ W] - ftemp[NW] - ftemp[SW]
                + ftemp[ E] + ftemp[NE] + ftemp[SE]);
          ftemp[BS] = ftemp[TN] - (1./6.)*c
            + 0.5*( ftemp[N ] + ftemp[NW] + ftemp[NE]
                - ftemp[S ] - ftemp[SW] - ftemp[SE]);
          ftemp[BN] = ftemp[TS] - (1./6.)*c
            - 0.5*( ftemp[N ] + ftemp[NW] + ftemp[NE]
                - ftemp[S ] - ftemp[SW] - ftemp[SE]);
#if 1 // Confirm density

          rho = ftemp[ C]
            + ftemp[ E] + ftemp[ W]
            + ftemp[ N] + ftemp[ S]
            + ftemp[ T] + ftemp[ B]
            + ftemp[NE] + ftemp[NW] + ftemp[SE] + ftemp[SW]
            + ftemp[TN] + ftemp[TS] + ftemp[BN] + ftemp[BS]
            + ftemp[TE] + ftemp[TW] + ftemp[BE] + ftemp[BW];
          if( rho_in - rho > 1e-6)
          {
            printf("%s %d >> ERROR: pressure_top_in FAIL! "
                "rho_in = %f, rho = %f\n",__FILE__,__LINE__, rho_in, rho);
          }
#endif

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.pressure_t_in[subs] ) */


    //********************************************************************************


    // P R E S S U R E   B O T T O M   O U T F L O W   B C  (peter's attempt?)
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // pressure bottom outflow
    //  -- Pressure boundary on bottom side using outflow pressure condition.

    if((id ==0) && !lattice->param.GZL && lattice->param.pressure_b_out[subs] )
    {
      k = 0;

      switch(NUM_FLUID_COMPONENTS)
      {
        case 1:
          {  rho_out = lattice->param.rho_out; break; }
        case 2:
          {
            if(subs==0) { rho_out = lattice->param.rho_out;}
            if(subs==1) { rho_out = lattice->param.rho_A[subs];}
            break;
          }
      }

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          ftemp = lattice->pdf[subs][n].ftemp;

          // Bottom, rho_out
          if( lattice->param.incompressible)
          {
            u_z = rho_out
              - ( ftemp[C]
                  + ftemp[N] + ftemp[S] + ftemp[E] + ftemp[W]
                  + 2.*( ftemp[B]
                    + ftemp[BW] + ftemp[BE] + ftemp[BS] + ftemp[BN]));
            c = u_z;
          }
          else // compressible
          {
            u_z =  1.
              - ( ftemp[C]
                  + ftemp[N ] + ftemp[S ] + ftemp[ E] + ftemp[ W]
                  + ftemp[NW] + ftemp[NE] + ftemp[SW] + ftemp[SE]
                  + 2.*( ftemp[B ]
                    + ftemp[BW] + ftemp[BE] + ftemp[BS] + ftemp[BN]))
              / rho_out;
            c = u_z*rho_out;
          }

          ftemp[T ] = ftemp[B ] + (1./3.)*c;
          ftemp[TW] = ftemp[BE] + (1./6.)*c
            + 0.5*(-ftemp[ W] - ftemp[NW] - ftemp[SW]
                + ftemp[ E] + ftemp[NE] + ftemp[SE]);
          ftemp[TE] = ftemp[BW] + (1./6.)*c
            - 0.5*(-ftemp[ W] - ftemp[NW] - ftemp[SW]
                + ftemp[ E] + ftemp[NE] + ftemp[SE]);
          ftemp[TS] = ftemp[BN] + (1./6.)*c
            + 0.5*( ftemp[N ] + ftemp[NW] + ftemp[NE]
                - ftemp[S ] - ftemp[SW] - ftemp[SE]);
          ftemp[TN] = ftemp[BS] + (1./6.)*c
            - 0.5*( ftemp[N ] + ftemp[NW] + ftemp[NE]
                - ftemp[S ] - ftemp[SW] - ftemp[SE]);

#if 1 // Confirm density
          rho = ftemp[ C]
            + ftemp[ E] + ftemp[ W]
            + ftemp[ N] + ftemp[ S]
            + ftemp[ T] + ftemp[ B]
            + ftemp[NE] + ftemp[NW] + ftemp[SE] + ftemp[SW]
            + ftemp[TN] + ftemp[TS] + ftemp[BN] + ftemp[BS]
            + ftemp[TE] + ftemp[TW] + ftemp[BE] + ftemp[BW];
          if( rho_out - rho > 1e-6)
          {
            printf("%s %d >> ERROR: pressure_bottom_out FAIL! "
                "rho_out = %f, rho = %f\n",__FILE__,__LINE__, rho_out, rho);
          }
#endif

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.pressure_b_out[subs] ) */


    //********************************************************************************
    // V E L O C I T Y   T O P   I N F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // velocity top inflow
    //  -- Velocity boundary on top side using inflow velocity condition.
    if( id == get_num_procs( lattice) - 1 && !lattice->param.GZL && 
        lattice->param.velocity_t_in[subs])
    {

      //printf("%s %d >> BOOM!",__FILE__,__LINE__);
      if(subs==0) u = lattice->param.uz_in;    //lattice->param.uz_in; Non-wetting
      if(subs==1) u = lattice->param.uz_in;    //lattice->param.uz_in;

      k = nk-1;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          ftemp = lattice->pdf[subs][n].ftemp;

          if( !lattice->solids[subs][n].is_solid)
          {
            rho = ( ftemp[C ]
                + ftemp[W ] + ftemp[E ] + ftemp[N ] + ftemp[S ]
                + ftemp[NW] + ftemp[NE] + ftemp[SW] + ftemp[SE]
                + 2.*( ftemp[T ]
                  + ftemp[TW] + ftemp[TE] + ftemp[TS] + ftemp[TN])
                )
              / ( 1. + u);
            c = rho*u;

            ftemp[B ] = ftemp[T ] - (1./3.)*c;
            ftemp[BW] = ftemp[TE] - (1./6.)*c + 0.5* (-ftemp[W ]- ftemp[NW]- ftemp[SW]
                + ftemp[E] +ftemp[NE] + ftemp[SE]);
            ftemp[BE] = ftemp[TW] - (1./6.)*c - 0.5* (-ftemp[W ]- ftemp[NW]- ftemp[SW]
                + ftemp[E] +ftemp[NE] + ftemp[SE]);
            ftemp[BS] = ftemp[TN] - (1./6.)*c + 0.5* (ftemp[N ] + ftemp[NW] + ftemp[NE]
                -ftemp[S ] -ftemp[SW] - ftemp[SE]);
            ftemp[BN] = ftemp[TS] - (1./6.)*c - 0.5* (ftemp[N ] + ftemp[NW] + ftemp[NE]
                -ftemp[S ] -ftemp[SW] - ftemp[SE]);

          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_t_in[subs]) */

    // V E L O C I T Y   B O T T O M   O U T F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // velocity bottom inflow
    //  -- Velocity boundary on bottom side using outflow velocity condition.
    if( id == 0 && !lattice->param.GZL && lattice->param.velocity_b_out[subs])
    {
      //printf("%s %d >> BOOM!",__FILE__,__LINE__);
      if(subs==0) u = lattice->param.uz_out;    //lattice->param.uz_out;  Non-wetting
      if(subs==1) u = lattice->param.uz_out;    //lattice->param.uz_out;

      k = 0;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          ftemp = lattice->pdf[subs][n].ftemp;

          if( !lattice->solids[subs][n].is_solid)
          {
            rho = ( ftemp[C ]
                + ftemp[W ] + ftemp[E ] + ftemp[N ] + ftemp[S ]
                + ftemp[NW] + ftemp[NE] + ftemp[SW] + ftemp[SE]
                + 2.*( ftemp[B ]
                  + ftemp[BW] + ftemp[BE] + ftemp[BS] + ftemp[BN])
                )
              / ( 1. - u);
            c = rho*u;


            ftemp[T ] = ftemp[B ] + (1./3.)*c;
            ftemp[TW] = ftemp[BE] + (1./6.)*c + 0.5* (-ftemp[W ]- ftemp[NW]- ftemp[SW]
                + ftemp[E] +ftemp[NE] + ftemp[SE]);
            ftemp[TE] = ftemp[BW] + (1./6.)*c - 0.5* (-ftemp[W ]- ftemp[NW]- ftemp[SW]
                + ftemp[E] +ftemp[NE] + ftemp[SE]);
            ftemp[TS] = ftemp[BN] + (1./6.)*c  + 0.5* (ftemp[N ] + ftemp[NW] + ftemp[NE]
                -ftemp[S ] -ftemp[SW] - ftemp[SE]);
            ftemp[TN] = ftemp[BS] + (1./6.)*c  - 0.5* (ftemp[N ] + ftemp[NW] + ftemp[NE]
                -ftemp[S ] -ftemp[SW] - ftemp[SE]);


          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_b_in[subs]) */

    //********************************************************************************


    // GZL P R E S S U R E // V E L O C I T Y   T O P   I N F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //  -- Velocity boundary on top side using inflow velocity condition.
#if PARALLEL
    id = get_proc_id(lattice);
#else
    id = get_num_procs(lattice)-1;
#endif

    if( (id ==get_num_procs(lattice)-1) && lattice->param.GZL && lattice->param.PressureBC)//lattice->param.velocity_t_in[subs])
    {
      printf("%s %d >> BOOM!",__FILE__,__LINE__);
      k = nk-1;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          f = lattice->pdf[subs][n].f;

#if STORE_UEQ
          //       u_x = &lattice->ueq[n].u[0] = 0.02;
          //       lattice->ueq[n].u[1] = 0.;
          //       lattice->ueq[n].u[2] = 0.;
#else /* !( STORE_UEQ) */
          //       lattice->macro_vars[subs][n].u[0] =  0.02 ;
          //       lattice->macro_vars[subs][n].u[1] =  0. ;
          //       lattice->macro_vars[subs][n].u[2] =  0. ;
#endif /* STORE_UEQ */

          n1 = IJK2N( i, j, k-1, ni, nj);
          f1 = lattice->pdf[subs][n1].f;
          feq1 = lattice->pdf[subs][n1].feq;

          //  rhoo = &( lattice->macro_vars[subs][n].rho);
          switch(NUM_FLUID_COMPONENTS)
          {case 1:
            {
              rho =  /*lattice->macro_vars[subs][n1].rho;*/lattice->param.rho_in;
              u_x =  0.;//lattice->macro_vars[subs][n1].u[0] ;
              u_y =  0.;//lattice->macro_vars[subs][n1].u[1] ;
              u_z =  lattice->macro_vars[subs][n1].u[2] ;
              lattice->macro_vars[subs][n].u[2] = u_z;
              break;
            }
            case 2:
            {
#if STORE_UEQ
              if(subs==0)
              {
                rho =lattice->param.rho_in;// 8+ 0.001*(real)lattice->time;// ;    subs= 0 Non-wetting
                if( lattice->time >1000)  {rho = lattice->param.rho_in;}
              }
              if(subs==1) rho = 0.000 ;    //lattice->param.uz_in; subs= 1 Wetting

              //          rho =  lattice->param.rho_A[subs];
              u_x =  0.;// lattice->ueq[n1].u[0] ;
              u_y =  0.;// lattice->ueq[n1].u[1] ;
              u_z =  lattice->ueq[n1].u[2] ;
              lattice->macro_vars[subs][n].u[2] = u_z;
              lattice->ueq[n].u[2] = u_z;

#endif
              break;
            }
          }
          lattice->macro_vars[subs][n].rho = rho;
          lattice->macro_vars[subs][n].u[0] = u_x;
          lattice->macro_vars[subs][n].u[1] = u_y;


          if( !lattice->solids[subs][n].is_solid)
          {
            compute_a_feq( feq, rho, u_x, u_y, u_z );

            for (a= 0; a<Q; a++)
            {
              temp[a] =  f1[a] - feq1[a];
              f[a]=feq[a]+ temp[a];
            }


          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_GZL) */

    //GZL P R E S S U R E // V E L O C I T Y   B O T T O M   O U T F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // velocity bottom inflow
    //  -- Velocity boundary on bottom side using inflow velocity condition.
#if PARALLEL
    id = get_proc_id(lattice);
#else
    id = 0;
#endif
    if( (id ==0) && lattice->param.GZL && lattice->param.PressureBC)//lattice->param.velocity_b_in[subs])
    {
      printf("%s %d >> BOOM!",__FILE__,__LINE__);
      k = 0;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          f = lattice->pdf[subs][n].f;

          n1 = IJK2N( i, j, k+1, ni, nj);
          f1   = lattice->pdf[subs][n1].f;
          feq1 = lattice->pdf[subs][n1].feq;

          switch(NUM_FLUID_COMPONENTS)
          {case 1:
            {
              rho =  /*lattice->macro_vars[subs][n1].rho;*/lattice->param.rho_out;
              u_x =  0.;//lattice->macro_vars[subs][n1].u[0] ;
              u_y =  0.;//lattice->macro_vars[subs][n1].u[1] ;
              u_z =  lattice->macro_vars[subs][n1].u[2] ;
              lattice->macro_vars[subs][n].u[2] = u_z;
              break;
            }
            case 2:
            {
#if STORE_UEQ
              if(subs==0) rho = 0.000;
              if(subs==1)
              {
                rho =lattice->param.rho_out; //8.0- 0.001*(real)lattice->time;//
                if( lattice->time >6000)   {rho = lattice->param.rho_out;}
              }
              //          rho =  lattice->param.rho_B[subs];
              u_x =  0.;//lattice->ueq[n1].u[0] ;
              u_y =  0.;//lattice->ueq[n1].u[1] ;
              u_z =  lattice->ueq[n1].u[2] ;
              lattice->macro_vars[subs][n].u[2] = u_z;
              lattice->ueq[n].u[2] = u_z;
#endif
              break;
            }
          }
          lattice->macro_vars[subs][n].rho = rho;
          lattice->macro_vars[subs][n].u[0] = u_x;
          lattice->macro_vars[subs][n].u[1] = u_y;



          if( !lattice->solids[subs][n].is_solid)
          {
            compute_a_feq( feq, rho, u_x, u_y, u_z );

            for (a= 0; a<Q; a++)
            {
              temp[a] =f1[a] - feq1[a];
              f[a] =  feq[a]+ temp[a];
            }

          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_GZL) */



    //PAN   V E L O C I T Y   T O P   I N F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //  -- Velocity boundary on top side using inflow velocity condition.
    if( lattice->param.GZL && (!lattice->param.PressureBC))//lattice->param.velocity_t_in[subs])
    {
      printf("%s %d >> BOOM!",__FILE__,__LINE__);
      k = nk-1;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          f = lattice->pdf[subs][n].f;
          //  feq =  lattice->pdf[subs][n].feq;


          n1 = IJK2N( i, j, k-1, ni, nj);
          f1 = lattice->pdf[subs][n1].f;
          feq1 = lattice->pdf[subs][n1].feq;

          rhoo = &( lattice->macro_vars[subs][n].rho);
          switch(NUM_FLUID_COMPONENTS)
          {case 1:
            {
              rho =  /*lattice->macro_vars[subs][n1].rho;*/lattice->param.rho_in;
              u_x =  0.;//lattice->macro_vars[subs][n1].u[0] ;
              u_y =  0.;//lattice->macro_vars[subs][n1].u[1] ;
              u_z =  -0.0;//lattice->macro_vars[subs][n1].u[2] ;
              break;
            }
            case 2:
            {
#if STORE_UEQ
              rho    =  lattice->param.rho_A[subs];
              u_x = 0.;// lattice->ueq[n1].u[0] ;
              u_y = 0.;// lattice->ueq[n1].u[1] ;
              u_z = 0.;// lattice->ueq[n1].u[2] ;
#endif
              break;
            }
          }
          lattice->macro_vars[subs][n].rho = rho;
          lattice->macro_vars[subs][n].u[0] = u_x;
          lattice->macro_vars[subs][n].u[1] = u_y;
          lattice->macro_vars[subs][n].u[2] = u_z;


          if( !lattice->solids[subs][n].is_solid)
          {
            compute_a_feq( feq, rho, u_x, u_y, u_z );

            for (a= 0; a<Q; a++)
            {f[a] = feq[a];}//+(f1[a] - feq1[a]);


          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_GZL) */

    //PAN V E L O C I T Y   B O T T O M   I N F L O W   B C
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // velocity bottom inflow
    //  -- Velocity boundary on bottom side using inflow velocity condition.
    if(  lattice->param.GZL && (!lattice->param.PressureBC))//lattice->param.velocity_b_in[subs])
    {
      printf("%s %d >> BOOM!",__FILE__,__LINE__);
      k = 0;

      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          f = lattice->pdf[subs][n].f;
          //  feq =  lattice->pdf[subs][n].feq;

          n1 = IJK2N( i, j, k+1, ni, nj);
          f1 = lattice->pdf[subs][n1].f;
          feq1 = lattice->pdf[subs][n1].feq;


          switch(NUM_FLUID_COMPONENTS)
          {case 1:
            {
              rho =  /*lattice->macro_vars[subs][n1].rho;*/lattice->param.rho_out;
              u_x =  0.;//lattice->macro_vars[subs][n1].u[0] ;
              u_y =  0.;//lattice->macro_vars[subs][n1].u[1] ;
              u_z =  -0.0;//lattice->macro_vars[subs][n1].u[2] ;
              break;
            }
            case 2:
            {
#if STORE_UEQ
              rho =  lattice->param.rho_B[subs];
              u_x =  0.;//lattice->ueq[n1].u[0] ;
              u_y =  0.;//lattice->ueq[n1].u[1] ;
              u_z =  0.;//lattice->ueq[n1].u[2] ;
#endif
              break;
            }
          }

          lattice->macro_vars[subs][n].rho = rho;
          lattice->macro_vars[subs][n].u[0] = u_x;
          lattice->macro_vars[subs][n].u[1] = u_y;
          lattice->macro_vars[subs][n].u[2] = u_z;

          if( !lattice->solids[subs][n].is_solid)
          {
            compute_a_feq( feq, rho, u_x, u_y, u_z );
            if(subs==1)
              compute_a_feq( feq1, 2.0- 0.0003*(real)lattice->time, u_x, u_y, u_z );

            for (a= 0; a<Q; a++)
            {
              f[a] = feq[a];

            }// +(f1[a] - feq1[a]);

          } /* if( !is_solid) */

        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */

    } /* if( lattice->param.velocity_GZL) */


  } /* for( subs=0; subs<(NUM_FLUID_COMPONENTS)-(INAMURO_SIGMA_COMPONENT); subs++) */


#endif
}



#endif

void compute_a_feq( real *feq, real rho, real u_x, real u_y, real u_z)
{
#if 0
  int  a;

  real tau;

  real W0,   W1,   W2;
  real ux,   uy,  uz,
       uxsq, uysq, uzsq, usq;
  real c;
  real udotx;

  //  real rho, u_x, u_y, u_z;

  W0 = 1./3.;
  W1 = 1./18.;
  W2 = 1./36.;

  tau       =    1.0;

  usq = u_x*u_x + u_y*u_y + u_z*u_z;

  feq[0] = W0 * rho*(1. - 1.5*usq);

  for( a=1; a<=6; a++)
  {
    udotx = ((real)vx[a]*u_x+(real)vy[a]*u_y+(real)vz[a]*u_z);
    feq[a] = W1*rho*(1. + 3.*udotx + 4.5 *udotx*udotx - 1.5*usq);
  }

  for( a=7; a<Q; a++)
  {
    udotx = ((real)vx[a]*u_x+(real)vy[a]*u_y+(real)vz[a]*u_z);
    feq[a]  = W2*rho*(1. + 3.*udotx + 4.5*udotx*udotx - 1.5*usq);
  }


#else
  // TODO
#endif
}




