//##############################################################################
//
// compute.c
//
//  - Routines for computing on the lattice:
//
//    - compute_rho_and_u
//    - compute_feq
//    - compute_big_u
//    - compute_gforce
//    - compute_fluid_fluid_force
//    - etc...
//

// void compute_macro_vars( struct lattice_struct *lattice)
//##############################################################################
//
// C O M P U T E   M A C R O   V A R S
//
//  - Compute macroscopic variables.

#if INAMURO_SIGMA_COMPONENT
// TODO:
#else /* !( INAMURO_SIGMA_COMPONENT) */
void compute_macro_vars( struct lattice_struct *lattice)
{
#if 0
  int a, n, k, i,j, ni, nj, nk;

  real *rho[ NUM_FLUID_COMPONENTS],
         *u_x[ NUM_FLUID_COMPONENTS],
         *u_y[ NUM_FLUID_COMPONENTS],
         *u_z[ NUM_FLUID_COMPONENTS];

  real ux_sum, uy_sum, uz_sum;

  real *ueq;

  real *ftemp;

  int    is_solid;

  int    subs;

  real tau0,
         tau1;

  int mpierr;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  for(subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {

    for( n=0; n<lattice->NumNodes; n++)
    {
      rho[subs] = &( lattice->macro_vars[subs][n].rho);
      u_x[subs] = &( lattice->macro_vars[subs][n].u[0]);
      u_y[subs] = &( lattice->macro_vars[subs][n].u[1]);
      u_z[subs] = &( lattice->macro_vars[subs][n].u[2]);

      ftemp     =   lattice->pdf[subs][n].ftemp;
      is_solid  = ( lattice->solids[subs][n].is_solid);

      if(!(is_solid))
      {
        *rho[subs] = 0.;
        *u_x[subs] = 0.;
        *u_y[subs] = 0.;
        *u_z[subs] = 0.;

        for(a=0; a<Q; a++)
        {
          (*rho[subs]) +=       (ftemp[a]);
          (*u_x[subs]) += vx[a]*(ftemp[a]);
          (*u_y[subs]) += vy[a]*(ftemp[a]);
          (*u_z[subs]) += vz[a]*(ftemp[a]);
        }

      } /*if(!(is_solid))*/
    } /*for( n=0; n<lattice->NumNodes; n++) */
  } /* for(subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

  if( NUM_FLUID_COMPONENTS==2)
  {
    tau0 = lattice->param.tau[0];
    tau1 = lattice->param.tau[1];

    for( n=0; n<lattice->NumNodes; n++)
    {
      for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
      {
        rho[subs]   = &( lattice->macro_vars[subs][n].rho);
        u_x[subs]   = &( lattice->macro_vars[subs][n].u[0]);
        u_y[subs]   = &( lattice->macro_vars[subs][n].u[1]);
        u_z[subs]   = &( lattice->macro_vars[subs][n].u[2]);
      }

#if STORE_UEQ
      ueq = lattice->ueq[n].u;
#endif /* STORE_UEQ */

      is_solid  = ( lattice->solids[0][n].is_solid);

      if( !( is_solid))
      {
        ux_sum =  *u_x[0]/tau0 + *u_x[1]/tau1;
        uy_sum =  *u_y[0]/tau0 + *u_y[1]/tau1;
        uz_sum =  *u_z[0]/tau0 + *u_z[1]/tau1;

#if STORE_UEQ
        if( *rho[0] + *rho[1] != 0.)
        {
          //        ueq[0] = ( ux_sum) / ( *rho[0]/tau0 + *rho[1]/tau1);
          //        ueq[1] = ( uy_sum) / ( *rho[0]/tau0 + *rho[1]/tau1);
          //        ueq[2] = ( uz_sum) / ( *rho[0]/tau0 + *rho[1]/tau1);
          *ueq++ = ( ux_sum) / ( *rho[0]/tau0 + *rho[1]/tau1);
          *ueq++ = ( uy_sum) / ( *rho[0]/tau0 + *rho[1]/tau1);
          *ueq++ = ( uz_sum) / ( *rho[0]/tau0 + *rho[1]/tau1);
        }
        else
        {
          *ueq++ = 0.;
          *ueq++ = 0.;
          *ueq++ = 0.;
        }
#endif /* STORE_UEQ */


        if( ux_sum != 0.)
        {
          if( *rho[0] != 0.) { *u_x[0] = *u_x[0] / *rho[0]; }
          else {             *u_x[0] = 0.; }
          if( *rho[1] != 0.) { *u_x[1] = *u_x[1] / *rho[1]; }
          else {             *u_x[1] = 0.; }
        }
        else { *u_x[0] = 0.; *u_x[1] = 0.; }

        if( uy_sum != 0.)
        {
          if( *rho[0] != 0.) { *u_y[0] = *u_y[0] / *rho[0]; }
          else {             *u_y[0] = 0.; }
          if( *rho[1] != 0.) { *u_y[1] = *u_y[1] / *rho[1]; }
          else {             *u_y[1] = 0.; }
        }
        else { *u_y[0] = 0.; *u_y[1] = 0.; }

        if( uz_sum != 0.)
        {
          if( *rho[0] != 0.) { *u_z[0] = *u_z[0] / *rho[0]; }
          else {             *u_z[0] = 0.; }
          if( *rho[1] != 0.) { *u_z[1] = *u_z[1] / *rho[1]; }
          else {             *u_z[1] = 0.; }
        }
        else { *u_z[0] = 0.; *u_z[1] = 0.; }

      } /* if( !( is_solid)) */

    } /* for( n=0; n<lattice->NumNodes; n++) */

  } /* if( NUM_FLUID_COMPONENTS==2) */

  else if( NUM_FLUID_COMPONENTS == 1)
  {
    for( n=0; n<lattice->NumNodes; n++)
    {
      rho[0]      = &( lattice->macro_vars[0][n].rho);
      u_x[0]      = &( lattice->macro_vars[0][n].u[0]);
      u_y[0]      = &( lattice->macro_vars[0][n].u[1]);
      u_z[0]      = &( lattice->macro_vars[0][n].u[2]);
      is_solid    =  ( lattice->solids[0][n].is_solid);


      if( !( is_solid) )
      {
        if( *rho[0] != 0. && *u_x[0] != 0.)
        {
          *u_x[0] = *u_x[0] / *rho[0];
        }
        else
        {
          *u_x[0] = 0.;
        }

        if( *rho[0] != 0. && *u_y[0] != 0.)
        {
          *u_y[0] = *u_y[0] / *rho[0];
        }
        else
        {
          *u_y[0] = 0.;
        }

        if( *rho[0] != 0. && *u_z[0] != 0.)
        {
          *u_z[0] = *u_z[0] / *rho[0];
        }
        else
        {
          *u_z[0] = 0.;
        }

      } /* if( !( is_solid)) */

    } /* for( n=0; n<lattice->NumNodes; n++) */
  }
  else
  {
    printf(
        "compute_macro_vars() -- "
        "Unhandled case "
        "NUM_FLUID_COMPONENTS = %d . "
        "Exiting!\n",
        NUM_FLUID_COMPONENTS);
    process_exit(1);
  }
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if PARALLEL

  rho_send_recv_begin(  lattice, 0);		//these calls deal with BOTH subs 0 and 1
  rho_send_recv_end(  lattice, 0);

  solid_send_recv_begin(  lattice, 0);		//these calls deal with subs 0, which is all that is required
  solid_send_recv_end(  lattice, 0);

#endif
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if NON_LOCAL_FORCES
  switch( NUM_FLUID_COMPONENTS)
  {
    case 1:
      {
        compute_phase_force( lattice, 0);
        compute_single_fluid_solid_force( lattice, 0);
        break;
      }

    case 2:
      {
        compute_fluid_fluid_force( lattice);
        compute_real_fluid_solid_force( lattice);
        break;
      }

    default:
      {
        printf("%s %d >> ERROR: Unhandled case %d\n",
            __FILE__,__LINE__,
            NUM_FLUID_COMPONENTS);
        process_exit(1);
        break;
      }
  }
#endif /* NON_LOCAL_FORCES */

#else
// TODO
#endif
} /* void compute_macro_vars( struct lattice_struct *lattice) */
#endif /* INAMURO_SIGMA_COMPONENT */

// void compute_feq( struct lattice_struct *lattice)
//##############################################################################
//
// C O M P U T E   F E Q
//
//  - Compute equilibrium distribution function, feq.
//
#if SAVE_MEMO

#else
void compute_feq( struct lattice_struct *lattice)
{
#if 0
  int n, a;

  real *feq;
  int    subs;

 for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {

    for( n=0; n<lattice->NumNodes; n++)
    {
      feq       =    lattice->pdf[subs][n].feq;
      compute_single_feq(  lattice, n, subs, feq);

    }/*for( n=0; n<lattice->NumNodes; n++)*/
  }/*for( subs=0; subs<(NUM_FLUID_COMPONENTS)-(INAMURO_SIGMA_COMPONENT); ... */

#else
// TODO
#endif
} /* void compute_feq( struct lattice_struct *lattice) */
#endif

void compute_single_feq( struct lattice_struct *lattice, int n, int subs, real *feq)
{
#if 0
  int  a;

  real tau;

  real W0,   W1,   W2;
  real ux,   uy,  uz,
         uxsq, uysq, uzsq, usq;
  real c; /*sound velocity*/
  real udotx;

  real rho, u_x, u_y, u_z;


#if STORE_UEQ
  real *ueq;
#endif /* STORE_UEQ */
  int    is_solid;



  /*Refers to Klaus Lglburger Thesis*/
  W0 = 1./3.;
  W1 = 1./18.;
  W2 = 1./36.;

      tau       =    lattice->param.tau[subs];
      rho       =    lattice->macro_vars[subs][n].rho;
      is_solid  =    lattice->solids[subs][n].is_solid;

//      if( /*fcountone*/  1 || !(is_solid))
      if( !lattice->time || !(is_solid))
      {
#if STORE_UEQ
        // Start with the composite macroscopic velocities.
        u_x       =    lattice->ueq[n].u[0];
        u_y       =    lattice->ueq[n].u[1];
        u_z       =    lattice->ueq[n].u[2];
        //  u_x       = *ueq;     *ueq++;
        //  u_y       = *ueq;     *ueq++;
        //  u_z       = *ueq;     *ueq++;
#else /* !( STORE_UEQ) */
        // Start with the individual component velocities.
        u_x       =    lattice->macro_vars[subs][n].u[0];
        u_y       =    lattice->macro_vars[subs][n].u[1];
        u_z       =    lattice->macro_vars[subs][n].u[2];
#endif /* STORE_UEQ */

        u_x       = u_x
#if NON_LOCAL_FORCES
          +tau*lattice->force[ subs][n].force[0]
          +tau*lattice->force[ subs][n].sforce[0]
#endif /* NON_LOCAL_FORCES */
          +tau*lattice->param.gforce[ subs][0]
          ;
        u_y       = u_y
#if NON_LOCAL_FORCES
          +tau*lattice->force[ subs][n].force[1]
          +tau*lattice->force[ subs][n].sforce[1]
#endif /* NON_LOCAL_FORCES */
          +tau*lattice->param.gforce[ subs][1]
          ;
        u_z       = u_z
#if NON_LOCAL_FORCES
          +tau*lattice->force[ subs][n].force[2]
          +tau*lattice->force[ subs][n].sforce[2]
#endif /* NON_LOCAL_FORCES */
          + tau*lattice->param.gforce[ subs][2]
          ;

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

      }/*if(!is_solid)*/
      else
      {
        for( a=0; a<Q; a++)
        {
          feq[a] = 0.;
        }
      }

#else
// TODO
#endif
} /* void compute_single_feq( struct lattice_struct *lattice) */


#if NON_LOCAL_FORCES
void compute_phase_force( lattice_ptr lattice, int subs)
{
#if 0
  real ***psi; //psi[LZ][LY][LX];
  real psi_temp;

  real rho;

  real *force;

  int    a;

  int    i,  j,  k,
         in, jn, kn,
         ip, jp, kp,
         ia, ja, ka;

  int    n,
         ni, nj, nk;

//printf("%s %d >> compute_phase_force() -- Hi!\n",
//    __FILE__,__LINE__);

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  psi =    ( real***)malloc( nk*sizeof(real**));
  for( k=0; k<nk; k++)
  {
    psi[k] = ( real**)malloc( nj*sizeof(real*));
  }
  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      psi[k][j] = ( real*)malloc( ni*sizeof(real));
    }
  }

  // Initialize psi.
  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        psi[k][j][i] = 0.;
      }
    }
  }

  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        n = IJK2N( i, j, k, ni, nj);
        rho = lattice->macro_vars[subs][n].rho;

        if( rho!=0 && !( lattice->solids[subs][n].is_solid))
        {
          psi[k][j][i] = 4.*exp(-200./(rho));
        }
        else
        {
          psi[k][j][i] = 0.;
        }
      } /* for( i=0; i<ni; i++) */
    } /* for( j=0; j<nj; j++) */
  } /* for( k=0; k<nk; k++) */

  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        n = IJK2N( i, j, k, ni, nj);
        force = lattice->force[subs][n].force;

        force[0] = 0.;
        force[1] = 0.;
        force[2] = 0.;

        if( !( lattice->solids[subs][ n].is_solid ))
        {
          for( a=1; a<=6; a++)
          {
            ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
            ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
            ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
            if( !( lattice->solids[subs][
                  IJK2N( ia, ja, ka, ni, nj)].is_solid ))
            {
              force[0] += WM*vx[a]*psi[ ka][ ja][ ia];
              force[1] += WM*vy[a]*psi[ ka][ ja][ ia];
              force[2] += WM*vz[a]*psi[ ka][ ja][ ia];
            }
          }

          for( a=7; a<Q; a++)
          {
            ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
            ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
            ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
            if( !( lattice->solids[subs][
                  IJK2N( ia, ja, ka, ni, nj)].is_solid ))
            {
              force[0] += WD*vx[a]*psi[ ka][ ja][ ia];
              force[1] += WD*vy[a]*psi[ ka][ ja][ ia];
              force[2] += WD*vz[a]*psi[ ka][ ja][ ia];
            }
          }

          force[0] = -lattice->param.big_V0*psi[k][j][i]*( force[0]);
          force[1] = -lattice->param.big_V0*psi[k][j][i]*( force[1]);
          force[2] = -lattice->param.big_V0*psi[k][j][i]*( force[2]);

        } /* if( !( lattice->solids[subs][ j*ni+i].is_solid )) */

        else
        {
            *( force  ) = 0.;
            *( force+1) = 0.;
            *( force+2) = 0.;
        }

      } /* for( i=0; i<ni; i++) */
    } /* for( j=0; j<nj; j++) */
  } /* for( k=0; k<nk; k++) */

 for( k=0; k<nk; k++)
 {
   for( j=0; j<nj; j++)
   {
     free( psi[k][j]);
   }
 }
 for( k=0; k<nk; k++)
 {
   free( psi[k]);
 }
 free( psi);

//printf("%s %d >> compute_phase_force() -- Bi!\n",
//    __FILE__,__LINE__);

#else
// TODO
#endif
} /* void compute_phase_force( lattice_ptr lattice) */

void compute_fluid_fluid_force( lattice_ptr lattice)
{
#if 0
  real psi_temp;

  real rho, rho0, rho1;

  real *force[2];

  real sumx[2], sumy[2], sumz[2];

  int    a, subs;
  int    id;

  int    i,  j,  k, k1, k2,
         in, jn, kn,
         ip, jp, kp,
         ia, ja, ka;

  int    n, na,
         ni, nj, nk;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

#if PARALLEL
        id = get_proc_id(lattice);
//  printf("%d\n",id);

  { k1 =0;   k2 = nk;}
  if(id ==0 && (!lattice->param.AllBoundaryPeriodic))
  { k1 =3;   k2 = nk;}
  if((id ==get_num_procs(lattice)-1) && (!lattice->param.AllBoundaryPeriodic))
  { k1 =0;   k2 = nk-2;}
//  printf("k1=%d, k2=%d, #########\n", k1, k2);
#else
  if(lattice->param.AllBoundaryPeriodic)
  { k1 =0;   k2 = nk;}
  else
  { k1 =3;   k2 = nk-2;}

#endif


  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
// for fluid 1 injection into porous media filled with fluid 2, rev_Huang
#if PARALLEL
//     rho_send_recv_begin( lattice, subs);
#endif
    for( k=k1; k<k2; k++)
    {
      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);
          force[subs] = lattice->force[subs][n].force;

          force[subs][0] = 0.;
          force[subs][1] = 0.;
          force[subs][2] = 0.;

          if( !( lattice->solids[subs][ n].is_solid ))
          {
            for( a=1; a<=6; a++)
            {
              ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
              ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );

#if PARALLEL
        if( k+vz[a]<nk &&  k+vz[a]>=0 )
        {
        ka = k+vz[a];
        na = IJK2N( ia, ja, ka, ni, nj);
              if( !( lattice->solids[subs][ na].is_solid ))
               {
                rho = lattice->macro_vars[subs][na].rho;
                force[subs][0] += WM*vx[a]*rho;
                force[subs][1] += WM*vy[a]*rho;
                force[subs][2] += WM*vz[a]*rho;
               }

        }
        else if(k+vz[a]>=nk)
        {
              if( !( lattice->process.z_neg_solid_to_recv[ja*ni+ia] ))
               {
                rho = lattice->process.z_neg_rho_to_recv[subs*ni*nj +ja*ni+ia];
                force[subs][0] += WM*vx[a]*rho;
                force[subs][1] += WM*vy[a]*rho;
                force[subs][2] += WM*vz[a]*rho;
// lattice->process.z_neg_solid_to_recv[n];
         }

        }
        else
        {
              if( !( lattice->process.z_pos_solid_to_recv[ja*ni+ia] ))
               {
                rho = lattice->process.z_pos_rho_to_recv[subs*ni*nj +ja*ni+ia];
                force[subs][0] += WM*vx[a]*rho;
                force[subs][1] += WM*vy[a]*rho;
                force[subs][2] += WM*vz[a]*rho;
         }
        }

#else
        ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
        na = IJK2N( ia, ja, ka, ni, nj);
              if( !( lattice->solids[subs][ na].is_solid ))
              {
                rho = lattice->macro_vars[subs][na].rho;
                force[subs][0] += WM*vx[a]*rho;
                force[subs][1] += WM*vy[a]*rho;
                force[subs][2] += WM*vz[a]*rho;
              }

#endif
      }//a=0,... 6

            for( a=7; a<Q; a++)
            {
              ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
              ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
#if PARALLEL
        if( k+vz[a]<nk &&  k+vz[a]>=0 )
        {
        ka = k+vz[a];
        na = IJK2N( ia, ja, ka, ni, nj);
              if( !( lattice->solids[subs][ na].is_solid ))
               {
                rho = lattice->macro_vars[subs][na].rho;
                force[subs][0] += WD*vx[a]*rho;
                force[subs][1] += WD*vy[a]*rho;
                force[subs][2] += WD*vz[a]*rho;
               }

        }
        else if(k+vz[a]>=nk)
        {
              if( !( lattice->process.z_neg_solid_to_recv[ja*ni+ia] ))
               {
                rho = lattice->process.z_neg_rho_to_recv[subs*ni*nj +ja*ni+ia];
                force[subs][0] += WD*vx[a]*rho;
                force[subs][1] += WD*vy[a]*rho;
                force[subs][2] += WD*vz[a]*rho;
// lattice->process.z_neg_solid_to_recv[n];
         }

        }
        else
        {
              if( !( lattice->process.z_pos_solid_to_recv[ja*ni+ia] ))
               {
                rho = lattice->process.z_pos_rho_to_recv[subs*ni*nj +ja*ni+ia];
                force[subs][0] += WD*vx[a]*rho;
                force[subs][1] += WD*vy[a]*rho;
                force[subs][2] += WD*vz[a]*rho;
         }
        }


#else
              ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
              na = IJK2N( ia, ja, ka, ni, nj);
              if( !( lattice->solids[subs][ na].is_solid ))
              {
                rho = lattice->macro_vars[subs][na].rho;
                force[subs][0] += WD*vx[a]*rho;
                force[subs][1] += WD*vy[a]*rho;
                force[subs][2] += WD*vz[a]*rho;
              }
#endif

            }//a=7....18

          } /* if( !( lattice->solids[subs][ n].is_solid )) */
        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */
    } /* for( k=0; k<nk; k++) */

#if PARALLEL
//    rho_send_recv_end( lattice, subs);
#endif
  } /* for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++) */

    for( k=0; k<nk; k++)
    {
      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          n = IJK2N( i, j, k, ni, nj);

          if( !( lattice->solids[0][ n].is_solid ))
          {
            force[0] = lattice->force[0][n].force;
            force[1] = lattice->force[1][n].force;
            rho0  = lattice->macro_vars[0][n].rho;
            rho1  = lattice->macro_vars[1][n].rho;

            sumx[/*subs*/0] = force[/*subs*/0][/*x*/0];
            sumx[/*subs*/1] = force[/*subs*/1][/*x*/0];

//            force[0][0] = -lattice->param.big_V0*rho0*sumx[1];
//            force[1][0] = -lattice->param.big_V0*rho1*sumx[0];
            force[0][0] = -lattice->param.big_V0*sumx[1];
            force[1][0] = -lattice->param.big_V0*sumx[0];

            sumy[/*subs*/0] = force[/*subs*/0][/*y*/1];
            sumy[/*subs*/1] = force[/*subs*/1][/*y*/1];

//            force[0][1] = -lattice->param.big_V0*rho0*sumy[1];
//            force[1][1] = -lattice->param.big_V0*rho1*sumy[0];
            force[0][1] = -lattice->param.big_V0*sumy[1];
            force[1][1] = -lattice->param.big_V0*sumy[0];

            sumz[/*subs*/0] = force[/*subs*/0][/*z*/2];
            sumz[/*subs*/1] = force[/*subs*/1][/*z*/2];

//            force[0][2] = -lattice->param.big_V0*rho0*sumz[1];
//            force[1][2] = -lattice->param.big_V0*rho1*sumz[0];
            force[0][2] = -lattice->param.big_V0*sumz[1];
            force[1][2] = -lattice->param.big_V0*sumz[0];

#if 0
            if( force[0][0]+ force[1][0]+ force[0][1]+ force[1][1]+ force[0][2]+ force[1][2] != 0.)
            {
printf("%s %d >> %f %f, %f %f, %f %f\n",
  __FILE__,__LINE__,
    force[0][0],
    force[1][0],
    force[0][1],
    force[1][1],
    force[0][2],
    force[1][2]
    );
            }
#endif

          } /* if( !( lattice->solids[subs][ n].is_solid )) */
        } /* for( i=0; i<ni; i++) */
      } /* for( j=0; j<nj; j++) */
    } /* for( k=0; k<nk; k++) */


#else
// TODO
#endif
} /* void compute_fluid_fluid_force( lattice_ptr lattice) */


void compute_single_fluid_solid_force( lattice_ptr lattice, int subs)
{
#if 0
  real ***psi; //psi[ NUM_FLUID_COMPONENTS][LX][LY];
  real psi_temp;

  real rho;

  real *sforce;

  int    a;

  int    i,  j,  k,
         in, jn, kn,
         ip, jp, kp,
         ia, ja, ka;

  int    n,
         ni, nj, nk;

//printf("%s %d >> compute_single_fluid_solid_force() -- Hi!\n",
//    __FILE__,__LINE__);

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;

  psi =    ( real***)malloc( nk*sizeof(real**));
  for( k=0; k<nk; k++)
  {
    psi[k] = ( real**)malloc( nj*sizeof(real*));
  }
  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      psi[k][j] = ( real*)malloc( ni*sizeof(real));
    }
  }

  // Initialize psi.
  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        psi[k][j][i] = 0.;
      }
    }
  }

  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        n = IJK2N( i, j, k, ni, nj);
        rho = lattice->macro_vars[subs][n].rho;

        if( rho!=0 && !( lattice->solids[subs][n].is_solid))
        {
          psi[k][j][i] = 4.*exp(-200./(rho));
        }
        else
        {
          psi[k][j][i] = 0.;
        }
      } /* for( i=0; i<ni; i++) */
    } /* for( j=0; j<nj; j++) */
  } /* for( k=0; k<nk; k++) */

  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        n = IJK2N( i, j, k, ni, nj);
        sforce = lattice->force[subs][n].sforce;

        sforce[0] = 0.;
        sforce[1] = 0.;
        sforce[2] = 0.;

        if( !( lattice->solids[subs][ n].is_solid ))
        {
          for( a=1; a<=6; a++)
          {
            ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
            ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
            ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
            if(  ( lattice->solids[subs][
                  IJK2N( ia, ja, ka, ni, nj)].is_solid ))
            {
              sforce[0] += WM*vx[a];
              sforce[1] += WM*vy[a];
              sforce[2] += WM*vz[a];
            }
          }

          for( a=7; a<Q; a++)
          {
            ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
            ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
            ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
            if(  ( lattice->solids[subs][
                  IJK2N( ia, ja, ka, ni, nj)].is_solid ))
            {
              sforce[0] += WD*vx[a];
              sforce[1] += WD*vy[a];
              sforce[2] += WD*vz[a];
            }
          }

          sforce[0] =
            -lattice->param.big_V0_solid[subs]*psi[k][j][i]*( sforce[0]);
          sforce[1] =
            -lattice->param.big_V0_solid[subs]*psi[k][j][i]*( sforce[1]);
          sforce[2] =
            -lattice->param.big_V0_solid[subs]*psi[k][j][i]*( sforce[2]);

        } /* if( !( lattice->solids[subs][ j*ni+i].is_solid )) */

        else
        {
            *( sforce  ) = 0.;
            *( sforce+1) = 0.;
            *( sforce+2) = 0.;
        }

      } /* for( i=0; i<ni; i++) */
    } /* for( j=0; j<nj; j++) */
  } /* for( k=0; k<nk; k++) */

 for( k=0; k<nk; k++)
 {
   for( j=0; j<nj; j++)
   {
     free( psi[k][j]);
   }
 }
 for( k=0; k<nk; k++)
 {
   free( psi[k]);
 }
 free( psi);

//printf("%s %d >> compute_single_fluid_solid_force() -- Bi!\n",
//    __FILE__,__LINE__);

#else
// TODO
#endif
} /* void compute_single_fluid_solid_sforce( lattice_ptr lattice, int subs) */

void compute_real_fluid_solid_force( lattice_ptr lattice)
{
#if 0
  real rho;

  real *sforce[ /*NUM_FLUID_COMPONENTS*/ 2];

  real sum_x, sum_y, sum_z;

  int    a;

  int    i,  j,  k,
         in, jn, kn,
         ip, jp, kp,
         ia, ja, ka;

  int    n,
         ni, nj, nk, na;

  ni = lattice->param.LX;
  nj = lattice->param.LY;
  nk = lattice->param.LZ;
#if PARALLEL
//     solid_send_recv_begin( lattice, 0);
#endif
  for( k=0; k<nk; k++)
  {
    for( j=0; j<nj; j++)
    {
      for( i=0; i<ni; i++)
      {
        n = IJK2N( i, j, k, ni, nj);
        sforce[0] = lattice->force[0][n].sforce;
        sforce[1] = lattice->force[1][n].sforce;

        sum_x = 0.;
        sum_y = 0.;
        sum_z = 0.;

        if( !( lattice->solids[0][ n].is_solid ))
        {
          for( a=1; a<=6; a++)
          {
            ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
            ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
#if PARALLEL
        if( k+vz[a]<nk &&  k+vz[a]>=0 )
        {
        ka = k+vz[a];
        na = IJK2N( ia, ja, ka, ni, nj);
              if( ( lattice->solids[0][ na].is_solid ))
               {
                sum_x += WM*vx[a];
                sum_y += WM*vy[a];
                sum_z += WM*vz[a];
               }

        }
        else if(k+vz[a]>=nk)
        {
              if( ( lattice->process.z_neg_solid_to_recv[ja*ni+ia] ))
               {
                sum_x += WM*vx[a];
                sum_y += WM*vy[a];
                sum_z += WM*vz[a];
         }

        }
        else
        {
              if( ( lattice->process.z_pos_solid_to_recv[ja*ni+ia] ))
               {
                sum_x += WM*vx[a];
                sum_y += WM*vy[a];
                sum_z += WM*vz[a];
         }
        }

#else
            ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
            if(  ( lattice->solids[0][
                  IJK2N( ia, ja, ka, ni, nj)].is_solid ))
            {
              sum_x += WM*vx[a];
              sum_y += WM*vy[a];
              sum_z += WM*vz[a];
            }
#endif

          } // a = 0....6

          for( a=7; a<Q; a++)
          {
            ia = ( i+vx[a]<ni)?( ( i+vx[a]>=0)?( i+vx[a]):( ni-1)):( 0   );
            ja = ( j+vy[a]<nj)?( ( j+vy[a]>=0)?( j+vy[a]):( nj-1)):( 0   );
#if PARALLEL
        if( k+vz[a]<nk &&  k+vz[a]>=0 )
        {
        ka = k+vz[a];
        na = IJK2N( ia, ja, ka, ni, nj);
              if( ( lattice->solids[0][ na].is_solid ))
               {
                sum_x += WD*vx[a];
                sum_y += WD*vy[a];
                sum_z += WD*vz[a];
               }

        }
        else if(k+vz[a]>=nk)
        {
              if( ( lattice->process.z_neg_solid_to_recv[ja*ni+ia] ))
               {
                sum_x += WD*vx[a];
                sum_y += WD*vy[a];
                sum_z += WD*vz[a];
         }

        }
        else
        {
              if( ( lattice->process.z_pos_solid_to_recv[ja*ni+ia] ))
               {
          sum_x += WD*vx[a];
                sum_y += WD*vy[a];
                sum_z += WD*vz[a];
         }
        }

#else
            ka = ( k+vz[a]<nk)?( ( k+vz[a]>=0)?( k+vz[a]):( nk-1)):( 0   );
            if(  ( lattice->solids[0][
                  IJK2N( ia, ja, ka, ni, nj)].is_solid ))
            {
              sum_x += WD*vx[a];
              sum_y += WD*vy[a];
              sum_z += WD*vz[a];
            }
#endif
          } //a = 7,....18

#if 0
    if( lattice->param.big_V0_solid[0]+lattice->param.big_V0_solid[1] != 0.)
    {
       printf("%s %d >> %f %f\n",
         __FILE__,__LINE__,
         lattice->param.big_V0_solid[0],
         lattice->param.big_V0_solid[1] );
    }
#endif

//          sforce[0][0] = -lattice->param.big_V0_solid[0]*( sum_x) *lattice->macro_vars[0][n].rho *lattice->param.rhow;
//          sforce[0][1] = -lattice->param.big_V0_solid[0]*( sum_y) *lattice->macro_vars[0][n].rho *lattice->param.rhow;
//          sforce[0][2] = -lattice->param.big_V0_solid[0]*( sum_z) *lattice->macro_vars[0][n].rho *lattice->param.rhow;
          sforce[0][0] = -lattice->param.big_V0_solid[0]*( sum_x)  *lattice->param.rhow;
          sforce[0][1] = -lattice->param.big_V0_solid[0]*( sum_y)  *lattice->param.rhow;
          sforce[0][2] = -lattice->param.big_V0_solid[0]*( sum_z)  *lattice->param.rhow;

//          sforce[1][0] = -lattice->param.big_V0_solid[1]*( sum_x) *lattice->macro_vars[1][n].rho *lattice->param.rhow;
//          sforce[1][1] = -lattice->param.big_V0_solid[1]*( sum_y) *lattice->macro_vars[1][n].rho *lattice->param.rhow;
//          sforce[1][2] = -lattice->param.big_V0_solid[1]*( sum_z) *lattice->macro_vars[1][n].rho *lattice->param.rhow;
          sforce[1][0] = -lattice->param.big_V0_solid[1]*( sum_x)  *lattice->param.rhow;
          sforce[1][1] = -lattice->param.big_V0_solid[1]*( sum_y)  *lattice->param.rhow;
          sforce[1][2] = -lattice->param.big_V0_solid[1]*( sum_z)  *lattice->param.rhow;

        } /* if( !( lattice->solids[0][ j*ni+i].is_solid )) */

        else
        {
          sforce[0][0] = 0.;
          sforce[0][1] = 0.;
          sforce[0][2] = 0.;

          sforce[1][0] = 0.;
          sforce[1][1] = 0.;
          sforce[1][2] = 0.;
        }
      } /* for( i=0; i<ni; i++) */
    } /* for( j=0; j<nj; j++) */
  } /* for( k=0; k<nk; k++) */

#if PARALLEL
//  solid_send_recv_end( lattice, 0);
#endif

#else
// TODO
#endif
} /* void compute_real_fluid_solid_force( lattice_ptr lattice) */

#endif /* NON_LOCAL_FORCES */

void compute_max_rho( lattice_ptr lattice, int subs)
{
#if 0
  int n;
  real rho;

  *max_rho = 0.;

  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      rho = fabs( get_rho( lattice, subs, n));
      if( rho > *max_rho) { *max_rho = rho;}
    }
  }

  process_allreduce_real_max( lattice, max_rho);

#else
  int n;
  real *max_rho = get_max_rho_ptr(lattice,subs);
  real rho;
  *max_rho = 0.;
  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      rho = fabs( get_rho( lattice, subs, n));
      if( rho > *(max_rho)) { *max_rho = rho;}
    }
  }

  process_allreduce_real_max( lattice, max_rho);
#endif
} /* void compute_max_rho( lattice_ptr lattice, real *max_rho, int subs) */

void compute_min_rho( lattice_ptr lattice, int subs)
{
#if 0
  int n;
  real rho;

  compute_max_rho( lattice, min_rho, subs);

  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      rho = fabs( get_rho( lattice, subs, n));
      if( rho < *min_rho) { *min_rho = rho;}
    }
  }

  process_allreduce_real_min( lattice, min_rho);

#else
  int n;
  real *min_rho = get_min_rho_ptr(lattice,subs);
  real rho;

  compute_max_rho( lattice, subs);

  *min_rho = get_max_rho(lattice,subs);

  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      rho = fabs( get_rho( lattice, subs, n));
      if( rho < *(min_rho)) { *min_rho = rho;}
    }
  }

  process_allreduce_real_min( lattice, min_rho);
#endif
} /* void compute_min_rho( lattice_ptr lattice, real *min_rho, int subs) */

void compute_ave_rho( lattice_ptr lattice, int subs)
{
#if 0
  int n;
  int nn; // Number of non-solid nodes.
  *ave_rho = 0.;
  nn = 0;
  for( n=0; n<get_NumNodes(lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      *ave_rho += fabs( get_rho( lattice, subs, n));
      nn++;
    }
  }

  process_allreduce_real_sum( lattice, ave_rho);
  process_allreduce_int_sum( lattice, &nn);

  if( nn != 0) { *ave_rho = (*ave_rho) / nn;}

#else
  int n;
  int nn; // Number of non-solid nodes.
  real *ave_rho = get_ave_rho_ptr(lattice,subs);
  *ave_rho = 0.;
  nn = 0;
  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      *ave_rho += get_rho( lattice, subs, n);
      nn++;
    }
  }

  process_allreduce_real_sum( lattice, ave_rho);
  process_allreduce_int_sum( lattice, &nn);

  if( nn != 0)
  {
    *ave_rho = (*ave_rho)/nn;
  }
#endif
} /* void compute_ave_rho( lattice_ptr lattice, real *ave_rho, int subs) */

void compute_max_u( lattice_ptr lattice, int subs)
{
  int n;
  real *max_ux = get_max_ux_ptr(lattice,subs);
  real *max_uy = get_max_uy_ptr(lattice,subs);
  real *max_uz = get_max_uz_ptr(lattice,subs);
  real u;
  *max_ux = 0.;
  *max_uy = 0.;
  *max_uz = 0.;
  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      u = fabs( get_ux( lattice, subs, n));
      if( u > *(max_ux)) { *max_ux = u;}

      u = fabs( get_uy( lattice, subs, n));
      if( u > *(max_uy)) { *max_uy= u;}

      u = fabs( get_uz( lattice, subs, n));
      if( u > *(max_uz)) { *max_uz= u;}
    }
  }

  process_allreduce_real_max( lattice, max_ux);
  process_allreduce_real_max( lattice, max_uy);
  process_allreduce_real_max( lattice, max_uz);

} /* void compute_max_u( lattice_ptr lattice, real *max_u, int subs) */

void compute_min_u( lattice_ptr lattice, int subs)
{
  int n;
  real *min_ux = get_min_ux_ptr(lattice,subs);
  real *min_uy = get_min_uy_ptr(lattice,subs);
  real *min_uz = get_min_uz_ptr(lattice,subs);
  real u;

  compute_max_u( lattice, subs);

  *min_ux = get_max_ux(lattice,subs);
  *min_uy = get_max_uy(lattice,subs);
  *min_uz = get_max_uz(lattice,subs);

  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      u = fabs( get_ux( lattice, subs, n));
      if( u < *(min_ux)) { *min_ux = u;}

      u = fabs( get_uy( lattice, subs, n));
      if( u < *(min_uy)) { *min_uy = u;}

      u = fabs( get_uz( lattice, subs, n));
      if( u < *(min_uz)) { *min_uz = u;}
    }
  }

  process_allreduce_real_min( lattice, min_ux);
  process_allreduce_real_min( lattice, min_uy);
  process_allreduce_real_min( lattice, min_uz);

} /* void compute_min_u( lattice_ptr lattice, real *min_u) */

void compute_ave_u( lattice_ptr lattice, int subs)
{
#if 0
  int n;
  int nn; // Number of non-solid nodes.
  real *ave_ux = ave_u+0;
  real *ave_uy = ave_u+1;
  real *ave_uz = ave_u+2;
  *ave_ux = 0.;
  *ave_uy = 0.;
  *ave_uz = 0.;
  nn = 0;
  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      *ave_ux += get_ux( lattice, subs, n);
      *ave_uy += get_uy( lattice, subs, n);
      *ave_uz += get_uz( lattice, subs, n);
      nn++;
    }
  }

  process_allreduce_real_sum( lattice, ave_ux);
  process_allreduce_real_sum( lattice, ave_uy);
  process_allreduce_real_sum( lattice, ave_uz);
  process_allreduce_int_sum( lattice, &nn);

  if( nn != 0)
  {
    *ave_ux = (*ave_ux)/nn;
    *ave_uy = (*ave_uy)/nn;
    *ave_uz = (*ave_uz)/nn;
  }
#else
  int n;
  int nn; // Number of non-solid nodes.
  real *ave_ux = get_ave_ux_ptr(lattice,subs);
  real *ave_uy = get_ave_uy_ptr(lattice,subs);
  real *ave_uz = get_ave_uz_ptr(lattice,subs);
  *ave_ux = 0.;
  *ave_uy = 0.;
  *ave_uz = 0.;
  nn = 0;
  for( n=0; n<get_NumNodes( lattice); n++)
  {
    if( is_not_solid( lattice, n))
    {
      *ave_ux += get_ux( lattice, subs, n);
      *ave_uy += get_uy( lattice, subs, n);
      *ave_uz += get_uz( lattice, subs, n);
      nn++;
    }
  }

  process_allreduce_real_sum( lattice, ave_ux);
  process_allreduce_real_sum( lattice, ave_uy);
  process_allreduce_real_sum( lattice, ave_uz);
  process_allreduce_int_sum( lattice, &nn);

  if( nn != 0)
  {
    *ave_ux = (*ave_ux)/nn;
    *ave_uy = (*ave_uy)/nn;
    *ave_uz = (*ave_uz)/nn;
  }
#endif
} /* void compute_ave_u( lattice_ptr lattice, real *ave_u, int subs) */

void compute_flux( lattice_ptr lattice, int subs)
{
#if 0
  int n, nn;
  real rho, u_x, u_y, u_z;
  *(flux+0) = 0.;
  *(flux+1) = 0.;
  *(flux+2) = 0.;
  *(flux+3) = 0.;
  nn = 0;
  for( n=0; n<lattice->NumNodes; n++)
  {
    rho = lattice->macro_vars[subs][n].rho;
    u_x = lattice->macro_vars[subs][n].u[0];
    u_y = lattice->macro_vars[subs][n].u[1];
    u_z = lattice->macro_vars[subs][n].u[2];

    if( is_not_solid(lattice,subs))
    {
      *(flux+0) += rho*sqrt(u_x*u_x+u_y*u_y+u_z*u_z);
      *(flux+1) += rho*u_x;
      *(flux+2) += rho*u_y;
      *(flux+3) += rho*u_z;
      nn++;
    }
  }
  if( nn != 0)
  {
    *(flux+0) = (*(flux+0))/nn;
    *(flux+1) = (*(flux+1))/nn;
    *(flux+2) = (*(flux+2))/nn;
    *(flux+3) = (*(flux+3))/nn;
  }

  process_allreduce_real_sum( lattice, flux+0);
  process_allreduce_real_sum( lattice, flux+1);
  process_allreduce_real_sum( lattice, flux+2);
  process_allreduce_real_sum( lattice, flux+3);

  flux[0] /= get_num_procs(lattice);
  flux[1] /= get_num_procs(lattice);
  flux[2] /= get_num_procs(lattice);
  flux[3] /= get_num_procs(lattice);

#else
// TODO
#endif
} /* void compute_flux( lattice_ptr lattice, real *flux, int subs) */

#if STORE_UEQ
//void compute_max_ueq( lattice_ptr lattice, real *max_u)
void compute_max_ueq( lattice_ptr lattice, real *max_u)
{
#if 0
  int n;
  *max_u = 0.;
  *(max_u+1) = 0.;
  *(max_u+2) = 0.;
  for( n=0; n<lattice->NumNodes; n++)
  {
    if( !( lattice->solids[0][n].is_solid))
    {
      if( fabs( lattice->ueq[n].u[0]) > *(max_u))
      {
        *max_u = fabs( lattice->ueq[n].u[0]);
      }
      if( fabs( lattice->ueq[n].u[1]) > *(max_u+1))
      {
        *(max_u+1) = fabs( lattice->ueq[n].u[1]);
      }
      if( fabs( lattice->ueq[n].u[2]) > *(max_u+2))
      {
        *(max_u+2) = fabs( lattice->ueq[n].u[2]);
      }
    }
  }

  process_allreduce_real_max( lattice, max_u+0);
  process_allreduce_real_max( lattice, max_u+1);
  process_allreduce_real_max( lattice, max_u+2);

#else
// TODO
#endif
} /* void compute_max_ueq( lattice_ptr lattice, real *max_u) */

//void compute_min_ueq( lattice_ptr lattice, real *min_u)
void compute_min_ueq( lattice_ptr lattice, real *min_u)
{
#if 0
  int n;
  compute_max_ueq( lattice, min_u);
  for( n=0; n<lattice->NumNodes; n++)
  {
    if( !( lattice->solids[0][n].is_solid))
    {
      if( fabs( lattice->ueq[n].u[0]) < *(min_u))
      {
        *min_u = fabs( lattice->ueq[n].u[0]);
      }
      if( fabs( lattice->ueq[n].u[1]) < *(min_u+1))
      {
        *(min_u+1) = fabs( lattice->ueq[n].u[1]);
      }
      if( fabs( lattice->ueq[n].u[2]) < *(min_u+2))
      {
        *(min_u+2) = fabs( lattice->ueq[n].u[2]);
      }
    }
  }

  process_allreduce_real_min( lattice, min_u+0);
  process_allreduce_real_min( lattice, min_u+1);
  process_allreduce_real_min( lattice, min_u+2);

#else
// TODO
#endif
} /* void compute_min_ueq( lattice_ptr lattice, real *min_u) */

//void compute_ave_ueq( lattice_ptr lattice, real *ave_u)
void compute_ave_ueq( lattice_ptr lattice, real *ave_u)
{
#if 0
  int n, nn;
  *ave_u = 0.;
  *(ave_u+1) = 0.;
  *(ave_u+2) = 0.;
  nn = 0;
  for( n=0; n<lattice->NumNodes; n++)
  {
    if( !( lattice->solids[0][n].is_solid))
    {
      *ave_u += fabs( lattice->ueq[n].u[0]);
      *(ave_u+1) += fabs( lattice->ueq[n].u[1]);
      *(ave_u+2) += fabs( lattice->ueq[n].u[2]);
      nn++;
    }
  }

  process_allreduce_real_sum( lattice, ave_u+0);
  process_allreduce_real_sum( lattice, ave_u+1);
  process_allreduce_real_sum( lattice, ave_u+2);
  process_allreduce_int_sum( lattice, &nn);

  if( nn != 0)
  {
    *ave_u     = (*ave_u)/nn;
    *(ave_u+1) = (*(ave_u+1))/nn;
    *(ave_u+2) = (*(ave_u+2))/nn;
  }

#else
// TODO
#endif
} /* void compute_ave_ueq( lattice_ptr lattice, real *ave_u) */

#endif /* STORE_UEQ */

