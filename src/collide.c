//##############################################################################
//
// collide.c
//

//##############################################################################
//
// void collide( lattice_ptr lattice)
void collide( lattice_ptr lattice)
{
  int n;
  int subs;
  int ni = get_ni(lattice);
  int nj = get_nj(lattice);
  int nk = get_nk(lattice);
  int i, j, k, a;

  real rho;
  real ux;
  real uy;
  real uz;

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

  real temp;

  real** fptr;

  fptr = (real**)malloc( sizeof(real*)*get_NumVelDirs(lattice));

  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
    n=0;
    for( k=0; k<nk; k++)
    {
      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          fptr[C ] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, C );
          fptr[E ] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, E );
          fptr[W ] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, W );
          fptr[N ] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, N );
          fptr[S ] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, S );
          fptr[NE] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, NE);
          fptr[SW] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, SW);
          fptr[NW] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, NW);
          fptr[SE] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, SE);
         if( get_NumDims(lattice)==3)
         {
          fptr[T ] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, T );
          fptr[B ] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, B );
          fptr[TE] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, TE);
          fptr[BW] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, BW);
          fptr[TW] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, TW);
          fptr[BE] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, BE);
          fptr[TN] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, TN);
          fptr[BS] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, BS);
          fptr[TS] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, TS);
          fptr[BN] = get_fptr(lattice,subs, i, j, k, 0, 0, 0, BN);
         }

          rho = 0.;
          ux = 0.;
          uy = 0.;
          uz = 0.;
          rho+=*(fptr[0]);
          for( a=1; a<get_NumVelDirs(lattice); a++)
          {
            rho+=*(fptr[a]);
            ux += (*(fptr[a]))*vx[a+((a%2)?(1):(-1))];
            uy += (*(fptr[a]))*vy[a+((a%2)?(1):(-1))];
            if( get_NumDims(lattice)==3)
            {
              uz += (*(fptr[a]))*vz[a+((a%2)?(1):(-1))];
            }
          }
          ux /= rho;
          uy /= rho;
          uz /= rho;

          if( !skip_collision_step( lattice))
          {
            if( !skip_body_force_term( lattice))
            {
              ux += get_gaccel_ux( lattice, subs);
              uy += get_gaccel_uy( lattice, subs);
              uz += get_gaccel_uz( lattice, subs);
            }

            usq = ux*ux + uy*uy + uz*uz;

#if 1
            *(fptr[0]) = *(fptr[0])*(1-1./get_tau(lattice,subs))
                      + ( /*feq[a]*/
                          W0 * rho*(1. - 1.5*usq)
                        ) / get_tau(lattice,subs);

            for( a=1; a<=4; a++)
            {
              udotx = ((real)vx[a+((a%2)?(1):(-1))]*ux
                      +(real)vy[a+((a%2)?(1):(-1))]*uy
                      +(real)vz[a+((a%2)?(1):(-1))]*uz);

              *(fptr[a]) = *(fptr[a])*(1-1./get_tau(lattice,subs))
                      + ( /*feq[a]*/
                          W1*rho*(1. + 3.*udotx + 4.5 *udotx*udotx - 1.5*usq)
                        ) / get_tau(lattice,subs);
            }

            for( ; a<=8; a++)
            {
              udotx = ((real)vx[a+((a%2)?(1):(-1))]*ux
                      +(real)vy[a+((a%2)?(1):(-1))]*uy
                      +(real)vz[a+((a%2)?(1):(-1))]*uz);

              *(fptr[a]) = *(fptr[a])*(1-1./get_tau(lattice,subs))
                      + ( /*feq[a]*/
                          W2*rho*(1. + 3.*udotx + 4.5 *udotx*udotx - 1.5*usq)
                        ) / get_tau(lattice,subs);
            }

            if( get_NumDims(lattice)==3)
            {
              for( ; a<=10; a++)
              {
                udotx = ((real)vx[a+((a%2)?(1):(-1))]*ux
                        +(real)vy[a+((a%2)?(1):(-1))]*uy
                        +(real)vz[a+((a%2)?(1):(-1))]*uz);

                *(fptr[a]) = *(fptr[a])*(1-1./get_tau(lattice,subs))
                        + ( /*feq[a]*/
                            W1*rho*(1. + 3.*udotx + 4.5 *udotx*udotx - 1.5*usq)
                          ) / get_tau(lattice,subs);
              }

              for( ; a<get_NumVelDirs(lattice); a++)
              {
                udotx = ((real)vx[a+((a%2)?(1):(-1))]*ux
                        +(real)vy[a+((a%2)?(1):(-1))]*uy
                        +(real)vz[a+((a%2)?(1):(-1))]*uz);

                *(fptr[a]) = *(fptr[a])*(1-1./get_tau(lattice,subs))
                        + ( /*feq[a]*/
                            W2*rho*(1. + 3.*udotx + 4.5 *udotx*udotx - 1.5*usq)
                          ) / get_tau(lattice,subs);
              }
            }
#else
            // Just assign the weighted rhos for debugging.
            *(fptr[0]) = W0*rho;

            for( a=1; a<=4; a++)
            {
              udotx = ((real)vx[a+((a%2)?(1):(-1))]*ux
                      +(real)vy[a+((a%2)?(1):(-1))]*uy
                      +(real)vz[a+((a%2)?(1):(-1))]*uz);

              *(fptr[a]) = W1*rho;
            }

            for( ; a<=8; a++)
            {
              udotx = ((real)vx[a+((a%2)?(1):(-1))]*ux
                      +(real)vy[a+((a%2)?(1):(-1))]*uy
                      +(real)vz[a+((a%2)?(1):(-1))]*uz);

              *(fptr[a]) = W2*rho;
            }

            if( get_NumDims(lattice)==3)
            {
              for( ; a<=10; a++)
              {
                udotx = ((real)vx[a+((a%2)?(1):(-1))]*ux
                        +(real)vy[a+((a%2)?(1):(-1))]*uy
                        +(real)vz[a+((a%2)?(1):(-1))]*uz);

                *(fptr[a]) = W1*rho;
              }

              for( ; a<get_NumVelDirs(lattice); a++)
              {
                udotx = ((real)vx[a+((a%2)?(1):(-1))]*ux
                        +(real)vy[a+((a%2)?(1):(-1))]*uy
                        +(real)vz[a+((a%2)?(1):(-1))]*uz);

                *(fptr[a]) = W2*rho;
              }
            }
#endif
          }

          // Swap opposing fs to undo the swap in stream_collide_stream.
          for( a=1; a<get_NumVelDirs(lattice); a+=2)
          {
            temp = *(fptr[a]);
            *(fptr[a]) = *(fptr[a+1]);
            *(fptr[a+1]) = temp;
          }

          rho = 0.;
          ux = 0.;
          uy = 0.;
          uz = 0.;
          for( a=0; a<get_NumVelDirs(lattice); a++)
          {
            rho+=*(fptr[a]);
            ux += (*(fptr[a]))*vx[a];
            uy += (*(fptr[a]))*vy[a];
            if( get_NumDims(lattice)==3)
            {
              uz += (*(fptr[a]))*vz[a];
            }
          }

          ux /= rho;
          uy /= rho;
          uz /= rho;

          set_rho(lattice,subs,n,rho);
          set_ux(lattice,subs,n,ux);
          set_uy(lattice,subs,n,uy);
          set_uz(lattice,subs,n,uz);

          n++;
        }
      }
    }
  }

  free( fptr);

} /* void collide( lattice_ptr lattice) */

