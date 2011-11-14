
void stream_collide_stream( lattice_ptr lattice)
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

  if( !fptr)
  {
    printf("%s %d BOOM!\n",__FILE__,__LINE__);
    process_exit(1);
  }

  for( subs=0; subs<NUM_FLUID_COMPONENTS; subs++)
  {
    n = 0;
    for( k=0; k<nk; k++)
    {
      for( j=0; j<nj; j++)
      {
        for( i=0; i<ni; i++)
        {
          fptr[C ] = get_fptr(lattice,subs, i,j,k,-vx[C ],-vy[C ],-vz[C ],C );
          fptr[E ] = get_fptr(lattice,subs, i,j,k,-vx[E ],-vy[E ],-vz[E ],E );
          fptr[W ] = get_fptr(lattice,subs, i,j,k,-vx[W ],-vy[W ],-vz[W ],W );
          fptr[N ] = get_fptr(lattice,subs, i,j,k,-vx[N ],-vy[N ],-vz[N ],N );
          fptr[S ] = get_fptr(lattice,subs, i,j,k,-vx[S ],-vy[S ],-vz[S ],S );
          fptr[NE] = get_fptr(lattice,subs, i,j,k,-vx[NE],-vy[NE],-vz[NE],NE);
          fptr[SW] = get_fptr(lattice,subs, i,j,k,-vx[SW],-vy[SW],-vz[SW],SW);
          fptr[NW] = get_fptr(lattice,subs, i,j,k,-vx[NW],-vy[NW],-vz[NW],NW);
          fptr[SE] = get_fptr(lattice,subs, i,j,k,-vx[SE],-vy[SE],-vz[SE],SE);
         if( get_NumDims(lattice)==3)
         {
          fptr[T ] = get_fptr(lattice,subs, i,j,k,-vx[T ],-vy[T ],-vz[T ],T );
          fptr[B ] = get_fptr(lattice,subs, i,j,k,-vx[B ],-vy[B ],-vz[B ],B );
          fptr[TE] = get_fptr(lattice,subs, i,j,k,-vx[TE],-vy[TE],-vz[TE],TE);
          fptr[BW] = get_fptr(lattice,subs, i,j,k,-vx[BW],-vy[BW],-vz[BW],BW);
          fptr[TW] = get_fptr(lattice,subs, i,j,k,-vx[TW],-vy[TW],-vz[TW],TW);
          fptr[BE] = get_fptr(lattice,subs, i,j,k,-vx[BE],-vy[BE],-vz[BE],BE);
          fptr[TN] = get_fptr(lattice,subs, i,j,k,-vx[TN],-vy[TN],-vz[TN],TN);
          fptr[BS] = get_fptr(lattice,subs, i,j,k,-vx[BS],-vy[BS],-vz[BS],BS);
          fptr[TS] = get_fptr(lattice,subs, i,j,k,-vx[TS],-vy[TS],-vz[TS],TS);
          fptr[BN] = get_fptr(lattice,subs, i,j,k,-vx[BN],-vy[BN],-vz[BN],BN);
         }

#if 1
          rho = 0.;
          ux = 0.;
          uy = 0.;
          uz = 0.;
          for( a=0; a<get_NumVelDirs(lattice); a++)
          {
            rho+= (*(fptr[a]));
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

          if( !skip_collision_step( lattice))
          {
            ux += get_gaccel_ux( lattice, subs);
            uy += get_gaccel_uy( lattice, subs);
            uz += get_gaccel_uz( lattice, subs);

            usq = ux*ux + uy*uy + uz*uz;

#if 1
            *(fptr[0]) = *(fptr[0])*(1-1./get_tau(lattice,subs))
                      + ( /*feq[a]*/
                          W0 * rho*(1. - 1.5*usq)
                        ) / get_tau(lattice,subs);

            for( a=1; a<=4; a++)
            {
              udotx = ((real)vx[a]*ux+(real)vy[a]*uy+(real)vz[a]*uz);

              *(fptr[a]) = *(fptr[a])*(1-1./get_tau(lattice,subs))
                      + ( /*feq[a]*/
                          W1*rho*(1. + 3.*udotx + 4.5 *udotx*udotx - 1.5*usq)
                        ) / get_tau(lattice,subs);
            }

            for( ; a<=8; a++)
            {
              udotx = ((real)vx[a]*ux+(real)vy[a]*uy+(real)vz[a]*uz);

              *(fptr[a]) = *(fptr[a])*(1-1./get_tau(lattice,subs))
                      + ( /*feq[a]*/
                          W2*rho*(1. + 3.*udotx + 4.5 *udotx*udotx - 1.5*usq)
                        ) / get_tau(lattice,subs);
            }

            if( get_NumDims(lattice)==3)
            {
              for( ; a<=10; a++)
              {
                udotx = ((real)vx[a]*ux+(real)vy[a]*uy+(real)vz[a]*uz);

                *(fptr[a]) = *(fptr[a])*(1-1./get_tau(lattice,subs))
                        + ( /*feq[a]*/
                            W1*rho*(1. + 3.*udotx + 4.5 *udotx*udotx - 1.5*usq)
                          ) / get_tau(lattice,subs);
              }

              for( ; a<get_NumVelDirs(lattice); a++)
              {
                udotx = ((real)vx[a]*ux+(real)vy[a]*uy+(real)vz[a]*uz);

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
              udotx = ((real)vx[a]*ux+(real)vy[a]*uy+(real)vz[a]*uz);

              *(fptr[a]) = W1*rho;
            }

            for( ; a<=8; a++)
            {
              udotx = ((real)vx[a]*ux+(real)vy[a]*uy+(real)vz[a]*uz);

              *(fptr[a]) = W2*rho;
            }

            if( get_NumDims(lattice)==3)
            {
              for( ; a<=10; a++)
              {
                udotx = ((real)vx[a]*ux+(real)vy[a]*uy+(real)vz[a]*uz);

                *(fptr[a]) = W1*rho;
              }

              for( ; a<get_NumVelDirs(lattice); a++)
              {
                udotx = ((real)vx[a]*ux+(real)vy[a]*uy+(real)vz[a]*uz);

                *(fptr[a]) = W2*rho;
              }
            }
#endif
          }

          // Swap fs of opposing neighbors, thus streaming into the correct
          // neighbor nodes although storing in the opposite direction.
          // The collide step that follows this stream_collide_stream step
          // will know that the fs are all pointing backwards.
          for( a=1; a<get_NumVelDirs(lattice); a+=2)
          {
            temp = *(fptr[a]);
            *(fptr[a]) = *(fptr[a+1]);
            *(fptr[a+1]) = temp;
          }
#else
#if 0
          if( i==1 && j==1)
          {
            *(fptr[C ]) = 0.12345678;
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
#endif

          n++;
        }
      }
    }
  }

  free( fptr);
}
