
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
          fptr[C ] = get_fptr(lattice,subs, i-vx[C ], j-vy[C ], k-vz[C ], C );
          fptr[E ] = get_fptr(lattice,subs, i-vx[E ], j-vy[E ], k-vz[E ], E );
          fptr[W ] = get_fptr(lattice,subs, i-vx[W ], j-vy[W ], k-vz[W ], W );
          fptr[N ] = get_fptr(lattice,subs, i-vx[N ], j-vy[N ], k-vz[N ], N );
          fptr[S ] = get_fptr(lattice,subs, i-vx[S ], j-vy[S ], k-vz[S ], S );
          fptr[NE] = get_fptr(lattice,subs, i-vx[NE], j-vy[NE], k-vz[NE], NE);
          fptr[SW] = get_fptr(lattice,subs, i-vx[SW], j-vy[SW], k-vz[SW], SW);
          fptr[NW] = get_fptr(lattice,subs, i-vx[NW], j-vy[NW], k-vz[NW], NW);
          fptr[SE] = get_fptr(lattice,subs, i-vx[SE], j-vy[SE], k-vz[SE], SE);
         if( get_NumDims(lattice)==3)
         {
          fptr[T ] = get_fptr(lattice,subs, i-vx[T ], j-vy[T ], k-vz[T ], T );
          fptr[B ] = get_fptr(lattice,subs, i-vx[B ], j-vy[B ], k-vz[B ], B );
          fptr[TE] = get_fptr(lattice,subs, i-vx[TE], j-vy[TE], k-vz[TE], TE);
          fptr[BW] = get_fptr(lattice,subs, i-vx[BW], j-vy[BW], k-vz[BW], BW);
          fptr[TW] = get_fptr(lattice,subs, i-vx[TW], j-vy[TW], k-vz[TW], TW);
          fptr[BE] = get_fptr(lattice,subs, i-vx[BE], j-vy[BE], k-vz[BE], BE);
          fptr[TN] = get_fptr(lattice,subs, i-vx[TN], j-vy[TN], k-vz[TN], TN);
          fptr[BS] = get_fptr(lattice,subs, i-vx[BS], j-vy[BS], k-vz[BS], BS);
          fptr[TS] = get_fptr(lattice,subs, i-vx[TS], j-vy[TS], k-vz[TS], TS);
          fptr[BN] = get_fptr(lattice,subs, i-vx[BN], j-vy[BN], k-vz[BN], BN);
         }

          if( 1)//get_time(lattice)>1)
          {
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
          }
          else
          {
            // For first time step, use initial values of macroscopic
            // variables.
            rho = get_rho(lattice,subs,n);
            ux  = get_ux (lattice,subs,n);
            uy  = get_uy (lattice,subs,n);
            uz  = get_uz (lattice,subs,n);
          }

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

          n++;
        }
      }
    }
  }

  free( fptr);
}
