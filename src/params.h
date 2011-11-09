//##############################################################################
//
// params.h
//
//  - Routines for reading and writing params.
//
//  - Reads from params.in
//
//  - Writes to params.dat
//

// void skip_label( FILE *in)
void skip_label( FILE *in)
{
  char c;

  c = fgetc( in);

  if( ( c >= '0' && c <='9') || c == '-' || c == '.')
  {
    // Digit or Sign or Decimal
    ungetc( c, in);
  }
  else
  {
    while( ( c = fgetc( in)) != ' ');
    while( ( c = fgetc( in)) == ' ');
    ungetc( c, in);
  }

} /* void skip_label( FILE *in) */

// void read_params( struct lattice_struct *lattice)
//##############################################################################
//
// R E A D   P A R A M S
//
//  - Read the problem parameters from a file.
//
void read_params( lattice_ptr lattice, const char *infile)
{
  FILE   *in;
  int    blank;
  real dblank;

  // Store different scanf format specifier for single and double precision.
  const char* rspec = ((sizeof(real)==sizeof(float))?("%f"):("%lf"));

//default
  lattice->param.GZL = 0;
  lattice->param.PressureBC = 0;
//default

  if( !( in = fopen( infile, "r")))
  {
    printf("%s %d %04d>> ERROR: "
      "fopen(\"%s\",\"r\") = NULL.  Bye, bye!\n",
      __FILE__, __LINE__,
      get_proc_id( lattice),
      infile);
    process_exit(1);
  }

  skip_label( in); fscanf( in, "%d", &( lattice->param.LX)             );
  skip_label( in); fscanf( in, "%d", &( lattice->param.LY)             );
  skip_label( in); fscanf( in, "%d", &( lattice->param.LZ)             );
  skip_label( in); fscanf( in, "%d", &( lattice->param.length_scale    ));
  skip_label( in); fscanf( in, "%d", &( lattice->param.NumFrames)      );
  skip_label( in); fscanf( in, "%d", &( lattice->param.FrameRate)      );
  skip_label( in); fscanf( in, rspec,&( lattice->param.tau[0])         );
  skip_label( in); fscanf( in, rspec,   lattice->param.gforce[0]       );
  skip_label( in); fscanf( in, rspec,   lattice->param.gforce[0]+1     );
  skip_label( in); fscanf( in, rspec,   lattice->param.gforce[0]+2     );
  skip_label( in); fscanf( in, "%d",    lattice->param.end_grav        );
  if( NUM_FLUID_COMPONENTS==2)
  {
  skip_label( in); fscanf( in, rspec,&( lattice->param.tau[1])         );
  skip_label( in); fscanf( in, rspec,   lattice->param.gforce[1]       );
  skip_label( in); fscanf( in, rspec,   lattice->param.gforce[1]+1     );
  skip_label( in); fscanf( in, rspec,   lattice->param.gforce[1]+2     );
  skip_label( in); fscanf( in, "%d",    lattice->param.end_grav+1      );
  }
  else if( NUM_FLUID_COMPONENTS==1)
  {
  skip_label( in); fscanf( in, rspec,   &dblank                        );
  skip_label( in); fscanf( in, rspec,   &dblank                        );
  skip_label( in); fscanf( in, rspec,   &dblank                        );
  skip_label( in); fscanf( in, rspec,   &dblank                        );
  skip_label( in); fscanf( in, "%d",    &blank                         );
  }
  else
  {
    printf(
      "read_params() -- "
      "Unhandled case "
      "NUM_FLUID_COMPONENTS = %d .  "
      "Exiting!\n",
      NUM_FLUID_COMPONENTS);
    process_exit(1);
  }
    skip_label( in); fscanf( in, "%d ", &(lattice->param.buoyancy      ) );
    skip_label( in); fscanf( in, "%d ", &(lattice->param.incompressible) );
    skip_label( in); fscanf( in, "%d ", &(lattice->param.simple_diffusion) );
//  skip_label( in); fscanf( in, rspec, &(lattice->param.rhow)           );
//    lattice->param.rhow = 1.0;
    skip_label( in); fscanf( in, rspec, &(  lattice->param.rho_A[0])        );

printf("rho_A = %f\n",lattice->param.rho_A[0]);

    skip_label( in); fscanf( in, rspec, &(  lattice->param.rho_B[0])        );
printf("rho_B = %f\n",lattice->param.rho_B[0]);

#if INAMURO_SIGMA_COMPONENT
    skip_label( in); fscanf( in, rspec,&( lattice->param.rho_sigma)      );
    skip_label( in); fscanf( in, rspec,&( lattice->param.rho_sigma_in)   );
    skip_label( in); fscanf( in, rspec,&( lattice->param.rho_sigma_out)  );
    skip_label( in); fscanf( in, rspec,&( lattice->param.u_sigma)        );
    skip_label( in); fscanf( in, rspec,&( lattice->param.u_sigma_in)     );
    skip_label( in); fscanf( in, rspec,&( lattice->param.u_sigma_out)    );
    skip_label( in); fscanf( in, "%d ",&( lattice->param.sigma_start)    );
    skip_label( in); fscanf( in, "%d ",&( lattice->param.sigma_stop )    );
    skip_label( in); fscanf( in, "%d ",&( lattice->param.sigma_btc_rate ));
    skip_label( in); fscanf( in, "%d ",&( lattice->param.sigma_btc_spot ));
#else /* !( INAMURO_SIGMA_COMPONENT) */
    skip_label( in); fscanf( in, rspec,   &dblank                        );
    skip_label( in); fscanf( in, rspec,   &dblank                        );
    skip_label( in); fscanf( in, rspec,   &dblank                        );
    skip_label( in); fscanf( in, rspec,   &dblank                        );
    skip_label( in); fscanf( in, rspec,   &dblank                        );
    skip_label( in); fscanf( in, rspec,   &dblank                        );
    skip_label( in); fscanf( in, "%d ",   &blank                         );
    skip_label( in); fscanf( in, "%d ",   &blank                         );
    skip_label( in); fscanf( in, "%d ",   &blank                         );
    skip_label( in); fscanf( in, "%d ",   &blank                         );
#endif /* INAMURO_SIGMA_COMPONENT */
    skip_label( in); fscanf( in, "%d ",&( lattice->param.GZL   )         );
    skip_label( in); fscanf( in, "%d ",&( lattice->param.PressureBC)     );
    skip_label( in); fscanf( in, "%d ",&( lattice->param.AllBoundaryPeriodic)) ;

    skip_label( in); fscanf( in, rspec,&( lattice->param.rho_in)         );
    skip_label( in); fscanf( in, rspec,&( lattice->param.rho_out)        );
    skip_label( in); fscanf( in, rspec,&( lattice->param.ux_in)          );
    skip_label( in); fscanf( in, rspec,&( lattice->param.ux_out)         );
    skip_label( in); fscanf( in, rspec,&( lattice->param.uy_in)          );
    skip_label( in); fscanf( in, rspec,&( lattice->param.uy_out)         );
    skip_label( in); fscanf( in, rspec,&( lattice->param.uz_in)          );
    skip_label( in); fscanf( in, rspec,&( lattice->param.uz_out)         );
    skip_label( in); fscanf( in, rspec,&( lattice->param.big_V0)         );
    skip_label( in); fscanf( in, rspec,   lattice->param.big_V0_solid+0  );

  if( NUM_FLUID_COMPONENTS==2)
  {
    skip_label( in); fscanf( in, rspec,   lattice->param.big_V0_solid+1  );
  }
  else if( NUM_FLUID_COMPONENTS==1)
  {
    skip_label( in); fscanf( in, rspec,  &dblank                         );
  }
  else
  {
    printf(
      "read_params() -- "
      "Unhandled case "
      "NUM_FLUID_COMPONENTS = %d .  "
      "Exiting!\n",
      NUM_FLUID_COMPONENTS);
    process_exit(1);
  }
    skip_label( in); fscanf( in, "%d", &( lattice->param.ns_flag       ) );
    skip_label( in); fscanf( in, rspec,&( lattice->param.ns            ) );
    skip_label( in); fscanf( in, "%d", &( lattice->param.ic_poisseuille) );
    skip_label( in); fscanf( in, "%d", &( lattice->param.bc_poisseuille) );
    skip_label( in); fscanf( in, "%d", &( lattice->param.bc_slip_north ) );
#if INAMURO_SIGMA_COMPONENT
    skip_label( in); fscanf( in, "%d", &( lattice->param.bc_sigma_slip ) );
#else /* !( INAMURO_SIGMA_COMPONENT) */
    skip_label( in); fscanf( in, "%d",   &blank                          );
#endif /* INAMURO_SIGMA_COMPONENT */
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_t_in+0 );
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_b_in+0 );
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_t_out+0);
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_b_out+0);
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_t_in+0 );
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_b_in+0 );
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_t_out+0);
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_b_out+0);
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_n_in+0 );
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_s_in+0 );
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_n_out+0);
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_s_out+0);
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_n_in+0 );
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_s_in+0 );
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_n_out+0);
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_s_out+0);
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_e_in+0 );
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_w_in+0 );
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_e_out+0);
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_w_out+0);
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_e_in+0 );
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_w_in+0 );
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_e_out+0);
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_w_out+0);
  if( NUM_FLUID_COMPONENTS==2)
  {

    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_t_in+1 );
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_b_in+1 );
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_t_out+1);
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_b_out+1);
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_t_in+1 );
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_b_in+1 );
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_t_out+1);
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_b_out+1);
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_n_in+1 );
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_s_in+1 );
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_n_out+1);
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_s_out+1);
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_n_in+1 );
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_s_in+1 );
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_n_out+1);
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_s_out+1);
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_e_in+1 );
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_w_in+1 );
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_e_out+1);
    skip_label( in); fscanf( in, "%d",    lattice->param.pressure_w_out+1);
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_e_in+1 );
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_w_in+1 );
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_e_out+1);
    skip_label( in); fscanf( in, "%d",    lattice->param.velocity_w_out+1);
  }
  else if( NUM_FLUID_COMPONENTS==1)
  {
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
  }
  else
  {
    printf(
      "read_params() -- "
      "Unhandled case "
      "NUM_FLUID_COMPONENTS = %d .  "
      "Exiting!\n",
      NUM_FLUID_COMPONENTS);
    process_exit(1);
  }
  if( INAMURO_SIGMA_COMPONENT)
  {
    skip_label( in); fscanf( in, "%d",    &lattice->param.constcon_t_in     );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constcon_b_in     );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constcon_t_out    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constcon_b_out    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constflx_t_in     );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constflx_b_in     );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constflx_t_out    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constflx_b_out    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constcon_n_in     );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constcon_s_in     );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constcon_n_out    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constcon_s_out    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constflx_n_in     );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constflx_s_in     );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constflx_n_out    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constflx_s_out    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constcon_e_in     );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constcon_w_in     );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constcon_e_out    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constcon_w_out    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constflx_e_in     );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constflx_w_in     );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constflx_e_out    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.constflx_w_out    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.zeroconcgrad_t    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.zeroconcgrad_b    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.zeroconcgrad_n    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.zeroconcgrad_s    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.zeroconcgrad_e    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.zeroconcgrad_w    );
    skip_label( in); fscanf( in, "%d",    &lattice->param.zeroconcgrad_full );
  }
  else
  {
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
    skip_label( in); fscanf( in, "%d",    &blank                         );
  }
  skip_label( in); fscanf( in, "%d", &( lattice->param.plot_scale_dynamic)   );
  skip_label( in); fscanf( in, "%d", &( lattice->param.use_colormap      )   );
  skip_label( in); fscanf( in, "%d", &( lattice->param.initial_condition )   );
#if INTEGER_IC_BOUND
  skip_label( in); fscanf( in, "%d",&( lattice->param.x0                )   );
  skip_label( in); fscanf( in, "%d",&( lattice->param.y0                )   );
  skip_label( in); fscanf( in, "%d",&( lattice->param.z0                )   );
  skip_label( in); fscanf( in, "%d",&( lattice->param.r0                )   );
#else
  skip_label( in); fscanf( in, rspec,&( lattice->param.x0                )   );
  skip_label( in); fscanf( in, rspec,&( lattice->param.y0                )   );
  skip_label( in); fscanf( in, rspec,&( lattice->param.z0                )   );
  skip_label( in); fscanf( in, rspec,&( lattice->param.r0                )   );
#endif
  skip_label( in); fscanf( in, rspec,&( lattice->param.cut               )   );
#if INTEGER_IC_BOUND
  skip_label( in); fscanf( in, "%d",&( lattice->param.x1                )   );
  skip_label( in); fscanf( in, "%d",&( lattice->param.x2                )   );
  skip_label( in); fscanf( in, "%d",&( lattice->param.y1                )   );
  skip_label( in); fscanf( in, "%d",&( lattice->param.y2                )   );
  skip_label( in); fscanf( in, "%d",&( lattice->param.z1                )   );
  skip_label( in); fscanf( in, "%d",&( lattice->param.z2                )   );
#else
  skip_label( in); fscanf( in, rspec,&( lattice->param.x1                )   );
  skip_label( in); fscanf( in, rspec,&( lattice->param.x2                )   );
  skip_label( in); fscanf( in, rspec,&( lattice->param.y1                )   );
  skip_label( in); fscanf( in, rspec,&( lattice->param.y2                )   );
  skip_label( in); fscanf( in, rspec,&( lattice->param.z1                )   );
  skip_label( in); fscanf( in, rspec,&( lattice->param.z2                )   );
#endif
  skip_label( in); fscanf( in, rspec,&( lattice->param.rel_x1            )   );
  skip_label( in); fscanf( in, rspec,&( lattice->param.rel_x2            )   );
  skip_label( in); fscanf( in, rspec,&( lattice->param.rel_y1            )   );
  skip_label( in); fscanf( in, rspec,&( lattice->param.rel_y2            )   );
  skip_label( in); fscanf( in, rspec,&( lattice->param.rel_z1            )   );
  skip_label( in); fscanf( in, rspec,&( lattice->param.rel_z2            )   );
  skip_label( in); fscanf( in, "%d", &( lattice->param.dump_rho          )   );
  skip_label( in); fscanf( in, "%d", &( lattice->param.dump_u            )   );
  skip_label( in); fscanf( in, "%d", &( lattice->param.dump_force        )   );
  skip_label( in); fscanf( in, "%d", &( lattice->param.dump_vor          )   );

  if( NUM_DIMENSIONS==2 && get_nk(lattice)!=1)
  {
    set_nk(lattice,1);
    printf("WARNING: Setting nk=1 because NumDims=2.\n");
  }

  lattice->NumNodes =
    lattice->param.LX*lattice->param.LY*lattice->param.LZ;

  if( (lattice->param.FrameRate)%2)
  {
    printf("%s %d >> read_params() -- "
           "WARNING: "
           "Frame rate should be even because time steps are performed "
           "in pairs in this implementation. "
           "Increasing frame rate from %d to %d."
           "\n"
          , __FILE__,__LINE__
          , lattice->param.FrameRate
          , lattice->param.FrameRate+1
          );
    lattice->param.FrameRate++;
  }

  lattice->NumTimeSteps =
    lattice->param.NumFrames * lattice->param.FrameRate;

  if( NUM_FLUID_COMPONENTS==2)
  {
    lattice->param.rho_A[1] = lattice->param.rho_B[0];
    lattice->param.rho_B[1] = lattice->param.rho_A[0];
  }

 // Set default values for x0, y0 and r0 if they are negative.
#if INTEGER_IC_BOUND
  if( lattice->param.x0 < 0)
  {
    lattice->param.x0 = lattice->param.LX/2;
  }
  if( lattice->param.y0 < 0)
  {
    lattice->param.y0 = lattice->param.LY/2;
  }
  if( lattice->param.z0 < 0)
  {
    lattice->param.z0 = lattice->param.LZ/2;
  }

  if( lattice->param.r0 < 0)
  {
    if( lattice->param.LX < lattice->param.LY)
    {
      if( lattice->param.LX < lattice->param.LZ)
      {
        lattice->param.r0 = lattice->param.LX;
      }
      else
      {
        lattice->param.r0 = lattice->param.LZ;
      }
    }
    else
    {
      if( lattice->param.LY < lattice->param.LZ)
      {
        lattice->param.r0 = lattice->param.LY;
      }
      else
      {
        lattice->param.r0 = lattice->param.LZ;
      }
    }
  }
#else
  if( lattice->param.x0 < 0.)
  {
    lattice->param.x0 = lattice->param.LX/2.;
  }
  if( lattice->param.y0 < 0.)
  {
    lattice->param.y0 = lattice->param.LY/2.;
  }
  if( lattice->param.z0 < 0.)
  {
    lattice->param.z0 = lattice->param.LZ/2.;
  }

  if( lattice->param.r0 < 0.)
  {
    lattice->param.r0 = (lattice->param.LX+lattice->param.LY+lattice->param.LZ)/12.;
  }
#endif

  // Set default value for cut if it is negative.
  if( lattice->param.cut < 0.)
  {
    lattice->param.cut = 0.3;
  }

  // Set default values for x1, x2, y1, y2, z1, z2 if they are negative.
#if INTEGER_IC_BOUND
  if( lattice->param.x1 < 0)
  {
    if( lattice->param.rel_x1 < 0.)
    {
      lattice->param.x1 = lattice->param.LX/4;
    }
    else
    {
      lattice->param.x1 = (int) llround(lattice->param.rel_x1*(real)lattice->param.LX);
    }
  }
  if( lattice->param.x2 < 0)
  {
    if( lattice->param.rel_x2 < 0.)
    {
      lattice->param.x2 = lattice->param.LX/2;
    }
    else
    {
      lattice->param.x2 = (int) llround(lattice->param.rel_x2*(real)lattice->param.LX);
    }
  }
  if( lattice->param.y1 < 0)
  {
    if( lattice->param.rel_y1 < 0.)
    {
      lattice->param.y1 = lattice->param.LY/4;
    }
    else
    {
      lattice->param.y1 = (int) llround(lattice->param.rel_y1*(real)lattice->param.LY);
    }
  }
  if( lattice->param.y2 < 0)
  {
    if( lattice->param.rel_y2 < 0.)
    {
      lattice->param.y2 = lattice->param.LY/2;
    }
    else
    {
      lattice->param.y2 = (int) llround(lattice->param.rel_y2*(real)lattice->param.LY);
    }
  }
  if( lattice->param.z1 < 0)
  {
    if( lattice->param.rel_z1 < 0.)
    {
      lattice->param.z1 = lattice->param.LZ/4;
    }
    else
    {
      lattice->param.z1 = (int) llround(lattice->param.rel_z1*(real)lattice->param.LZ);
    }
  }
  if( lattice->param.z2 < 0)
  {
    if( lattice->param.rel_z2 < 0.)
    {
      lattice->param.z2 = lattice->param.LZ/2;
    }
    else
    {
      lattice->param.z2 = (int) llround(lattice->param.rel_z2*(real)lattice->param.LZ);
    }
  }
#else
  if( lattice->param.x1 < 0.)
  {
    if( lattice->param.rel_x1 < 0.)
    {
      lattice->param.x1 = lattice->param.LX/4.-1;
    }
    else
    {
      lattice->param.x1 = lattice->param.rel_x1*lattice->param.LX-1;
    }
  }
  if( lattice->param.x2 < 0.)
  {
    if( lattice->param.rel_x2 < 0.)
    {
      lattice->param.x2 = lattice->param.LX/2.-1;
    }
    else
    {
      lattice->param.x2 = lattice->param.rel_x2*lattice->param.LX-1;
    }
  }
  if( lattice->param.y1 < 0.)
  {
    if( lattice->param.rel_y1 < 0.)
    {
      lattice->param.y1 = lattice->param.LY/4.-1;
    }
    else
    {
      lattice->param.y1 = lattice->param.rel_y1*lattice->param.LY-1;
    }
  }
  if( lattice->param.y2 < 0.)
  {
    if( lattice->param.rel_y2 < 0.)
    {
      lattice->param.y2 = lattice->param.LY/2.-1;
    }
    else
    {
      lattice->param.y2 = lattice->param.rel_y2*lattice->param.LY-1;
    }
  }
  if( lattice->param.z1 < 0.)
  {
    if( lattice->param.rel_z1 < 0.)
    {
      lattice->param.z1 = lattice->param.LZ/4.-1;
    }
    else
    {
      lattice->param.z1 = lattice->param.rel_z1*lattice->param.LZ-1;
    }
  }
  if( lattice->param.z2 < 0.)
  {
    if( lattice->param.rel_z2 < 0.)
    {
      lattice->param.z2 = lattice->param.LZ/2.-1;
    }
    else
    {
      lattice->param.z2 = lattice->param.rel_z2*lattice->param.LZ-1;
    }
  }
#endif


  fclose(in);

#if VERBOSITY_LEVEL >= 1
  printf("%s %d >> read_params() -- Done reading file \"params.in\".\n",
      __FILE__,__LINE__);
#endif /* VERBOSITY_LEVEL >= 1 */

printf("rho_A = %f\n",lattice->param.rho_A);

} /* void read_params( struct lattice_struct *lattice) */

// void dump_params( struct lattice_struct *lattice)
//##############################################################################
//
// D U M P   P A R A M S
//
//  - Output the problem parameters to a file.
//
void dump_params( struct lattice_struct *lattice)
{
  FILE *o;
  char filename[1024];

  if( get_num_procs(lattice) > 1)
  {
    sprintf( filename, "./out/params_proc%04d.dat", get_proc_id( lattice));
  }
  else
  {
    sprintf( filename, "./out/params.dat");
  }

  if( !( o = fopen(filename,"w+")))
  {
    printf("%s %d >> ERROR: fopen(\"%s\",\"w+\") = NULL.  Bye, bye!\n",
        __FILE__, __LINE__, filename);
    process_exit(1);
  }

  fprintf( o, "LX                   %d\n", lattice->param.LX             );
  fprintf( o, "LY                   %d\n", lattice->param.LY             );
  fprintf( o, "LY                   %d\n", lattice->param.LZ             );
  fprintf( o, "length_scale         %d\n", lattice->param.length_scale   );
  fprintf( o, "NumNodes             %d\n", lattice->NumNodes             );
  fprintf( o, "NumFrames            %d\n", lattice->param.NumFrames      );
  fprintf( o, "FrameRate            %d\n", lattice->param.FrameRate      );
  fprintf( o, "NumTimeSteps         %d\n", lattice->NumTimeSteps         );
  fprintf( o, "tau[0]               %f\n", lattice->param.tau[0]         );
  fprintf( o, "gforce[0][0]         %f\n", lattice->param.gforce[0][0]   );
  fprintf( o, "gforce[0][1]         %f\n", lattice->param.gforce[0][1]   );
  fprintf( o, "gforce[0][2]         %f\n", lattice->param.gforce[0][2]   );

  fprintf( o, "end_grav[0]          %d\n", lattice->param.end_grav[0]    );
  if( NUM_FLUID_COMPONENTS==2)
  {
  fprintf( o, "tau[1]               %f\n", lattice->param.tau[1]         );
  fprintf( o, "gforce[1][0]         %f\n", lattice->param.gforce[1][0]   );
  fprintf( o, "gforce[1][1]         %f\n", lattice->param.gforce[1][1]   );
  fprintf( o, "gforce[1][2]         %f\n", lattice->param.gforce[1][2]   );

  fprintf( o, "end_grav[1]          %d\n", lattice->param.end_grav[1]    );
  }
  else if( NUM_FLUID_COMPONENTS==1)
  {
  fprintf( o, "tau[1]               %s\n", "--"                          );
  fprintf( o, "gforce[1][0]         %s\n", "--"                          );
  fprintf( o, "gforce[1][1]         %s\n", "--"                          );
  fprintf( o, "gforce[1][2]         %s\n", "--"                          );
  fprintf( o, "end_grav[1]          %s\n", "--"                          );
  }
  else
  {
    printf(
      "read_params() -- "
      "Unhandled case "
      "NUM_FLUID_COMPONENTS = %d .  "
      "Exiting!\n",
      NUM_FLUID_COMPONENTS);
    process_exit(1);
  }
  fprintf( o, "buoyancy             %d\n", lattice->param.buoyancy        );
  fprintf( o, "incompressible       %d\n", lattice->param.incompressible  );
  fprintf( o, "simple_diffusion     %d\n", lattice->param.simple_diffusion);
  fprintf( o, "rho_A[0]             %f\n", lattice->param.rho_A[0]       );
  fprintf( o, "rho_B[0]             %f\n", lattice->param.rho_B[0]       );
  if( NUM_FLUID_COMPONENTS==2)
  {
  fprintf( o, "rho_A[1]             %f\n", lattice->param.rho_A[1]       );
  fprintf( o, "rho_B[1]             %f\n", lattice->param.rho_B[1]       );
  }
  else if( NUM_FLUID_COMPONENTS==1)
  {
  fprintf( o, "rho_A[1]             %s\n", "--"                          );
  fprintf( o, "rho_B[1]             %s\n", "--"                          );
  }
  else
  {
    printf(
      "read_params() -- "
      "Unhandled case "
      "NUM_FLUID_COMPONENTS = %d .  "
      "Exiting!\n",
      NUM_FLUID_COMPONENTS);
    process_exit(1);
  }
  //fprintf( o, "rho_w                %f\n", lattice->param.rhow          );

#if INAMURO_SIGMA_COMPONENT
  fprintf( o, "rho_sigma            %f\n", lattice->param.rho_sigma      );
  fprintf( o, "rho_sigma_in         %f\n", lattice->param.rho_sigma_in   );
  fprintf( o, "rho_sigma_out        %f\n", lattice->param.rho_sigma_out  );
  fprintf( o, "u_sigma              %f\n", lattice->param.u_sigma        );
  fprintf( o, "u_sigma_in           %f\n", lattice->param.u_sigma_in     );
  fprintf( o, "u_sigma_out          %f\n", lattice->param.u_sigma_out    );
  fprintf( o, "sigma_start          %d\n", lattice->param.sigma_start    );
  fprintf( o, "sigma_stop           %d\n", lattice->param.sigma_stop     );
  fprintf( o, "sigma_btc_rate       %d\n", lattice->param.sigma_btc_rate );
  fprintf( o, "sigma_btc_spot       %d\n", lattice->param.sigma_btc_spot );
#else /* !( INAMURO_SIGMA_COMPONENT) */
  fprintf( o, "rho_sigma            %s\n", "--"                          );
  fprintf( o, "rho_sigma_in         %s\n", "--"                          );
  fprintf( o, "rho_sigma_out        %s\n", "--"                          );
  fprintf( o, "u_sigma              %s\n", "--"                          );
  fprintf( o, "u_sigma_in           %s\n", "--"                          );
  fprintf( o, "u_sigma_out          %s\n", "--"                          );
  fprintf( o, "sigma_start          %s\n", "--"                          );
  fprintf( o, "sigma_stop           %s\n", "--"                          );
  fprintf( o, "sigma_btc_rate       %s\n", "--"                          );
  fprintf( o, "sigma_btc_spot       %s\n", "--"                          );
#endif /* INAMURO_SIGMA_COMPONENT */
  fprintf( o, "GZL                  %d\n", lattice->param.GZL            );
  fprintf( o, "PressureBC           %d\n", lattice->param.PressureBC     );
  fprintf( o, "AllBoundaryPeriodic           %d\n", lattice->param.AllBoundaryPeriodic     );
  fprintf( o, "rho_in               %f\n", lattice->param.rho_in         );
  fprintf( o, "rho_out              %f\n", lattice->param.rho_out        );
  fprintf( o, "ux_in                %f\n", lattice->param.ux_in          );
  fprintf( o, "ux_out               %f\n", lattice->param.ux_out         );
  fprintf( o, "uy_in                %f\n", lattice->param.uy_in          );
  fprintf( o, "uy_out               %f\n", lattice->param.uy_out         );
  fprintf( o, "uz_in                %f\n", lattice->param.uz_in          );
  fprintf( o, "uz_out               %f\n", lattice->param.uz_out         );
  fprintf( o, "big_V0               %f\n", lattice->param.big_V0         );
  fprintf( o, "big_V0_solid[0]      %f\n", lattice->param.big_V0_solid[0]);
  if( NUM_FLUID_COMPONENTS==2)
  {
  fprintf( o, "big_V0_solid[1]      %f\n", lattice->param.big_V0_solid[1]);
  }
  else if( NUM_FLUID_COMPONENTS==1)
  {
  fprintf( o, "big_V0_solid[1]      %s\n", "--"                          );
  }
  else
  {
    printf(
      "read_params() -- "
      "Unhandled case "
      "NUM_FLUID_COMPONENTS = %d .  "
      "Exiting!\n",
      NUM_FLUID_COMPONENTS);
    process_exit(1);
  }
  fprintf( o, "periodic_x[0]        %d\n", lattice->periodic_x[0]        );
  fprintf( o, "periodic_y[0]        %d\n", lattice->periodic_y[0]        );
  if( NUM_FLUID_COMPONENTS==2)
  {
  fprintf( o, "periodic_x[1]        %d\n", lattice->periodic_x[1]        );
  fprintf( o, "periodic_y[1]        %d\n", lattice->periodic_y[1]        );
  }
  else if( NUM_FLUID_COMPONENTS==1)
  {
  fprintf( o, "periodic_x[1]        %s\n", "--"                          );
  fprintf( o, "periodic_y[1]        %s\n", "--"                          );
  }
  else
  {
    printf(
      "read_params() -- "
      "Unhandled case "
      "NUM_FLUID_COMPONENTS = %d .  "
      "Exiting!\n",
      NUM_FLUID_COMPONENTS);
    process_exit(1);
  }
  fprintf( o, "ns_flag              %d\n", lattice->param.ns_flag        );
  fprintf( o, "ns                   %f\n", lattice->param.ns             );
  fprintf( o, "ic_poisseuille       %d\n", lattice->param.ic_poisseuille );
  fprintf( o, "bc_poisseuille       %d\n", lattice->param.bc_poisseuille );
  fprintf( o, "bc_slip_north        %d\n", lattice->param.bc_slip_north  );
#if INAMURO_SIGMA_COMPONENT
  fprintf( o, "bc_sigma_slip        %d\n", lattice->param.bc_sigma_slip  );
#else /* !( INAMURO_SIGMA_COMPONENT) */
  fprintf( o, "bc_sigma_slip        %s\n", "--"                          );
#endif /* INAMURO_SIGMA_COMPONENT */
  fprintf( o, "pressure_t_in[0]     %d\n", lattice->param.pressure_t_in[0]  );
  fprintf( o, "pressure_b_in[0]     %d\n", lattice->param.pressure_b_in[0]  );
  fprintf( o, "pressure_t_out[0]    %d\n", lattice->param.pressure_t_out[0] );
  fprintf( o, "pressure_b_out[0]    %d\n", lattice->param.pressure_b_out[0] );
  fprintf( o, "velocity_t_in[0]     %d\n", lattice->param.velocity_t_in[0]  );
  fprintf( o, "velocity_b_in[0]     %d\n", lattice->param.velocity_b_in[0]  );
  fprintf( o, "velocity_t_out[0]    %d\n", lattice->param.velocity_t_out[0] );
  fprintf( o, "velocity_b_out[0]    %d\n", lattice->param.velocity_b_out[0] );
  fprintf( o, "pressure_n_in[0]     %d\n", lattice->param.pressure_n_in[0]  );
  fprintf( o, "pressure_s_in[0]     %d\n", lattice->param.pressure_s_in[0]  );
  fprintf( o, "pressure_n_out[0]    %d\n", lattice->param.pressure_n_out[0] );
  fprintf( o, "pressure_s_out[0]    %d\n", lattice->param.pressure_s_out[0] );
  fprintf( o, "velocity_n_in[0]     %d\n", lattice->param.velocity_n_in[0]  );
  fprintf( o, "velocity_s_in[0]     %d\n", lattice->param.velocity_s_in[0]  );
  fprintf( o, "velocity_n_out[0]    %d\n", lattice->param.velocity_n_out[0] );
  fprintf( o, "velocity_s_out[0]    %d\n", lattice->param.velocity_s_out[0] );
  fprintf( o, "pressure_e_in[0]     %d\n", lattice->param.pressure_e_in[0]  );
  fprintf( o, "pressure_w_in[0]     %d\n", lattice->param.pressure_w_in[0]  );
  fprintf( o, "pressure_e_out[0]    %d\n", lattice->param.pressure_e_out[0] );
  fprintf( o, "pressure_w_out[0]    %d\n", lattice->param.pressure_w_out[0] );
  fprintf( o, "velocity_e_in[0]     %d\n", lattice->param.velocity_e_in[0]  );
  fprintf( o, "velocity_w_in[0]     %d\n", lattice->param.velocity_w_in[0]  );
  fprintf( o, "velocity_e_out[0]    %d\n", lattice->param.velocity_e_out[0] );
  fprintf( o, "velocity_w_out[0]    %d\n", lattice->param.velocity_w_out[0] );

  if( NUM_FLUID_COMPONENTS==2)
  {

  fprintf( o, "pressure_t_in[0]     %d\n", lattice->param.pressure_t_in[1]  );
  fprintf( o, "pressure_b_in[0]     %d\n", lattice->param.pressure_b_in[1]  );
  fprintf( o, "pressure_t_out[0]    %d\n", lattice->param.pressure_t_out[1] );
  fprintf( o, "pressure_b_out[0]    %d\n", lattice->param.pressure_b_out[1] );
  fprintf( o, "velocity_t_in[0]     %d\n", lattice->param.velocity_t_in[1]  );
  fprintf( o, "velocity_b_in[0]     %d\n", lattice->param.velocity_b_in[1]  );
  fprintf( o, "velocity_t_out[0]    %d\n", lattice->param.velocity_t_out[1] );
  fprintf( o, "velocity_b_out[0]    %d\n", lattice->param.velocity_b_out[1] );
  fprintf( o, "pressure_n_in[1]     %d\n", lattice->param.pressure_n_in[1]  );
  fprintf( o, "pressure_s_in[1]     %d\n", lattice->param.pressure_s_in[1]  );
  fprintf( o, "pressure_n_out[1]    %d\n", lattice->param.pressure_n_out[1] );
  fprintf( o, "pressure_s_out[1]    %d\n", lattice->param.pressure_s_out[1] );
  fprintf( o, "velocity_n_in[1]     %d\n", lattice->param.velocity_n_in[1]  );
  fprintf( o, "velocity_s_in[1]     %d\n", lattice->param.velocity_s_in[1]  );
  fprintf( o, "velocity_n_out[1]    %d\n", lattice->param.velocity_n_out[1] );
  fprintf( o, "velocity_s_out[1]    %d\n", lattice->param.velocity_s_out[1] );
  fprintf( o, "pressure_e_in[1]     %d\n", lattice->param.pressure_e_in[1]  );
  fprintf( o, "pressure_w_in[1]     %d\n", lattice->param.pressure_w_in[1]  );
  fprintf( o, "pressure_e_out[1]    %d\n", lattice->param.pressure_e_out[1] );
  fprintf( o, "pressure_w_out[1]    %d\n", lattice->param.pressure_w_out[1] );
  fprintf( o, "velocity_e_in[1]     %d\n", lattice->param.velocity_e_in[1]  );
  fprintf( o, "velocity_w_in[1]     %d\n", lattice->param.velocity_w_in[1]  );
  fprintf( o, "velocity_e_out[1]    %d\n", lattice->param.velocity_e_out[1] );
  fprintf( o, "velocity_w_out[1]    %d\n", lattice->param.velocity_w_out[1] );
  }
  else if( NUM_FLUID_COMPONENTS==1)
  {
  fprintf( o, "pressure_t_in[1]     %s\n", "--"                             );
  fprintf( o, "pressure_b_in[1]     %s\n", "--"                             );
  fprintf( o, "pressure_t_out[1]    %s\n", "--"                             );
  fprintf( o, "pressure_b_out[1]    %s\n", "--"                             );
  fprintf( o, "velocity_t_in[1]     %s\n", "--"                             );
  fprintf( o, "velocity_b_in[1]     %s\n", "--"                             );
  fprintf( o, "velocity_t_out[1]    %s\n", "--"                             );
  fprintf( o, "velocity_b_out[1]    %s\n", "--"                             );
  fprintf( o, "pressure_n_in[1]     %s\n", "--"                             );
  fprintf( o, "pressure_s_in[1]     %s\n", "--"                             );
  fprintf( o, "pressure_n_out[1]    %s\n", "--"                             );
  fprintf( o, "pressure_s_out[1]    %s\n", "--"                             );
  fprintf( o, "velocity_n_in[1]     %s\n", "--"                             );
  fprintf( o, "velocity_s_in[1]     %s\n", "--"                             );
  fprintf( o, "velocity_n_out[1]    %s\n", "--"                             );
  fprintf( o, "velocity_s_out[1]    %s\n", "--"                             );
  fprintf( o, "pressure_e_in[1]     %s\n", "--"                             );
  fprintf( o, "pressure_w_in[1]     %s\n", "--"                             );
  fprintf( o, "pressure_e_out[1]    %s\n", "--"                             );
  fprintf( o, "pressure_w_out[1]    %s\n", "--"                             );
  fprintf( o, "velocity_e_in[1]     %s\n", "--"                             );
  fprintf( o, "velocity_w_in[1]     %s\n", "--"                             );
  fprintf( o, "velocity_e_out[1]    %s\n", "--"                             );
  fprintf( o, "velocity_w_out[1]    %s\n", "--"                             );
  }
  else
  {
    printf(
      "read_params() -- "
      "Unhandled case "
      "NUM_FLUID_COMPONENTS = %d .  "
      "Exiting!\n",
      NUM_FLUID_COMPONENTS);
    process_exit(1);
  }
  if( INAMURO_SIGMA_COMPONENT)
  {
  fprintf( o, "constcon_t_in        %d\n", lattice->param.constcon_t_in  );
  fprintf( o, "constcon_b_in        %d\n", lattice->param.constcon_b_in  );
  fprintf( o, "constcon_t_out       %d\n", lattice->param.constcon_t_out );
  fprintf( o, "constcon_b_out       %d\n", lattice->param.constcon_b_out );
  fprintf( o, "constflx_t_in        %d\n", lattice->param.constflx_t_in  );
  fprintf( o, "constflx_b_in        %d\n", lattice->param.constflx_b_in  );
  fprintf( o, "constflx_t_out       %d\n", lattice->param.constflx_t_out );
  fprintf( o, "constflx_b_out       %d\n", lattice->param.constflx_b_out );
  fprintf( o, "constcon_n_in        %d\n", lattice->param.constcon_n_in  );
  fprintf( o, "constcon_s_in        %d\n", lattice->param.constcon_s_in  );
  fprintf( o, "constcon_n_out       %d\n", lattice->param.constcon_n_out );
  fprintf( o, "constcon_s_out       %d\n", lattice->param.constcon_s_out );
  fprintf( o, "constflx_n_in        %d\n", lattice->param.constflx_n_in  );
  fprintf( o, "constflx_s_in        %d\n", lattice->param.constflx_s_in  );
  fprintf( o, "constflx_n_out       %d\n", lattice->param.constflx_n_out );
  fprintf( o, "constflx_s_out       %d\n", lattice->param.constflx_s_out );
  fprintf( o, "constcon_e_in        %d\n", lattice->param.constcon_e_in  );
  fprintf( o, "constcon_w_in        %d\n", lattice->param.constcon_w_in  );
  fprintf( o, "constcon_e_out       %d\n", lattice->param.constcon_e_out );
  fprintf( o, "constcon_w_out       %d\n", lattice->param.constcon_w_out );
  fprintf( o, "constflx_e_in        %d\n", lattice->param.constflx_e_in  );
  fprintf( o, "constflx_w_in        %d\n", lattice->param.constflx_w_in  );
  fprintf( o, "constflx_e_out       %d\n", lattice->param.constflx_e_out );
  fprintf( o, "constflx_w_out       %d\n", lattice->param.constflx_w_out );
  fprintf( o, "zeroconcgrad_t       %d\n", lattice->param.zeroconcgrad_t );
  fprintf( o, "zeroconcgrad_b       %d\n", lattice->param.zeroconcgrad_b );
  fprintf( o, "zeroconcgrad_n       %d\n", lattice->param.zeroconcgrad_n );
  fprintf( o, "zeroconcgrad_s       %d\n", lattice->param.zeroconcgrad_s );
  fprintf( o, "zeroconcgrad_e       %d\n", lattice->param.zeroconcgrad_e );
  fprintf( o, "zeroconcgrad_w       %d\n", lattice->param.zeroconcgrad_w );
  fprintf( o, "zeroconcgrad_full    %d\n", lattice->param.zeroconcgrad_full );
  }
  else
  {
  fprintf( o, "constcon_t_in        %s\n", "--"                             );
  fprintf( o, "constcon_b_in        %s\n", "--"                             );
  fprintf( o, "constcon_t_out       %s\n", "--"                             );
  fprintf( o, "constcon_b_out       %s\n", "--"                             );
  fprintf( o, "constflx_t_in        %s\n", "--"                             );
  fprintf( o, "constflx_b_in        %s\n", "--"                             );
  fprintf( o, "constflx_t_out       %s\n", "--"                             );
  fprintf( o, "constflx_b_out       %s\n", "--"                             );
  fprintf( o, "constcon_n_in        %s\n", "--"                             );
  fprintf( o, "constcon_s_in        %s\n", "--"                             );
  fprintf( o, "constcon_n_out       %s\n", "--"                             );
  fprintf( o, "constcon_s_out       %s\n", "--"                             );
  fprintf( o, "constflx_n_in        %s\n", "--"                             );
  fprintf( o, "constflx_s_in        %s\n", "--"                             );
  fprintf( o, "constflx_n_out       %s\n", "--"                             );
  fprintf( o, "constflx_s_out       %s\n", "--"                             );
  fprintf( o, "constcon_e_in        %s\n", "--"                             );
  fprintf( o, "constcon_w_in        %s\n", "--"                             );
  fprintf( o, "constcon_e_out       %s\n", "--"                             );
  fprintf( o, "constcon_w_out       %s\n", "--"                             );
  fprintf( o, "constflx_e_in        %s\n", "--"                             );
  fprintf( o, "constflx_w_in        %s\n", "--"                             );
  fprintf( o, "constflx_e_out       %s\n", "--"                             );
  fprintf( o, "constflx_w_out       %s\n", "--"                             );
  fprintf( o, "zeroconcgrad_t       %s\n", "--"                             );
  fprintf( o, "zeroconcgrad_b       %s\n", "--"                             );
  fprintf( o, "zeroconcgrad_n       %s\n", "--"                             );
  fprintf( o, "zeroconcgrad_s       %s\n", "--"                             );
  fprintf( o, "zeroconcgrad_e       %s\n", "--"                             );
  fprintf( o, "zeroconcgrad_w       %s\n", "--"                             );
  fprintf( o, "zeroconcgrad_full    %s\n", "--"                             );
  }

  fprintf( o, "plot_scale_dynamic   %d\n", lattice->param.plot_scale_dynamic);
  fprintf( o, "use_colormap         %d\n", lattice->param.use_colormap      );
  fprintf( o, "initial_condition    %d\n", lattice->param.initial_condition );
  fprintf( o, "x0                   %f\n", lattice->param.x0                );
  fprintf( o, "y0                   %f\n", lattice->param.y0                );
  fprintf( o, "z0                   %f\n", lattice->param.z0                );
  fprintf( o, "r0                   %f\n", lattice->param.r0                );
  fprintf( o, "cut                  %f\n", lattice->param.cut               );
  fprintf( o, "x1                   %f\n", lattice->param.x1                );
  fprintf( o, "x2                   %f\n", lattice->param.x2                );
  fprintf( o, "y1                   %f\n", lattice->param.y1                );
  fprintf( o, "y2                   %f\n", lattice->param.y2                );
  fprintf( o, "z1                   %f\n", lattice->param.z1                );
  fprintf( o, "z2                   %f\n", lattice->param.z2                );
  fprintf( o, "rel_x1               %f\n", lattice->param.rel_x1            );
  fprintf( o, "rel_x2               %f\n", lattice->param.rel_x2            );
  fprintf( o, "rel_y1               %f\n", lattice->param.rel_y1            );
  fprintf( o, "rel_y2               %f\n", lattice->param.rel_y2            );
  fprintf( o, "rel_z1               %f\n", lattice->param.rel_z1            );
  fprintf( o, "rel_z2               %f\n", lattice->param.rel_z2            );
  fprintf( o, "dump_rho             %d\n", lattice->param.dump_rho          );
  fprintf( o, "dump_u               %d\n", lattice->param.dump_u            );
  fprintf( o, "dump_force           %d\n", lattice->param.dump_force        );
  fprintf( o, "dump_vor             %d\n", lattice->param.dump_vor          );

  fclose(o);

#if VERBOSITY_LEVEL >= 1
  printf("%s %d >> dump_params() -- Wrote file \"params.dat\".\n",
      __FILE__, __LINE__);
#endif /* VERBOSITY_LEVEL >= 1 */

} /* void dump_params( struct lattice_struct *lattice) */
