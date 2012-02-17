//##############################################################################
//
// lattice.h
//
//  - Lattice data structures for lbgpu_prime.
//

#ifndef LATTICE_H
#define LATTICE_H

// struct pdf_struct
//
//  - Structure to hold the particle distribution functions.
//
#if SAVE_MEMO
struct pdf_struct
{
  real ftemp[ /*NUM_DIRS*/ 19];
};
#else
struct pdf_struct
{
  real feq[   /*NUM_DIRS*/ 19];
  real f[     /*NUM_DIRS*/ 19];
  real ftemp[ /*NUM_DIRS*/ 19];

};
#endif
typedef struct pdf_struct *pdf_ptr;

// struct macro_vars_struct
//
//  - Structure to hold the macroscopic variables.
//
struct macro_vars_struct
{
  real rho;
  real u[   /*NUM_DIMS*/ 3];

};

#if STORE_UEQ
// struct ueq_struct
//
//  - Structure to hold the composite velocity used to compute feq.
//
struct ueq_struct
{
  real u[   /*NUM_DIMS*/ 3];

};
#endif /* STORE_UEQ */

// struct solids_struct
//
//  - Structure to hold solids information.
//
struct solids_struct
{
  unsigned char is_solid;

};
typedef struct solids_struct *solids_ptr;

#if DO_NOT_STORE_SOLIDS
// struct node_struct
//
//  - Structure to hold node information.
//
struct node_struct
{
  int i, j; // Lattice coordinates of node.

  int n; // Node index.

  int nn[18]; // Indices of neighboring nodes.

};
typedef struct node_struct *node_ptr;
#endif /* DO_NOT_STORE_SOLIDS */

#if NON_LOCAL_FORCES
struct force_struct
{
  real force[  /*NUM_DIMS*/3];
  real sforce[ /*NUM_DIMS*/3];
};
#endif /* NON_LOCAL_FORCES */
#if POROUS_MEDIA
struct ns_struct
{
  real ns;
};
#endif /* POROUS_MEDIA */

struct process_struct
{
  int id;
  int num_procs;
#if PARALLEL
  int g_LX;
  int g_LY;
  int g_LZ;
  int g_NumNodes;

  int g_SX, g_EX;
  int g_SY, g_EY;
  int g_SZ, g_EZ;
  int g_StartNode;

  real *pos_dir_pdf_to_send;
  real *pos_dir_pdf_to_recv;
  real *neg_dir_pdf_to_send;
  real *neg_dir_pdf_to_recv;

  real *pos_dir_rho_to_send;
  real *pos_dir_rho_to_recv;
  real *neg_dir_rho_to_send;
  real *neg_dir_rho_to_recv;
  unsigned char *pos_dir_solid_to_send;
  unsigned char *pos_dir_solid_to_recv;
  unsigned char *neg_dir_solid_to_send;
  unsigned char *neg_dir_solid_to_recv;


  MPI_Request send_req_0;
  MPI_Request recv_req_0;
  MPI_Request send_req_1;
  MPI_Request recv_req_1;
//rev Huang
  MPI_Request send_req_2;
  MPI_Request recv_req_2;
  MPI_Request send_req_3;
  MPI_Request recv_req_3;
  MPI_Request send_req_4;
  MPI_Request recv_req_4;
  MPI_Request send_req_5;
  MPI_Request recv_req_5;

  MPI_Status mpi_status;
  int mpierr;

#endif
};
typedef struct process_struct *process_ptr;

// struct param_struct
//
// Structure to hold parameters that specify the configuration of the problem
// to solve.
//
// Values for these parameters are read from file "./in/params.in". (Note that
// there are also some compile-time configuration options in file
// "./in/flags.h".)
//
// Each parameter is preceeded by documentation consisting of three kinds of
// fields: "Param", "Type", and "Comments". There may be multiple "Param"
// fields for sets of related parameters, e.g. LX and LY. There is only one
// occurance of the "Type" field and "Comments" field for each (set of)
// parameter(s). The "Param" field(s) and the "Type" field are each one line.
// The "Comments" field can be multiple lines.
//
// Distinct fluids are referred to both as components and as substances.
// Normally in conversation we talk about fluid components, but variables are
// indexed by "subs". The term "component" can be used to refer to a fluid
// component or a solute component (c.f. Inamuro) as well as the more broadly
// conventional mathematical usage of a vector component. The reader should be
// wary of context, therefore.
//
// In the "Param" field(s) of the documentation preceeding each parameter, the
// parameters that are arrays are stated with typical and/or suggestive index
// names between square brackets. For example, tau[subs] is written to suggest
// the fact that tau contains an array of values, one for each substance or
// fluid component, where "subs" is a typical index name for indexing fluid
// components. The commentary section ("Comments" field) for each parameter is
// explicit about the meaning of the indices.
//
// A common mistake made when configuring a problem is to specify in the file
// "params.in" an integer number for a type real parameter or vice versa.
// This can cause mystifying symptoms. Be careful to include a decimal point in
// all type real parameter values in "params.in" regardless of whether there
// is a fractional part to the value. Likewise, do not include decimal points
// in type integer parameter values. Remember! If things have been going
// smoothly and then suddenly go bad in weird ways (sorry I can't recollect any
// specific examples) be sure to check this issue in "params.in". Note that
// file "params.dat" is output by the code mostly as a sanity check: if values
// in "params.dat" do not match values in "params.in", the real versus
// integer parameter type issue may be the reason.
//
// Here are some useful Unix/Linux command lines for processing this
// documentation (the "params.txt" file):
//
//  >> grep Param params.txt
//
//  >> grep Comment params.txt
//
//  >> TODO
//
//  >> TODO
//
// Here is the command line that generates the "params.txt" file from the
// "data_structures.h" file:
//
//  >> cat ./src/data_structures.h | \
//     sed -n '/^\/\/ struct param_struct/,/\/\* struct param_struct/p' | \
//     grep "\/\/" | \
//     sed 's/[ ]*\/\/ //; s/[ ]*\/\/$//' > ./doc/params.txt
//
// If you are developing under Windows, I recommend Cygwin
//
//   http://www.cygwin.com
//
// or something comparable. The Unix/Linux command line is indispensable. See
//
//   http://www.cs.usfca.edu/~parrt/course/601/lectures/unix.util.html
//   http://www.student.northpark.edu/pemente/sed/sed1line.txt
//
// for example.
//
// -----------------------------------------------------------------------------
//
struct param_struct
{
  // Param: LX
  // Param: LY
  // Param: LZ
  // Type: int
  // Comments: LX is the number of lattice nodes in x-direction. LY is the
  // number of lattice nodes in y-direction. LZ is the number of lattice nodes
  // in z-direction.
  int    LX,
         LY,
         LZ;

  // Param BX, BY, BZ for X, Y and Z block sizes
  // These must be chosen by the user subject to a set of rules.
  //
  int    BX,
         BY,
         BZ;

 // int    NGPU;

  // Param: length_scale
  // Type: int
  // Comments: length_scale is a characteristic length. This is used for
  // computing Re.
  //
  int    length_scale;

  // Param: NumFrames
  // Type: int
  // Comments: NumFrames is the number of frames to output.
  //
  int    NumFrames;

  // Param: FrameRate
  // Type: int
  // Comments: FrameRate is the number of timesteps per frame.
  //
  int    FrameRate;

  // Param: NumTimeSteps
  // Type: int
  // Comments: NumTimeSteps is the number of time steps computed as
  // NumFrames*FrameRate.
  //
  int    NumTimeSteps;

  // Param: tau[subs]
  // Type: real*
  // Comments: tau is the relaxation time(s). If NUM_FLUID_COMPONENTS==1, then
  // tau[0] is the relaxation time. If NUM_FLUID_COMPONENTS==2, then tau[0] is
  // the relaxation time for component 0 and tau[1] is the relaxation time for
  // component 1.
  //
  real tau[ NUM_FLUID_COMPONENTS];

  // Param: gforce[subs][dir]
  // Type: real**
  // Comments: gforce[subs][dir] is the gravitational (body) force(s). If
  // NUM_FLUID_COMPONENTS==1, then gforce[0][dir] holds the gravitational force.
  // If NUM_FLUID_COMPONENTS==2, then gforce[0][dir] holds the gravitational
  // force for component 0 and gforce[1][dir] holds the gravitational force for
  // component 1. In either case, gforce[subs][0] is the gravitational force
  // in the x-direction and gforce[subs][1] is the gravitational force in the
  // y-direction, gforce[subs][2] is the gravitational force in the
  // z-direction.

  //
  real gforce[   NUM_FLUID_COMPONENTS][3];

  // Param: end_grav[subs]
  // Type: int
  // Comments: end_grav[subs] is the timestep at which to end the gravitational
  // force. If NUM_FLUID_COMPONENTS==1, then end_grav[0] holds the ending
  // timestep. If NUM_FLUID_COMPONENTS==2, then end_grav[0] holds the ending
  // timestep for fluid compenent 0 and end_grav[1] holds the ending timestep
  // for fluid component 1. This mechanism is only enforced when the
  // MANAGE_BODY_FORCE flag is turned on in "flags.h". It is implemented as a
  // call to the "manage_body_force()" function at each time step, so it is a
  // special feature that should generally be turned off. In the future this
  // mechanism should mature.
  //
  int    end_grav[ NUM_FLUID_COMPONENTS];

  // Param: buoyancy
  // Type: int
  // Comments: The "buoyancy" flag toggles buoyancy effects due to the density
  // of solute. This flag is only relevant when INAMURO_SIGMA_COMPONENT is
  // turned on in "flags.h".
  //
  int    buoyancy;

  // Param: incompressible
  // Type: int
  // Comments: The "incompressible" flag toggles an incompressible lattice
  // Boltzmann method (c.f. TODO: cite). The usual method that we experiment
  // with is the common "weakly" compressible BGK method. For some things
  // (e.g., TODO) an incompressible model is desired.
  //
  int    incompressible;

  // Param: simple_diffusion
  // Type: int Comments: The simple_diffusion flag toggles a simplified
  // diffusion computation. It simplifies the computation of feq for the solute
  // component. It is based on the observation that the process of diffusion
  // requires less symmetry than flow (c.f. TODO: cite).
  //
  int    simple_diffusion;

  // Param: rho_A[subs]
  // Param: rho_B[subs]
  // Type: real*
  // Comments: rho_A[subs] and rho_B[subs] are the initial density values. The
  // usage of these values varies with respect to the chosen initial condition
  // (see the initial_condition parameter below). If NUM_FLUID_COMPONENTS==1,
  // then rho_A[0] and/or rho_B[0] are used to initialize the fluid density.
  // If NUM_FLUID_COMPONENTS==2, then rho_A[0] and/or rho_B[0] are used to
  // initialize density of fluid component 0 and rho_A[1] and/or rho_B[1] are
  // used to initialize density of fluid component 1. All the values are not
  // necessarily required for all the initial conditions. See the documentation
  // related to the initial_condition parameter for specifics.
  //
  real rho_A[    NUM_FLUID_COMPONENTS];
  real rho_B[    NUM_FLUID_COMPONENTS];
  real rhow ;

#if INAMURO_SIGMA_COMPONENT
  // Initial/boundary condition values for solute component.
  real rho_sigma;
  real rho_sigma_in;
  real rho_sigma_out;
  real u_sigma;
  real u_sigma_in;
  real u_sigma_out;

  // Param: sigma_start;
  // Type: int    sigma_start;
  // Comments: Timestep at which solute is activated. If sigma_start == 0, solute will be activated from the beginning. The reason for this option is to make it easy to allow the fluid to equilibrate before introducing solute into the flow.
  //
  int    sigma_start;

  // Param: sigma_stop;
  // Type: int    sigma_stop;
  // Comments: Timestep at which solute is deactivated. If sigma_stop < 0 then the solute will never stop. WARNING: if 0 < sigma_stop < sigma_start, solute will never start. This option is useful, for instance, if it is desired to observe how a solute washes out of the domain.
  //
  int    sigma_stop;

  // Param: sigma_btc_rate;
  // Type: int    sigma_btc_rate;
  // Comments: Rate at which to accumulate temporal breakthrough curve. If sigma_btc_rate <= 0, no btc is stored.
  //
  int    sigma_btc_rate;

  // Param: sigma_btc_spot;
  // Type: int    sigma_btc_spot;
  // Comments: The spot in the domain (as distance from inflow) for measuring btc. If sigma_btc_spot is greater than the length of the domain or <0 then the btc will be measured at the outflow.
  //
  int    sigma_btc_spot;

  // Param: *sigma_btc;
  // Type: real *sigma_btc;
  // Comments: Pointer to array for storing break through curve. If sigma_btc_rate <= 0, this will remain unallocated.
  //
  real *sigma_btc;

#endif /* INAMURO_SIGMA_COMPONENT */

  // Inflow and outflow conditions for density and velocity.
  real rho_in,   rho_out;
  real ux_in,    ux_out;
  real uy_in,    uy_out;
  real uz_in,    uz_out;

  // Param: big_V0;
  // Type: real big_V0;
  // Comments: Interaction parameter.
  //
  real big_V0;

  // Param: big_V0_solid[NUM_FLUID_COMPONENTS];
  // Type: real big_V0_solid[NUM_FLUID_COMPONENTS];
  // Comments: Adsorption parameters.
  //
  real big_V0_solid[NUM_FLUID_COMPONENTS];

  // Param: ns_flag;
  // Type: int ns_flag;
  // Comments: Flag for how to use the solid density parameter.
  // If ns_flag = 0, then initialize the domain uniformly with ns.
  // If ns_flag = 1, then read ns values from the ns<LX>x<LY>.bmp file.
  // If ns_flag = 2, then read ns mask from ns<LX>x<LY>.bmp -- that is,
  // initialize the domain with ns where there are black pixels in the
  // ns<LX>x<LY>.bmp file and 0 otherwise.
  //
  int    ns_flag;

  // Param: ns;
  // Type: real ns;
  // Comments: Solid density parameter for porous media.
  //
  real ns;

  // Param: ic_poisseuille;
  // Param: bc_poisseuille;
  // Type: int    ic_poisseuille;
  // Comments: Flags for poisseuille flow initial condition (ic) and boundary condition (bc).
  //
  int    ic_poisseuille;
  int    bc_poisseuille;

  // Param: bc_slip_north;
  // Type: int    bc_slip_north;
  // Comments: Toggle slip condition on north wall.
  //
  int    bc_slip_north;
#if INAMURO_SIGMA_COMPONENT
  // Param: bc_sigma_slip;
  // Type: int    bc_sigma_slip;
  // Comments: Slip BC for solute on side walls. Will this make a difference on Taylor dispersion results?  NOTE: Only use this for flow through a channel. Slip BC is not implemented for arbitrary geometries.
  //
  int    bc_sigma_slip;
#endif /* INAMURO_SIGMA_COMPONENT */
  int  GZL;
  int  PressureBC;
  int AllBoundaryPeriodic;
  // Flags for pressure and velocity boundaries.
  int    pressure_t_in[  NUM_FLUID_COMPONENTS];
  int    pressure_b_in[  NUM_FLUID_COMPONENTS];
  int    pressure_t_out[ NUM_FLUID_COMPONENTS];
  int    pressure_b_out[ NUM_FLUID_COMPONENTS];
  int    velocity_t_in[  NUM_FLUID_COMPONENTS];
  int    velocity_b_in[  NUM_FLUID_COMPONENTS];
  int    velocity_t_out[ NUM_FLUID_COMPONENTS];
  int    velocity_b_out[ NUM_FLUID_COMPONENTS];
  int    pressure_n_in[  NUM_FLUID_COMPONENTS];
  int    pressure_s_in[  NUM_FLUID_COMPONENTS];
  int    pressure_n_out[ NUM_FLUID_COMPONENTS];
  int    pressure_s_out[ NUM_FLUID_COMPONENTS];
  int    velocity_n_in[  NUM_FLUID_COMPONENTS];
  int    velocity_s_in[  NUM_FLUID_COMPONENTS];
  int    velocity_n_out[ NUM_FLUID_COMPONENTS];
  int    velocity_s_out[ NUM_FLUID_COMPONENTS];
  int    pressure_e_in[  NUM_FLUID_COMPONENTS];
  int    pressure_w_in[  NUM_FLUID_COMPONENTS];
  int    pressure_e_out[ NUM_FLUID_COMPONENTS];
  int    pressure_w_out[ NUM_FLUID_COMPONENTS];
  int    velocity_e_in[  NUM_FLUID_COMPONENTS];
  int    velocity_w_in[  NUM_FLUID_COMPONENTS];
  int    velocity_e_out[ NUM_FLUID_COMPONENTS];
  int    velocity_w_out[ NUM_FLUID_COMPONENTS];

  // Flags for (inamuro) concentration boundaries.
  int    constcon_t_in;
  int    constcon_b_in;
  int    constcon_t_out;
  int    constcon_b_out;
  int    constflx_t_in;
  int    constflx_b_in;
  int    constflx_t_out;
  int    constflx_b_out;
  int    constcon_n_in;
  int    constcon_s_in;
  int    constcon_n_out;
  int    constcon_s_out;
  int    constflx_n_in;
  int    constflx_s_in;
  int    constflx_n_out;
  int    constflx_s_out;
  int    constcon_e_in;
  int    constcon_w_in;
  int    constcon_e_out;
  int    constcon_w_out;
  int    constflx_e_in;
  int    constflx_w_in;
  int    constflx_e_out;
  int    constflx_w_out;

  // Flags for zero concentration gradient boundaries.
  // These differ from inamuro's concentration flux boundaries
  // in that they allow advective flux of the concentration
  // across the boundary.
  int    zeroconcgrad_t;
  int    zeroconcgrad_b;
  int    zeroconcgrad_n;
  int    zeroconcgrad_s;
  int    zeroconcgrad_e;
  int    zeroconcgrad_w;

  // Param: zeroconcgrad_full;
  // Type: int    zeroconcgrad_full;
  // Comments: The zero concentration gradient can be computed by copying (from the adjacent interior neighbor) either just the unknown distribution functions after streaming or all the distributions after streaming. This flag toggles between those methods (where "full" denotes the latter method).
  //
  int    zeroconcgrad_full;

  // Param: use_colormap;
  // Type: int    use_colormap;
  // Comments: See the colormap routines in lbio.c to learn about colormap support.
  //
  int    use_colormap;

  // Param: plot_scale_dynamic;
  // Type: int    plot_scale_dynamic;
  // Comments: Toggle dynamic scaling of the plots. If this is on, color values to be plotted will be dynamically scaled to range between the minimum and maximum values of the data for the current plot. If this is off, color values will be scaled only by the initial value for that quantity (e.g. rho_sigma, rho_A[subs], etc...). NOTE: This is only implemented in rho2bmp() as of 10/15/2004.
  //
  int    plot_scale_dynamic;

  // Param: initial_condition;
  // Type: int    initial_condition;
  // Comments: See bottom of flags.h for listing/description of initial conditions.
  //
  int    initial_condition;

  // Center and radius of a circle. Used in the initial conditions for
  // initializing a bubble in the domain, for instance.

#if INTEGER_IC_BOUND
  int x0;
  int y0;
  int z0;
  int r0;
#else
  real x0;
  real y0;
  real z0;
  real r0;
#endif

  // Param: cut;
  // Type: real cut;
  // Comments: Cut-off value. Used in the initial conditions for generating static in the two-component case. Sometimes the static needs to be biased toward one component in order for it to evolve stabily.
  //
  real cut;

  // Corner coordinates for a rectangle. Used by initial conditions for
  // initializing a square-shaped region in the domain.

#if INTEGER_IC_BOUND
  int x1, x2, z1;
  int y1, y2, z2;
#else
  real x1, x2, z1;
  real y1, y2, z2;
#endif

  // Factors for setting coordinates relative to the domain size. If
  // the absolute coordinates are negative, use these values to set
  // the coordinates, e.g. x1 = rel_x1*LX. If these rel_* values are
  // also negative, then revert to some default absolute coordinates.
  real rel_x1, rel_x2, rel_z1;
  real rel_y1, rel_y2, rel_z2;

  // Param: dump_rho;
  // Type: int    dump_rho;
  // Comments: Flags to toggle output of the macroscopic density.
  //
  int    dump_rho;

  // Param: dump_u;
  // Type: int    dump_u;
  // Comments: Flags to toggle output of the macroscopic velocity.
  //
  int    dump_u;

  // Param: dump_force;
  // Type: int    dump_force;
  // Comments: Flags to toggle output of the interaction and adsorption forces.
  //
  int    dump_force;

  // Param: dump_vor;
  // Type: int    dump_vor;
  // Comments: Flags to toggle output of vorticity.
  //
  int    dump_vor;

}; /* struct param_struct */

struct vars_struct
{
  real* f_memblock;
  real** f1d;
  real* macrovars_memblock;
  real** macrovars1d;
};

struct attrib_struct
{
  real min_macrovars[NUM_FLUID_COMPONENTS][/*MAX_DIMENSIONS+1*/4]
     , max_macrovars[NUM_FLUID_COMPONENTS][/*MAX_DIMENSIONS+1*/4]
     , ave_macrovars[NUM_FLUID_COMPONENTS][/*MAX_DIMENSIONS+1*/4];
  real          flux[NUM_FLUID_COMPONENTS][/*MAX_DIMENSIONS+1*/4];
};


// struct bcs_in_struct
//
//  - Structure to hold input values for boundary conditions.
//
//  - These are activated by setting the corresponding flag to 2 in the
//    params.in file. There must be a file with the corresponding name and
//    suffix 'dat' in the 'in' folder, e.g., pressure_n_in0.in. The values
//    in the input file will be used for the boundary conditions, successive
//    values on successive time steps. If there are more timesteps than input
//    values, the input values will be cycled. This is particularly useful for
//    imposing temporally periodic boundary conditions, e.g., to simulate tidal
//    periods.
//
struct bcs_in_struct
{
  real*  pressure_n_in0;
  int      num_pressure_n_in0;

  real*  pressure_s_in0;
  int      num_pressure_s_in0;
#if 0
  double*  pressure_s_in0;
  double* pressure_n_out0;
  double* pressure_s_out0;
  double*  velocity_n_in0;
  double*  velocity_s_in0;
  double* velocity_n_out0;
  double* velocity_s_out0;
  double*  pressure_e_in0;
  double*  pressure_w_in0;
  double* pressure_e_out0;
  double* pressure_w_out0;
  double*  velocity_e_in0;
  double*  velocity_w_in0;
  double* velocity_e_out0;
  double* velocity_w_out0;
  double*  pressure_n_in1;
  double*  pressure_s_in1;
  double* pressure_n_out1;
  double* pressure_s_out1;
  double*  velocity_n_in1;
  double*  velocity_s_in1;
  double* velocity_n_out1;
  double* velocity_s_out1;
  double*  pressure_e_in1;
  double*  pressure_w_in1;
  double* pressure_e_out1;
  double* pressure_w_out1;
  double*  velocity_e_in1;
  double*  velocity_w_in1;
  double* velocity_e_out1;
  double* velocity_w_out1;
  double*  constcon_n_in;
  double*  constcon_s_in;
  double* constcon_n_out;
  double* constcon_s_out;
  double*  constflx_n_in;
  double*  constflx_s_in;
  double* constflx_n_out;
  double* constflx_s_out;
  double*  constcon_e_in;
  double*  constcon_w_in;
  double* constcon_e_out;
  double* constcon_w_out;
  double*  constflx_e_in;
  double*  constflx_w_in;
  double* constflx_e_out;
  double* constflx_w_out;
#endif
};

typedef struct bcs_in_struct* bcs_in_ptr;




// struct lattice_struct
//
//  - Structure with all the lattice information.
//
//  - Contains arrays of the previously defined structures.
//
struct lattice_struct
{
  int    NumSubs;
  int    NumDims;
  int    NumVelDirs;
  int    NumBoundDirs;
  int    NumUnboundDirs;
  int    EndBoundSize;
  int    NumNodes;
  int    NumTimeSteps;
  int    time;
  int    frame;
  real tic;
  real toc;
  int    periodic_x[ NUM_FLUID_COMPONENTS];
  int    periodic_y[ NUM_FLUID_COMPONENTS];

  int    SizeBTC; // Number of BTC measurements to store.
  int    FlowDir; // Direction {1,2} ==> {Horiz,Vert} of flow.

  struct param_struct      param;
  struct attrib_struct     attrib;
  struct process_struct    process;

  struct vars_struct* vars;
  unsigned char* solids_memblock;
  real* ns_memblock;

  // For the time varying north and south pressure boundaries
  struct bcs_in_struct      bcs_in[     NUM_FLUID_COMPONENTS];

#if 0 // TODO
#if STORE_UEQ
  struct ueq_struct        *ueq;
#endif /* STORE_UEQ */
#endif

};
typedef struct lattice_struct *lattice_ptr;

struct report_struct
{
  FILE *file;
  char name[1024];
};
typedef struct report_struct *report_ptr;

struct bitmap_file_header
{
  // 1 2 bfType 19778 must always be set to 'BM' to declare that this is
  // a .bmp-file.
  char bfType[2];
  // 3 4 bfSize ?? specifies the size of the file in bytes.
  char bfSize[4];
  // 7 2 bfReserved1 0 must always be set to zero.
  char bfReserved1[2];
  // 9 2 bfReserved2 0 must always be set to zero.
  char bfReserved2[2];
  // 11 4 bfOffBits 1078 specifies the offset from the beginning of the
  // file to the bitmap data.
  char bfOffBits[4];

};

struct bitmap_info_header
{
  // 15 4 biSize 40 specifies the size of the BITMAPINFOHEADER structure,
  // in bytes.
  char biSize[4];
  // 19 4 biWidth 100 specifies the width of the image, in pixels.
  char biWidth[4];
  // 23 4 biHeight 100 specifies the height of the image, in pixels.
  char biHeight[4];
  // 27 2 biPlanes 1 specifies the number of planes of the target device,
  // must be set to zero. [DT: Should be set to one, right? Not zero.]
  char biPlanes[2];
  // 29 2 biBitCount 8 specifies the number of bits per pixel.
  char biBitCount[2];
  // 31 4 biCompression 0 Specifies the type of compression, usually set
  // to zero (no compression).
  char biCompression[4];
  // 35 4 biSizeImage 0 specifies the size of the image data, in bytes.
  // If there is no compression, it is valid to set this member to zero.
  char biSizeImage[4];
  // 39 4 biXPelsPerMeter 0 specifies the the horizontal pixels per meter
  // on the designated targer device, usually set to zero.
  char biXPelsPerMeter[4];
  // 43 4 biYPelsPerMeter 0 specifies the the vertical pixels per meter
  // on the designated targer device, usually set to zero.
  char biYPelsPerMeter[4];
  // 47 4 biClrUsed 0 specifies the number of colors used in the bitmap,
  // if set to zero the number of colors is calculated using the biBitCount
  // member.
  char biClrUsed[4];
  // 51 4 biClrImportant 0 specifies the number of color that are
  // 'important' for the bitmap, if set to zero, all colors are important.
  char biClrImportant[4];

};

struct rgb_quad
{
  // 1 1 rgbBlue - specifies the blue part of the color.
  char Blue;
  // 2 1 rgbGreen - specifies the green part of the color.
  char Green;
  // 3 1 rgbRed - specifies the red part of the color.
  char Red;
  // 4 1 rgbReserved - must always be set to zero.
  char Reserved;
};


//##############################################################################
//
// Forward declarations for lattice data accessors.
//

int get_LX( lattice_ptr lattice);
int get_LY( lattice_ptr lattice);
int get_LZ( lattice_ptr lattice);
int get_BX( lattice_ptr lattice);
int get_BY( lattice_ptr lattice);
int get_BZ( lattice_ptr lattice);
int get_ni( lattice_ptr lattice);
int get_nj( lattice_ptr lattice);
int get_nk( lattice_ptr lattice);
int* get_ni_ptr( lattice_ptr lattice);
int* get_nj_ptr( lattice_ptr lattice);
int* get_nk_ptr( lattice_ptr lattice);
void set_LX( lattice_ptr lattice, const int arg_LX);
void set_LY( lattice_ptr lattice, const int arg_LY);
void set_LZ( lattice_ptr lattice, const int arg_LZ);

int get_FrameRate( lattice_ptr lattice);
int get_NumFrames( lattice_ptr lattice);
int get_Frame( lattice_ptr lattice);

void set_NumNodes( lattice_ptr lattice);
int get_NumNodes( lattice_ptr lattice);
int* get_NumNodes_ptr( lattice_ptr lattice);

void set_NumDims( lattice_ptr lattice, int value);
int get_NumDims( lattice_ptr lattice);
int* get_NumDims_ptr( lattice_ptr lattice);

void set_NumVelDirs( lattice_ptr lattice, int value);
int get_NumVelDirs( lattice_ptr lattice);
int* get_NumVelDirs_ptr( lattice_ptr lattice);

void set_NumBoundDirs( lattice_ptr lattice, int value);
int get_NumBoundDirs( lattice_ptr lattice);

void set_NumUnboundDirs( lattice_ptr lattice, int value);
int get_NumUnboundDirs( lattice_ptr lattice);

void set_EndBoundSize( lattice_ptr lattice, int value);
int get_EndBoundSize( lattice_ptr lattice);
int* get_EndBoundSize_ptr( lattice_ptr lattice);

int get_NumSubs( lattice_ptr lattice);
int* get_NumSubs_ptr( lattice_ptr lattice);


int is_solid( lattice_ptr lattice, const int n);
void set_is_solid(
  lattice_ptr lattice
, const int n
, const unsigned char val);
int is_not_solid( lattice_ptr lattice, const int n);

real get_rho( lattice_ptr lattice, const int subs, const int n);
real* get_rho_ptr( lattice_ptr lattice, const int subs, const int n);
real get_ux( lattice_ptr lattice, const int subs, const int n);
real get_uy( lattice_ptr lattice, const int subs, const int n);
real get_uz( lattice_ptr lattice, const int subs, const int n);
real* get_ux_ptr( lattice_ptr lattice, const int subs, const int n);
real* get_uy_ptr( lattice_ptr lattice, const int subs, const int n);
real* get_uz_ptr( lattice_ptr lattice, const int subs, const int n);

#endif
