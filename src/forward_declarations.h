#ifndef FORWARD_DECLARATIONS_H
#define FORWARD_DECLARATIONS_H
//##############################################################################
//
// forward_declarations.h
//
//  - Forward declarations of routines for lb_prime.
//

void dump_frame_summary( lattice_ptr lattice);
void dump_macro_vars( lattice_ptr lattice, int time);
void collide( lattice_ptr lattice);
void stream( lattice_ptr lattice);
void stream_collide_stream( lattice_ptr lattice);
void stream_save( lattice_ptr lattice);
void compute_macro_vars( lattice_ptr lattice);
void compute_feq( lattice_ptr lattice);
void compute_single_feq( lattice_ptr lattice, int n, int subs, real *feq);
void compute_a_feq(  real *feq, real rho, real u_x, real u_y, real u_z);
void compute_max_u( lattice_ptr lattice, int subs);
void compute_min_u( lattice_ptr lattice, int subs);
void compute_ave_u( lattice_ptr lattice, int subs);
void compute_max_rho( lattice_ptr lattice, int subs);
void compute_min_rho( lattice_ptr lattice, int subs);
void compute_ave_rho( lattice_ptr lattice, int subs);
void compute_flux( lattice_ptr lattice, int subs);
void compute_ave_ueq( lattice_ptr lattice, real *max_u);
void process_matrix( lattice_ptr lattice, int **matrix);
void process_bcs( lattice_ptr lattice, int subs);
void read_params( lattice_ptr lattice, const char *infile);
void construct_lattice( lattice_ptr *lattice, int argc, char **argv);
void init_problem( lattice_ptr lattice);
void destruct_lattice( lattice_ptr lattice);
void dump_north_pointing_pdfs(
       lattice_ptr lattice,
       const int subs,
       const int z_slice,
       char *comment_str,
       const int which_pdf);
void dump_south_pointing_pdfs(
       lattice_ptr lattice,
       const int subs,
       const int z_slice,
       char *comment_str,
       const int which_pdf);
int domain_is_too_big_to_display( lattice_ptr lattice);
int domain_is_not_too_big_to_display( lattice_ptr lattice);
void display_warning_about_contrived_data( lattice_ptr lattice);
void dump_macro_vars( lattice_ptr lattice, int time);
void dump_pdf( lattice_ptr lattice, int time);
void dump_lattice_info( lattice_ptr lattice);
#if DO_NOT_STORE_SOLIDS
void dump_node_info( lattice_ptr lattice);
#endif /* DO_NOT_STORE_SOLIDS */
void dump_checkpoint( lattice_ptr lattice, int time, char *fn);
void read_checkpoint( lattice_ptr lattice);
void read_solids( lattice_ptr lattice, char *filename);
void write_raw(
       lattice_ptr lattice,
       real *a,
       int     stride,
       real  a_max,
       real  a_min,
       char   *filename);
void write_raw_u(
       lattice_ptr lattice,
       real *ux,
       real *uy,
       real *uz,
       int     stride,
       real  ux_max,
       real  ux_min,
       real  uy_max,
       real  uy_min,
       real  uz_max,
       real  uz_min,
       char   *filename);
void write_dat(
       lattice_ptr lattice,
       real *a,
       int     stride,
       real  a_max,
       real  a_min,
       char   *filename);
void write_plt(
       lattice_ptr lattice,
       real *a,
       real *b,
       char   *filename, char   *filename2);
void write_plt_uvw(
       lattice_ptr lattice,
       real *a,
       real *b,
       char   *filename);
void write_plt_single(
       lattice_ptr lattice,
       real *a,
       char   *filename);
void write_rho_txt( lattice_ptr lattice);
void write_u_txt( lattice_ptr lattice);
void write_f_txt( lattice_ptr lattice);
void slice( lattice_ptr lattice);
void private_slice( lattice_ptr lattice, int i0, int j0, int i1, int j1);
#if NON_LOCAL_FORCES
void compute_phase_force( lattice_ptr lattice, int subs);
void compute_fluid_fluid_force( lattice_ptr lattice);
void compute_real_fluid_solid_force( lattice_ptr lattice);
void compute_single_fluid_solid_force( lattice_ptr lattice, int subs);
void dump_forces( lattice_ptr lattice);
#endif /* NON_LOCAL_FORCES */
void count_colormap( int *num_colors);
void allocate_colormap( real ***colormap, int num_colors);
void read_colormap( real **colormap, int num_colors);
void deallocate_colormap( real ***colormap, int num_colors);
void get_color(
       real **colormap, int num_colors,
       real c, char *r, char *g, char *b);
#if MANAGE_BODY_FORCE
inline void manage_body_force( lattice_ptr lattice);
#endif /* MANAGE_BODY_FORCE */
#if INAMURO_SIGMA_COMPONENT && STORE_BTC
void sigma_stuff( lattice_ptr lattice);
#endif /* INAMURO_SIGMA_COMPONENT && STORE_BTC */

void bmp_read_header( FILE *in, struct bitmap_info_header *bmih);
void bmp_read_entry(
  FILE *in,
  struct bitmap_info_header bmih,
  char *r, char *g, char *b);

void report_open( report_ptr report, char *name);
void report_integer_entry(
       report_ptr report, char *label, int value, char *units);
void report_ratio_entry(
       report_ptr report, char *label, real num, real den, char *units);
void report_entry( report_ptr report, char *entry_left, char *entry_right);
void report_close( report_ptr report);
void report_partition( report_ptr report);

int get_sizeof_lattice_structure( lattice_ptr lattice);
int get_sizeof_lattice(           lattice_ptr lattice);
int get_num_active_nodes(         lattice_ptr lattice);

#endif /* FORWARD_DECLARATIONS_H */
