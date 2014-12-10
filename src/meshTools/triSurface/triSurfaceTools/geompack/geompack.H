double d_epsilon ( void );
double d_max ( double x, double y );
double d_min ( double x, double y );
void d2vec_part_quick_a ( int n, double a[], int *l, int *r );
void d2vec_permute ( int n, double a[], int p[] );
int *d2vec_sort_heap_index_a ( int n, double a[] );
void d2vec_sort_quick_a ( int n, double a[] );
int diaedg ( double x0, double y0, double x1, double y1, double x2, double y2,
  double x3, double y3 );
void dmat_transpose_print ( int m, int n, double a[], const char *title );
void dmat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, const char *title );
void dmat_uniform ( int m, int n, double b, double c, int *seed, double r[] );
int dtris2 ( int point_num, double point_xy[], int *tri_num,
  int tri_vert[], int tri_nabe[] );
bool dvec_eq ( int n, double a1[], double a2[] );
bool dvec_gt ( int n, double a1[], double a2[] );
bool dvec_lt ( int n, double a1[], double a2[] );
void dvec_print ( int n, double a[], const char *title );
void dvec_swap ( int n, double a1[], double a2[] );
int i_max ( int i1, int i2 );
int i_min ( int i1, int i2 );
int i_modp ( int i, int j );
int i_sign ( int i );
int i_wrap ( int ival, int ilo, int ihi );
void imat_transpose_print ( int m, int n, int a[], const char *title );
void imat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, const char *title );
void ivec_heap_d ( int n, int a[] );
int *ivec_indicator ( int n );
void ivec_sort_heap_a ( int n, int a[] );
void ivec_sorted_unique ( int n, int a[], int *nuniq );
int lrline ( double xu, double yu, double xv1, double yv1, double xv2,
  double yv2, double dv );
bool perm_check ( int n, int p[] );
void perm_inv ( int n, int p[] );
int *points_delaunay_naive_2d ( int n, double p[], int *ntri );
int s_len_trim ( const char *s );
int swapec ( int i, int *top, int *btri, int *bedg, int point_num,
  double point_xy[], int tri_num, int tri_vert[], int tri_nabe[],
  int stack[] );
void timestamp ( void );
char *timestring ( void );
double *triangle_circumcenter_2d ( double t[] );
bool triangulation_plot_eps ( const char *file_out_name,
  int g_num, double g_xy[], int tri_num, int nod_tri[] );
void triangulation_print ( int point_num, double xc[], int tri_num,
  int tri_vert[], int tri_nabe[] );
void vbedg ( double x, double y, int point_num, double point_xy[], int tri_num,
  int tri_vert[], int tri_nabe[], int *ltri, int *ledg, int *rtri, int *redg );
