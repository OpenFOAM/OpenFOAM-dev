# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "geompack.H"

//******************************************************************************

double d_epsilon ( void )

//******************************************************************************
//
//  Purpose:
//
//    D_EPSILON returns the round off unit for double precision arithmetic.
//
//  Discussion:
//
//    D_EPSILON is a number R which is a power of 2 with the property that,
//    to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double D_EPSILON, the floating point round-off unit.
//
{
    double r = 1.0;

    while (1.0 < 1.0 + r)
    {
        r = r/2.0;
    }

    return 2.0*r;
}


//*********************************************************************

double d_max ( double x, double y )

//*********************************************************************
//
//  Purpose:
//
//    D_MAX returns the maximum of two real values.
//
//  Modified:
//
//    10 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double D_MAX, the maximum of X and Y.
//
{
  if ( y < x )
  {
    return x;
  }
  else
  {
    return y;
  }
}
//*********************************************************************

double d_min ( double x, double y )

//*********************************************************************
//
//  Purpose:
//
//    D_MIN returns the minimum of two real values.
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double D_MIN, the minimum of X and Y.
//
{
  if ( y < x )
  {
    return y;
  }
  else
  {
    return x;
  }
}
//******************************************************************************

void d2vec_part_quick_a ( int n, double a[], int *l, int *r )

//******************************************************************************
//
//  Purpose:
//
//    D2VEC_PART_QUICK_A reorders an R2 vector as part of a quick sort.
//
//  Discussion:
//
//    The routine reorders the entries of A.  Using A(1:2,1) as a
//    key, all entries of A that are less than or equal to the key will
//    precede the key, which precedes all entries that are greater than the key.
//
//  Example:
//
//    Input:
//
//      N = 8
//
//      A = ( (2,4), (8,8), (6,2), (0,2), (10,6), (10,0), (0,6), (4,8) )
//
//    Output:
//
//      L = 2, R = 4
//
//      A = ( (0,2), (0,6), (2,4), (8,8), (6,2), (10,6), (10,0), (4,8) )
//             -----------          ----------------------------------
//             LEFT          KEY    RIGHT
//
//  Modified:
//
//    01 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of A.
//
//    Input/output, double A[N*2].  On input, the array to be checked.
//    On output, A has been reordered as described above.
//
//    Output, int *L, *R, the indices of A that define the three segments.
//    Let KEY = the input value of A(1:2,1).  Then
//    I <= L                 A(1:2,I) < KEY;
//         L < I < R         A(1:2,I) = KEY;
//                 R <= I    A(1:2,I) > KEY.
//
{
  int i;
  int j;
  double key[2];
  int ll;
  int m;
  int rr;

  if ( n < 1 )
  {
    cout << "\n";
    cout << "D2VEC_PART_QUICK_A - Fatal error!\n";
    cout << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    *l = 0;
    *r = 2;
    return;
  }

  key[0] = a[2*0+0];
  key[1] = a[2*0+1];
  m = 1;
//
//  The elements of unknown size have indices between L+1 and R-1.
//
  ll = 1;
  rr = n + 1;

  for ( i = 2; i <= n; i++ )
  {
    if ( dvec_gt ( 2, a+2*ll, key ) )
    {
      rr = rr - 1;
      dvec_swap ( 2, a+2*(rr-1), a+2*ll );
    }
    else if ( dvec_eq ( 2, a+2*ll, key ) )
    {
      m = m + 1;
      dvec_swap ( 2, a+2*(m-1), a+2*ll );
      ll = ll + 1;
    }
    else if ( dvec_lt ( 2, a+2*ll, key ) )
    {
      ll = ll + 1;
    }

  }
//
//  Now shift small elements to the left, and KEY elements to center.
//
  for ( i = 0; i < ll - m; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      a[2*i+j] = a[2*(i+m)+j];
    }
  }

  ll = ll - m;

  for ( i = ll; i < ll+m; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      a[2*i+j] = key[j];
    }
  }

  *l = ll;
  *r = rr;

  return;
}
//******************************************************************************

void d2vec_permute ( int n, double a[], int p[] )

//******************************************************************************
//
//  Purpose:
//
//    D2VEC_PERMUTE permutes an R2 vector in place.
//
//  Discussion:
//
//    This routine permutes an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   2,    4,    5,    1,    3 )
//      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
//          (11.0, 22.0, 33.0, 44.0, 55.0 )
//
//    Output:
//
//      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
//             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
//
//  Modified:
//
//    19 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input/output, double A[2*N], the array to be permuted.
//
//    Input, int P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array.  P must be a legal permutation
//    of the integers from 1 to N, otherwise the algorithm will
//    fail catastrophically.
//
{
  double a_temp[2];
  int i;
  int iget;
  int iput;
  int istart;

  if ( !perm_check ( n, p ) )
  {
    cout << "\n";
    cout << "D2VEC_PERMUTE - Fatal error!\n";
    cout << "  The input array does not represent\n";
    cout << "  a proper permutation.\n";
    exit ( 1 );
  }
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = -p[istart-1];
      continue;
    }
    else
    {
      a_temp[0] = a[0+(istart-1)*2];
      a_temp[1] = a[1+(istart-1)*2];
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = -p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cout << "\n";
          cout << "D2VEC_PERMUTE - Fatal error!\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[0+(iput-1)*2] = a_temp[0];
          a[1+(iput-1)*2] = a_temp[1];
          break;
        }
        a[0+(iput-1)*2] = a[0+(iget-1)*2];
        a[1+(iput-1)*2] = a[1+(iget-1)*2];
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = -p[i];
  }

  return;
}
//******************************************************************************

int *d2vec_sort_heap_index_a ( int n, double a[] )

//******************************************************************************
//
//  Purpose:
//
//    D2VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of
//    an R2 vector.
//
//  Discussion:
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      A(1:2,INDX(I)), I = 1 to N is sorted,
//
//    or explicitly, by the call
//
//      call D2VEC_PERMUTE ( N, A, INDX )
//
//    after which A(1:2,I), I = 1 to N is sorted.
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[2*N], an array to be index-sorted.
//
//    Output, int D2VEC_SORT_HEAP_INDEX_A[N], the sort index.  The
//    I-th element of the sorted array is A(0:1,D2VEC_SORT_HEAP_INDEX_A(I-1)).
//
{
  double aval[2];
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  if ( n == 1 )
  {
    indx = new int[1];
    indx[0] = 1;
    return indx;
  }

  indx = ivec_indicator ( n );

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval[0] = a[0+(indxt-1)*2];
      aval[1] = a[1+(indxt-1)*2];
    }
    else
    {
      indxt = indx[ir-1];
      aval[0] = a[0+(indxt-1)*2];
      aval[1] = a[1+(indxt-1)*2];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }

    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if (   a[0+(indx[j-1]-1)*2] <  a[0+(indx[j]-1)*2] ||
             ( a[0+(indx[j-1]-1)*2] == a[0+(indx[j]-1)*2] &&
               a[1+(indx[j-1]-1)*2] <  a[1+(indx[j]-1)*2] ) )
        {
          j = j + 1;
        }
      }

      if (   aval[0] <  a[0+(indx[j-1]-1)*2] ||
           ( aval[0] == a[0+(indx[j-1]-1)*2] &&
             aval[1] <  a[1+(indx[j-1]-1)*2] ) )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return indx;
}
//*****************************************************************************

void d2vec_sort_quick_a ( int n, double a[] )

//*****************************************************************************
//
//  Purpose:
//
//    D2VEC_SORT_QUICK_A ascending sorts an R2 vector using quick sort.
//
//  Discussion:
//
//    The data structure is a set of N pairs of real numbers.
//    These values are stored in a one dimensional array, by pairs.
//
//  Modified:
//
//    01 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, double A[N*2].
//    On input, the array to be sorted.
//    On output, the array has been sorted.
//
{
# define LEVEL_MAX 25

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( n < 1 )
  {
    cout << "\n";
    cout << "D2VEC_SORT_QUICK_A - Fatal error!\n";
    cout << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[level-1] = n + 1;
  base = 1;
  n_segment = n;

  while ( 0 < n_segment )
  {
//
//  Partition the segment.
//
    d2vec_part_quick_a ( n_segment, a+2*(base-1)+0, &l_segment, &r_segment );
//
//  If the left segment has more than one element, we need to partition it.
//
    if ( 1 < l_segment )
    {
      if ( LEVEL_MAX < level )
      {
        cout << "\n";
        cout << "D2VEC_SORT_QUICK_A - Fatal error!\n";
        cout << "  Exceeding recursion maximum of " << LEVEL_MAX << "\n";
        exit ( 1 );
      }

      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }
//
//  The left segment and the middle segment are sorted.
//  Must the right segment be partitioned?
//
    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }
//
//  Otherwise, we back up a level if there is an earlier one.
//
    else
    {
      for ( ; ; )
      {
        if ( level <= 1 )
        {
          n_segment = 0;
          break;
        }

        base = rsave[level-1];
        n_segment = rsave[level-2] - rsave[level-1];
        level = level - 1;

        if ( 0 < n_segment )
        {
          break;
        }
      }
    }
  }
  return;
# undef LEVEL_MAX
}
//******************************************************************************

int diaedg ( double x0, double y0, double x1, double y1, double x2, double y2,
  double x3, double y3 )

//******************************************************************************
//
//  Purpose:
//
//    DIAEDG chooses a diagonal edge.
//
//  Discussion:
//
//    The routine determines whether 0--2 or 1--3 is the diagonal edge
//    that should be chosen, based on the circumcircle criterion, where
//    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
//    quadrilateral in counterclockwise order.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    Barry Joe,
//    Department of Computing Science,
//    University of Alberta,
//    Edmonton, Alberta, Canada  T6G 2H1
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double X0, Y0, X1, Y1, X2, Y2, X3, Y3, the coordinates of the
//    vertices of a quadrilateral, given in counter clockwise order.
//
//    Output, int DIAEDG, chooses a diagonal:
//    +1, if diagonal edge 02 is chosen;
//    -1, if diagonal edge 13 is chosen;
//     0, if the four vertices are cocircular.
//
{
  double ca;
  double cb;
  double dx10;
  double dx12;
  double dx30;
  double dx32;
  double dy10;
  double dy12;
  double dy30;
  double dy32;
  double s;
  double tol;
  double tola;
  double tolb;
  int value;

  tol = 100.0 * d_epsilon ( );

  dx10 = x1 - x0;
  dy10 = y1 - y0;
  dx12 = x1 - x2;
  dy12 = y1 - y2;
  dx30 = x3 - x0;
  dy30 = y3 - y0;
  dx32 = x3 - x2;
  dy32 = y3 - y2;

  tola = tol * d_max ( fabs ( dx10 ),
               d_max ( fabs ( dy10 ),
               d_max ( fabs ( dx30 ), fabs ( dy30 ) ) ) );

  tolb = tol * d_max ( fabs ( dx12 ),
               d_max ( fabs ( dy12 ),
               d_max ( fabs ( dx32 ), fabs ( dy32 ) ) ) );

  ca = dx10 * dx30 + dy10 * dy30;
  cb = dx12 * dx32 + dy12 * dy32;

  if ( tola < ca && tolb < cb )
  {
    value = -1;
  }
  else if ( ca < -tola && cb < -tolb )
  {
    value = 1;
  }
  else
  {
    tola = d_max ( tola, tolb );
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb
      + ( dx32 * dy12 - dx12 * dy32 ) * ca;

    if ( tola < s )
    {
      value = -1;
    }
    else if ( s < -tola )
    {
      value = 1;
    }
    else
    {
      value = 0;
    }

  }

  return value;
}
//******************************************************************************

void dmat_transpose_print ( int m, int n, double a[], const char *title )

//******************************************************************************
//
//  Purpose:
//
//    DMAT_TRANSPOSE_PRINT prints a real matrix, transposed.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, const char *TITLE, an optional title.
//
{
  dmat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//******************************************************************************

void dmat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, const char *title )

//******************************************************************************
//
//  Purpose:
//
//    DMAT_TRANSPOSE_PRINT_SOME prints some of a real matrix, transposed.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, const char *TITLE, an optional title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  for ( i2lo = i_max ( ilo, 1 ); i2lo <= i_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i_min ( i2hi, m );
    i2hi = i_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i_max ( jlo, 1 );
    j2hi = i_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << " ";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }
  cout << "\n";

  return;
# undef INCX
}
//******************************************************************************

void dmat_uniform ( int m, int n, double b, double c, int *seed, double r[] )

//******************************************************************************
//
//  Purpose:
//
//    DMAT_UNIFORM fills a double precision array with scaled
//    pseudorandom values.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Modified:
//
//    30 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double B, C, the limits of the pseudorandom values.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and D_UNIFORM will be 0.  On output, SEED has
//    been updated.
//
//    Output, double DMAT_UNIFORM[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int j;
  int k;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
      r[i+j*m] = b + ( c - b ) * double( *seed ) * 4.656612875E-10;
    }
  }

  return;
}
//******************************************************************************

int dtris2 ( int point_num, double point_xy[], int *tri_num,
  int tri_vert[], int tri_nabe[] )

//******************************************************************************
//
//  Purpose:
//
//    DTRIS2 constructs a Delaunay triangulation of 2D vertices.
//
//  Discussion:
//
//    The routine constructs the Delaunay triangulation of a set of 2D vertices
//    using an incremental approach and diagonal edge swaps.  Vertices are
//    first sorted in lexicographically increasing (X,Y) order, and
//    then are inserted one at a time from outside the convex hull.
//
//  Modified:
//
//    15 January 2004
//
//  Author:
//
//    Barry Joe,
//    Department of Computing Science,
//    University of Alberta,
//    Edmonton, Alberta, Canada  T6G 2H1
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of vertices.
//
//    Input/output, double POINT_XY[POINT_NUM*2], the coordinates of
//    the vertices.  On output, the vertices have been sorted into
//    dictionary order.
//
//    Output, int *TRI_NUM, the number of triangles in the triangulation;
//    TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the number
//    of boundary vertices.
//
//    Output, int TRI_VERT[TRI_NUM*3], the nodes that make up each triangle.
//    The elements are indices of POINT_XY.  The vertices of the triangles are
//    in counter clockwise order.
//
//    Output, int TRI_NABE[TRI_NUM*3], the triangle neighbor list.
//    Positive elements are indices of TIL; negative elements are used for links
//    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
//    where I, J = triangle, edge index; TRI_NABE[I,J] refers to
//    the neighbor along edge from vertex J to J+1 (mod 3).
//
//    Output, int DTRIS2, is 0 for no error.
{
  double cmax;
  int e;
  int error;
  int i;
  int *indx;
  int j;
  int k;
  int l;
  int ledg;
  int lr;
  int ltri;
  int m;
  int m1;
  int m2;
  int n;
  int redg;
  int rtri;
  int *stack;
  int t;
  double tol;
  int top;

  stack = new int[point_num];

  tol = 100.0 * d_epsilon ( );
//
//  Sort the vertices by increasing (x,y).
//
  indx = d2vec_sort_heap_index_a ( point_num, point_xy );

  d2vec_permute ( point_num, point_xy, indx );
//
//  Make sure that the data points are "reasonably" distinct.
//
  m1 = 1;

  for ( i = 2; i <= point_num; i++ )
  {
    m = m1;
    m1 = i;

    k = -1;

    for ( j = 0; j <= 1; j++ )
    {
      cmax = d_max ( fabs ( point_xy[2*(m-1)+j] ),
                     fabs ( point_xy[2*(m1-1)+j] ) );

      if ( tol * ( cmax + 1.0 )
           < fabs ( point_xy[2*(m-1)+j] - point_xy[2*(m1-1)+j] ) )
      {
        k = j;
        break;
      }

    }

    if ( k == -1 )
    {
      cout << "\n";
      cout << "DTRIS2 - Fatal error!\n";
      cout << "  Fails for point number I = " << i << "\n";
      cout << "  M =  " << m  << "\n";
      cout << "  M1 = " << m1 << "\n";
      cout << "  X,Y(M)  = " << point_xy[2*(m-1)+0] << "  "
                             << point_xy[2*(m-1)+1] << "\n";
      cout << "  X,Y(M1) = " << point_xy[2*(m1-1)+0] << "  "
                             << point_xy[2*(m1-1)+1] << "\n";
      delete [] stack;
      return 224;
    }

  }
//
//  Starting from points M1 and M2, search for a third point M that
//  makes a "healthy" triangle (M1,M2,M)
//
  m1 = 1;
  m2 = 2;
  j = 3;

  for ( ; ; )
  {
    if ( point_num < j )
    {
      cout << "\n";
      cout << "DTRIS2 - Fatal error!\n";
      delete [] stack;
      return 225;
    }

    m = j;

    lr = lrline ( point_xy[2*(m-1)+0], point_xy[2*(m-1)+1],
      point_xy[2*(m1-1)+0], point_xy[2*(m1-1)+1],
      point_xy[2*(m2-1)+0], point_xy[2*(m2-1)+1], 0.0 );

    if ( lr != 0 )
    {
      break;
    }

    j = j + 1;

  }
//
//  Set up the triangle information for (M1,M2,M), and for any other
//  triangles you created because points were collinear with M1, M2.
//
  *tri_num = j - 2;

  if ( lr == -1 )
  {
    tri_vert[3*0+0] = m1;
    tri_vert[3*0+1] = m2;
    tri_vert[3*0+2] = m;
    tri_nabe[3*0+2] = -3;

    for ( i = 2; i <= *tri_num; i++ )
    {
      m1 = m2;
      m2 = i+1;
      tri_vert[3*(i-1)+0] = m1;
      tri_vert[3*(i-1)+1] = m2;
      tri_vert[3*(i-1)+2] = m;
      tri_nabe[3*(i-2)+0] = -3 * i;
      tri_nabe[3*(i-2)+1] = i;
      tri_nabe[3*(i-1)+2] = i - 1;

    }

    tri_nabe[3*(*tri_num-1)+0] = -3 * (*tri_num) - 1;
    tri_nabe[3*(*tri_num-1)+1] = -5;
    ledg = 2;
    ltri = *tri_num;
  }
  else
  {
    tri_vert[3*0+0] = m2;
    tri_vert[3*0+1] = m1;
    tri_vert[3*0+2] = m;
    tri_nabe[3*0+0] = -4;

    for ( i = 2; i <= *tri_num; i++ )
    {
      m1 = m2;
      m2 = i+1;
      tri_vert[3*(i-1)+0] = m2;
      tri_vert[3*(i-1)+1] = m1;
      tri_vert[3*(i-1)+2] = m;
      tri_nabe[3*(i-2)+2] = i;
      tri_nabe[3*(i-1)+0] = -3 * i - 3;
      tri_nabe[3*(i-1)+1] = i - 1;
    }

    tri_nabe[3*(*tri_num-1)+2] = -3 * (*tri_num);
    tri_nabe[3*0+1] = -3 * (*tri_num) - 2;
    ledg = 2;
    ltri = 1;
  }
//
//  Insert the vertices one at a time from outside the convex hull,
//  determine visible boundary edges, and apply diagonal edge swaps until
//  Delaunay triangulation of vertices (so far) is obtained.
//
  top = 0;

  for ( i = j+1; i <= point_num; i++ )
  {
    m = i;
    m1 = tri_vert[3*(ltri-1)+ledg-1];

    if ( ledg <= 2 )
    {
      m2 = tri_vert[3*(ltri-1)+ledg];
    }
    else
    {
      m2 = tri_vert[3*(ltri-1)+0];
    }

    lr = lrline ( point_xy[2*(m-1)+0], point_xy[2*(m-1)+1],
      point_xy[2*(m1-1)+0], point_xy[2*(m1-1)+1],
      point_xy[2*(m2-1)+0], point_xy[2*(m2-1)+1], 0.0 );

    if ( 0 < lr )
    {
      rtri = ltri;
      redg = ledg;
      ltri = 0;
    }
    else
    {
      l = -tri_nabe[3*(ltri-1)+ledg-1];
      rtri = l / 3;
      redg = (l % 3) + 1;
    }

    vbedg ( point_xy[2*(m-1)+0], point_xy[2*(m-1)+1], point_num,
      point_xy, *tri_num, tri_vert, tri_nabe, &ltri, &ledg, &rtri, &redg );

    n = *tri_num + 1;
    l = -tri_nabe[3*(ltri-1)+ledg-1];

    for ( ; ; )
    {
      t = l / 3;
      e = ( l % 3 ) + 1;
      l = -tri_nabe[3*(t-1)+e-1];
      m2 = tri_vert[3*(t-1)+e-1];

      if ( e <= 2 )
      {
        m1 = tri_vert[3*(t-1)+e];
      }
      else
      {
        m1 = tri_vert[3*(t-1)+0];
      }

      *tri_num = *tri_num + 1;
      tri_nabe[3*(t-1)+e-1] = *tri_num;
      tri_vert[3*(*tri_num-1)+0] = m1;
      tri_vert[3*(*tri_num-1)+1] = m2;
      tri_vert[3*(*tri_num-1)+2] = m;
      tri_nabe[3*(*tri_num-1)+0] = t;
      tri_nabe[3*(*tri_num-1)+1] = *tri_num - 1;
      tri_nabe[3*(*tri_num-1)+2] = *tri_num + 1;
      top = top + 1;

      if ( point_num < top )
      {
        cout << "\n";
        cout << "DTRIS2 - Fatal error!\n";
        cout << "  Stack overflow.\n";
        delete [] stack;
        return 8;
      }

      stack[top-1] = *tri_num;

      if ( t == rtri && e == redg )
      {
        break;
      }

    }

    tri_nabe[3*(ltri-1)+ledg-1] = -3 * n - 1;
    tri_nabe[3*(n-1)+1] = -3 * (*tri_num) - 2;
    tri_nabe[3*(*tri_num-1)+2] = -l;
    ltri = n;
    ledg = 2;

    error = swapec ( m, &top, &ltri, &ledg, point_num, point_xy, *tri_num,
      tri_vert, tri_nabe, stack );

    if ( error != 0 )
    {
      cout << "\n";
      cout << "DTRIS2 - Fatal error!\n";
      cout << "  Error return from SWAPEC.\n";
      delete [] stack;
      return error;
    }

  }
//
//  Now account for the sorting that we did.
//
  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < *tri_num; j++ )
    {
      tri_vert[i+j*3] = indx [ tri_vert[i+j*3] - 1 ];
    }
  }

  perm_inv ( point_num, indx );

  d2vec_permute ( point_num, point_xy, indx );

  delete [] indx;
  delete [] stack;

  return 0;
}
//******************************************************************************

bool dvec_eq ( int n, double a1[], double a2[] )

//******************************************************************************
//
//  Purpose:
//
//    DVEC_EQ is true if every pair of entries in two vectors is equal.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], two vectors to compare.
//
//    Output, bool DVEC_EQ.
//    DVEC_EQ is TRUE if every pair of elements A1(I) and A2(I) are equal,
//    and FALSE otherwise.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] != a2[i] )
    {
      return false;
    }
  }
  return true;

}
//******************************************************************************

bool dvec_gt ( int n, double a1[], double a2[] )

//******************************************************************************
//
//  Purpose:
//
//    DVEC_GT == ( A1 > A2 ) for real vectors.
//
//  Discussion:
//
//    The comparison is lexicographic.
//
//    A1 > A2  <=>                              A1(1) > A2(1) or
//                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
//                 ...
//                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, double A1[N], A2[N], the vectors to be compared.
//
//    Output, bool DVEC_GT, is TRUE if and only if A1 > A2.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {

    if ( a2[i] < a1[i] )
    {
       return true;
    }
    else if ( a1[i] < a2[i] )
    {
      return false;
    }

  }

  return false;
}
//******************************************************************************

bool dvec_lt ( int n, double a1[], double a2[] )

//******************************************************************************
//
//  Purpose:
//
//    DVEC_LT == ( A1 < A2 ) for real vectors.
//
//  Discussion:
//
//    The comparison is lexicographic.
//
//    A1 < A2  <=>                              A1(1) < A2(1) or
//                 ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
//                 ...
//                 ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, double A1[N], A2[N], the vectors to be compared.
//
//    Output, bool DVEC_LT, is TRUE if and only if A1 < A2.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] < a2[i] )
    {
      return true;
    }
    else if ( a2[i] < a1[i] )
    {
      return false;
    }

  }

  return false;
}
//********************************************************************

void dvec_print ( int n, double a[], const char *title )

//********************************************************************
//
//  Purpose:
//
//    DVEC_PRINT prints a double precision vector.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, const char *TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = 0; i <= n-1; i++ )
  {
    cout << setw(6)  << i + 1 << "  "
         << setw(14) << a[i]  << "\n";
  }

  return;
}
//******************************************************************************

void dvec_swap ( int n, double a1[], double a2[] )

//******************************************************************************
//
//  Purpose:
//
//    DVEC_SWAP swaps the entries of two real vectors.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the arrays.
//
//    Input/output, double A1[N], A2[N], the vectors to swap.
//
{
  int i;
  double temp;

  for ( i = 0; i < n; i++ )
  {
    temp  = a1[i];
    a1[i] = a2[i];
    a2[i] = temp;
  }

  return;
}
//****************************************************************************

int i_max ( int i1, int i2 )

//****************************************************************************
//
//  Purpose:
//
//    I_MAX returns the maximum of two integers.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I_MAX, the larger of I1 and I2.
//
{
  if ( i2 < i1 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************

int i_min ( int i1, int i2 )

//****************************************************************************
//
//  Purpose:
//
//    I_MIN returns the smaller of two integers.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I_MIN, the smaller of I1 and I2.
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//*********************************************************************

int i_modp ( int i, int j )

//*********************************************************************
//
//  Purpose:
//
//    I_MODP returns the nonnegative remainder of integer division.
//
//  Formula:
//
//    If
//      NREM = I_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//  Comments:
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I_MODP(A,360) is between 0 and 360, always.
//
//  Examples:
//
//        I         J     MOD  I_MODP   I_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I_MODP, the nonnegative remainder when I is
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cout << "\n";
    cout << "I_MODP - Fatal error!\n";
    cout << "  I_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//********************************************************************

int i_sign ( int i )

//********************************************************************
//
//  Purpose:
//
//    I_SIGN returns the sign of an integer.
//
//  Discussion:
//
//    The sign of 0 and all positive integers is taken to be +1.
//    The sign of all negative integers is -1.
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer whose sign is desired.
//
//    Output, int I_SIGN, the sign of I.
{
  if ( i < 0 )
  {
    return (-1);
  }
  else
  {
    return 1;
  }

}
//******************************************************************************

int i_wrap ( int ival, int ilo, int ihi )

//******************************************************************************
//
//  Purpose:
//
//    I_WRAP forces an integer to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I  I_WRAP
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i_min ( ilo, ihi );
  jhi = i_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i_modp ( ival - jlo, wide );
  }

  return value;
}
//******************************************************************************

void imat_transpose_print ( int m, int n, int a[], const char *title )

//******************************************************************************
//
//  Purpose:
//
//    IMAT_TRANSPOSE_PRINT prints an integer matrix, transposed.
//
//  Modified:
//
//    31 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Input, const char *TITLE, a title to be printed.
//
{
  imat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//******************************************************************************

void imat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, const char *title )

//******************************************************************************
//
//  Purpose:
//
//    IMAT_TRANSPOSE_PRINT_SOME prints some of an integer matrix, transposed.
//
//  Modified:
//
//    09 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, const char *TITLE, a title for the matrix.
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i_min ( i2hi, m );
    i2hi = i_min ( i2hi, ihi );

    cout << "\n";
//
//  For each row I in the current range...
//
//  Write the header.
//
    cout << "  Row:    ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    j2lo = i_max ( jlo, 1 );
    j2hi = i_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
      cout << setw(5) << j << "  ";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }

  }

  cout << "\n";

  return;
# undef INCX
}
//********************************************************************

void ivec_heap_d ( int n, int a[] )

/*********************************************************************
//
//  Purpose:
//
//    IVEC_HEAP_D reorders an array of integers into a descending heap.
//
//  Definition:
//
//    A heap is an array A with the property that, for every index J,
//    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
//    2*J+1 and 2*J+2 are legal).
//
//  Diagram:
//
//                  A(0)
//                /      \
//            A(1)         A(2)
//          /     \        /  \
//      A(3)       A(4)  A(5) A(6)
//      /  \       /   \
//    A(7) A(8)  A(9) A(10)
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the input array.
//
//    Input/output, int A[N].
//    On input, an unsorted array.
//    On output, the array has been reordered into a heap.
*/
{
  int i;
  int ifree;
  int key;
  int m;
//
//  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
//
  for ( i = (n/2)-1; 0 <= i; i-- )
  {
//
//  Copy the value out of the parent node.
//  Position IFREE is now "open".
//
    key = a[i];
    ifree = i;

    for ( ;; )
    {
//
//  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
//  IFREE.  (One or both may not exist because they equal or exceed N.)
//
      m = 2 * ifree + 1;
//
//  Does the first position exist?
//
      if ( n <= m )
      {
        break;
      }
      else
      {
//
//  Does the second position exist?
//
        if ( m + 1 < n )
        {
//
//  If both positions exist, take the larger of the two values,
//  and update M if necessary.
//
          if ( a[m] < a[m+1] )
          {
            m = m + 1;
          }
        }
//
//  If the large descendant is larger than KEY, move it up,
//  and update IFREE, the location of the free position, and
//  consider the descendants of THIS position.
//
        if ( key < a[m] )
        {
          a[ifree] = a[m];
          ifree = m;
        }
        else
        {
          break;
        }

      }

    }
//
//  When you have stopped shifting items up, return the item you
//  pulled out back to the heap.
//
    a[ifree] = key;

  }

  return;
}
//******************************************************************************

int *ivec_indicator ( int n )

//******************************************************************************
//
//  Purpose:
//
//    IVEC_INDICATOR sets an integer vector to the indicator vector.
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, int IVEC_INDICATOR(N), the initialised array.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }

  return a;
}
//********************************************************************

void ivec_sort_heap_a ( int n, int a[] )

//********************************************************************
//
//  Purpose:
//
//    IVEC_SORT_HEAP_A ascending sorts an array of integers using heap sort.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, int A[N].
//    On input, the array to be sorted;
//    On output, the array has been sorted.
//
{
  int n1;
  int temp;

  if ( n <= 1 )
  {
    return;
  }
//
//  1: Put A into descending heap form.
//
  ivec_heap_d ( n, a );
//
//  2: Sort A.
//
//  The largest object in the heap is in A[0].
//  Move it to position A[N-1].
//
  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;
//
//  Consider the diminished heap of size N1.
//
  for ( n1 = n-1; 2 <= n1; n1-- )
  {
//
//  Restore the heap structure of the initial N1 entries of A.
//
    ivec_heap_d ( n1, a );
//
//  Take the largest object from A[0] and move it to A[N1-1].
//
    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;

  }

  return;
}
//******************************************************************************

void ivec_sorted_unique ( int n, int a[], int *nuniq )

//******************************************************************************
//
//  Purpose:
//
//    IVEC_SORTED_UNIQUE finds unique elements in a sorted integer array.
//
//  Modified:
//
//    02 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements in A.
//
//    Input/output, int A[N].  On input, the sorted
//    integer array.  On output, the unique elements in A.
//
//    Output, int *NUNIQ, the number of unique elements in A.
//
{
  int i;

  *nuniq = 0;

  if ( n <= 0 )
  {
    return;
  }

  *nuniq = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] != a[*nuniq] )
    {
      *nuniq = *nuniq + 1;
      a[*nuniq] = a[i];
    }

  }

  return;
}
//******************************************************************************

int lrline ( double xu, double yu, double xv1, double yv1, double xv2,
  double yv2, double dv )

//******************************************************************************
//
//  Purpose:
//
//    LRLINE determines where a point lies in relation to a directed line.
//
//  Discussion:
//
//    LRLINE determines whether a point is to the left of, right of,
//    or on a directed line parallel to a line through given points.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    Barry Joe,
//    Department of Computing Science,
//    University of Alberta,
//    Edmonton, Alberta, Canada  T6G 2H1
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double XU, YU, XV1, YV1, XV2, YV2, are vertex coordinates; the
//    directed line is parallel to and at signed distance DV to the left of
//    the directed line from (XV1,YV1) to (XV2,YV2); (XU,YU) is the vertex for
//    which the position relative to the directed line is to be determined.
//
//    Input, double DV, the signed distance, positive for left.
//
//    Output, int LRLINE, is +1, 0, or -1 depending on whether (XU,YU) is
//    to the right of, on, or left of the directed line.  LRLINE is 0 if
//    the line degenerates to a point.
//
{
  double dx;
  double dxu;
  double dy;
  double dyu;
  double t;
  double tol = 0.0000001;
  double tolabs;
  int value = 0;

  dx = xv2 - xv1;
  dy = yv2 - yv1;
  dxu = xu - xv1;
  dyu = yu - yv1;

  tolabs = tol * d_max ( fabs ( dx ),
                 d_max ( fabs ( dy ),
                 d_max ( fabs ( dxu ),
                 d_max ( fabs ( dyu ), fabs ( dv ) ) ) ) );

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy );

  if ( tolabs < t )
  {
    value = 1;
  }
  else if ( -tolabs <= t )
  {
    value = 0;
  }
  else if ( t < -tolabs )
  {
    value = -1;
  }

  return value;
}
//******************************************************************************

bool perm_check ( int n, int p[] )

//******************************************************************************
//
//  Purpose:
//
//    PERM_CHECK checks that a vector represents a permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from 1
//    to N occurs among the N entries of the permutation.
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the array to check.
//
//    Output, bool PERM_CHECK, is TRUE if the permutation is OK.
//
{
  bool found;
  int i;
  int seek;

  for ( seek = 1; seek <= n; seek++ )
  {
    found = false;

    for ( i = 0; i < n; i++ )
    {
      if ( p[i] == seek )
      {
        found = true;
        break;
      }
    }

    if ( !found )
    {
      return false;
    }

  }

  return true;
}
//******************************************************************************

void perm_inv ( int n, int p[] )

//******************************************************************************
//
//  Purpose:
//
//    PERM_INV inverts a permutation "in place".
//
//  Modified:
//
//    13 January 2004
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input/output, int P[N], the permutation, in standard index form.
//    On output, P describes the inverse permutation
//
{
  int i;
  int i0;
  int i1;
  int i2;
  int is;

  if ( n <= 0 )
  {
    cout << "\n";
    cout << "PERM_INV - Fatal error!\n";
    cout << "  Input value of N = " << n << "\n";
    exit ( 1 );
  }

  if ( !perm_check ( n, p ) )
  {
    cout << "\n";
    cout << "PERM_INV - Fatal error!\n";
    cout << "  The input array does not represent\n";
    cout << "  a proper permutation.\n";
    exit ( 1 );
  }

  is = 1;

  for ( i = 1; i <= n; i++ )
  {
    i1 = p[i-1];

    while ( i < i1 )
    {
      i2 = p[i1-1];
      p[i1-1] = -i2;
      i1 = i2;
    }

    is = - i_sign ( p[i-1] );
    p[i-1] = i_sign ( is ) * abs ( p[i-1] );
  }

  for ( i = 1; i <= n; i++ )
  {
    i1 = -p[i-1];

    if ( 0 <= i1 )
    {
      i0 = i;

      for ( ; ; )
      {
        i2 = p[i1-1];
        p[i1-1] = i0;

        if ( i2 < 0 )
        {
          break;
        }
        i0 = i1;
        i1 = i2;
      }
    }
  }

  return;
}
//********************************************************************

int *points_delaunay_naive_2d ( int n, double p[], int *ntri )

//********************************************************************
//
//  Purpose:
//
//    POINTS_DELAUNAY_NAIVE_2D computes the Delaunay triangulation in 2D.
//
//  Discussion:
//
//    A naive and inefficient (but extremely simple) method is used.
//
//    This routine is only suitable as a demonstration code for small
//    problems.  Its running time is of order N^4.  Much faster algorithms
//    are available.
//
//    Given a set of nodes in the plane, a triangulation is a set of
//    triples of distinct nodes, forming triangles, so that every
//    point with the convex hull of the set of  nodes is either one
//    of the nodes, or lies on an edge of one or more triangles,
//    or lies within exactly one triangle.
//
//  Modified:
//
//    05 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Joseph O'Rourke,
//    Computational Geometry,
//    Cambridge University Press,
//    Second Edition, 1998, page 187.
//
//  Parameters:
//
//    Input, int N, the number of nodes.  N must be at least 3.
//
//    Input, double P[2*N], the coordinates of the nodes.
//
//    Output, int *NTRI, the number of triangles.
//
//    Output, int POINTS_DELAUNAY_NAIVE_2D[3*NTRI], the indices of the
//    nodes making each triangle.
//
{
  int count;
  int flag;
  int i;
  int j;
  int k;
  int m;
  int pass;
  int *tri = NULL;
  double xn;
  double yn;
  double zn;
  double *z;

  count = 0;

  z = new double [ n ];

  for ( i = 0; i < n; i++ )
  {
    z[i] = p[0+i*2] * p[0+i*2] + p[1+i*2] * p[1+i*2];
  }
//
//  First pass counts triangles,
//  Second pass allocates triangles and sets them.
//
  for ( pass = 1; pass <= 2; pass++ )
  {
    if ( pass == 2 )
    {
      tri = new int[3*count];
    }
    count = 0;
//
//  For each triple (I,J,K):
//
    for ( i = 0; i < n - 2; i++ )
    {
      for ( j = i+1; j < n; j++ )
      {
        for ( k = i+1; k < n; k++ )
        {
          if ( j != k )
          {
            xn = ( p[1+j*2] - p[1+i*2] ) * ( z[k] - z[i] )
               - ( p[1+k*2] - p[1+i*2] ) * ( z[j] - z[i] );
            yn = ( p[0+k*2] - p[0+i*2] ) * ( z[j] - z[i] )
               - ( p[0+j*2] - p[0+i*2] ) * ( z[k] - z[i] );
            zn = ( p[0+j*2] - p[0+i*2] ) * ( p[1+k*2] - p[1+i*2] )
               - ( p[0+k*2] - p[0+i*2] ) * ( p[1+j*2] - p[1+i*2] );

            flag = ( zn < 0 );

            if ( flag )
            {
              for ( m = 0; m < n; m++ )
              {
                flag = flag && ( ( p[0+m*2] - p[0+i*2] ) * xn
                               + ( p[1+m*2] - p[1+i*2] ) * yn
                               + ( z[m] - z[i] ) * zn <= 0 );
              }
            }

            if ( flag )
            {
              if ( pass == 2 )
              {
                tri[0+count*3] = i;
                tri[1+count*3] = j;
                tri[2+count*3] = k;
              }
              count = count + 1;
            }

          }
        }
      }
    }
  }

  *ntri = count;
  delete [] z;

  return tri;
}
//******************************************************************************

int s_len_trim ( const char *s )

//******************************************************************************
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Modified:
//
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, const char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char* t;

  n = strlen ( s );
  t = const_cast<char*>(s) + n - 1;

  while ( 0 < n )
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
//******************************************************************************

int swapec ( int i, int *top, int *btri, int *bedg, int point_num,
  double point_xy[], int tri_num, int tri_vert[], int tri_nabe[],
  int stack[] )

//******************************************************************************
//
//  Purpose:
//
//    SWAPEC swaps diagonal edges until all triangles are Delaunay.
//
//  Discussion:
//
//    The routine swaps diagonal edges in a 2D triangulation, based on
//    the empty circumcircle criterion, until all triangles are Delaunay,
//    given that I is the index of the new vertex added to the triangulation.
//
//  Modified:
//
//    03 September 2003
//
//  Author:
//
//    Barry Joe,
//    Department of Computing Science,
//    University of Alberta,
//    Edmonton, Alberta, Canada  T6G 2H1
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int I, the index of the new vertex.
//
//    Input/output, int *TOP, the index of the top of the stack.
//    On output, TOP is zero.
//
//    Input/output, int *BTRI, *BEDG; on input, if positive, are the
//    triangle and edge indices of a boundary edge whose updated indices
//    must be recorded.  On output, these may be updated because of swaps.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double POINT_XY[POINT_NUM*2], the coordinates of the points.
//
//    Input, int TRI_NUM, the number of triangles.
//
//    Input/output, int TRI_VERT[TRI_NUM*3], the triangle incidence list.
//    May be updated on output because of swaps.
//
//    Input/output, int TRI_NABE[TRI_NUM*3], the triangle neighbor list;
//    negative values are used for links of the counter-clockwise linked
//    list of boundary edges;  May be updated on output because of swaps.
//
//      LINK = -(3*I + J-1) where I, J = triangle, edge index.
//
//    Workspace, int STACK[MAXST]; on input, entries 1 through TOP
//    contain the indices of initial triangles (involving vertex I)
//    put in stack; the edges opposite I should be in interior;  entries
//    TOP+1 through MAXST are used as a stack.
//
//    Output, int SWAPEC, is set to 8 for abnormal return.
//
{
  int a;
  int b;
  int c;
  int e;
  int ee;
  int em1;
  int ep1;
  int f;
  int fm1;
  int fp1;
  int l;
  int r;
  int s;
  int swap;
  int t;
  int tt;
  int u;
  double x;
  double y;
//
//  Determine whether triangles in stack are Delaunay, and swap
//  diagonal edge of convex quadrilateral if not.
//
  x = point_xy[2*(i-1)+0];
  y = point_xy[2*(i-1)+1];

  for ( ; ; )
  {
    if ( *top <= 0 )
    {
      break;
    }

    t = stack[(*top)-1];
    *top = *top - 1;

    if ( tri_vert[3*(t-1)+0] == i )
    {
      e = 2;
      b = tri_vert[3*(t-1)+2];
    }
    else if ( tri_vert[3*(t-1)+1] == i )
    {
      e = 3;
      b = tri_vert[3*(t-1)+0];
    }
    else
    {
      e = 1;
      b = tri_vert[3*(t-1)+1];
    }

    a = tri_vert[3*(t-1)+e-1];
    u = tri_nabe[3*(t-1)+e-1];

    if ( tri_nabe[3*(u-1)+0] == t )
    {
      f = 1;
      c = tri_vert[3*(u-1)+2];
    }
    else if ( tri_nabe[3*(u-1)+1] == t )
    {
      f = 2;
      c = tri_vert[3*(u-1)+0];
    }
    else
    {
      f = 3;
      c = tri_vert[3*(u-1)+1];
    }

    swap = diaedg ( x, y,
      point_xy[2*(a-1)+0], point_xy[2*(a-1)+1],
      point_xy[2*(c-1)+0], point_xy[2*(c-1)+1],
      point_xy[2*(b-1)+0], point_xy[2*(b-1)+1] );

    if ( swap == 1 )
    {
      em1 = i_wrap ( e - 1, 1, 3 );
      ep1 = i_wrap ( e + 1, 1, 3 );
      fm1 = i_wrap ( f - 1, 1, 3 );
      fp1 = i_wrap ( f + 1, 1, 3 );

      tri_vert[3*(t-1)+ep1-1] = c;
      tri_vert[3*(u-1)+fp1-1] = i;
      r = tri_nabe[3*(t-1)+ep1-1];
      s = tri_nabe[3*(u-1)+fp1-1];
      tri_nabe[3*(t-1)+ep1-1] = u;
      tri_nabe[3*(u-1)+fp1-1] = t;
      tri_nabe[3*(t-1)+e-1] = s;
      tri_nabe[3*(u-1)+f-1] = r;

      if ( 0 < tri_nabe[3*(u-1)+fm1-1] )
      {
        *top = *top + 1;
        stack[(*top)-1] = u;
      }

      if ( 0 < s )
      {
        if ( tri_nabe[3*(s-1)+0] == u )
        {
          tri_nabe[3*(s-1)+0] = t;
        }
        else if ( tri_nabe[3*(s-1)+1] == u )
        {
          tri_nabe[3*(s-1)+1] = t;
        }
        else
        {
          tri_nabe[3*(s-1)+2] = t;
        }

        *top = *top + 1;

        if ( point_num < *top )
        {
          return 8;
        }

        stack[(*top)-1] = t;
      }
      else
      {
        if ( u == *btri && fp1 == *bedg )
        {
          *btri = t;
          *bedg = e;
        }

        l = - ( 3 * t + e - 1 );
        tt = t;
        ee = em1;

        while ( 0 < tri_nabe[3*(tt-1)+ee-1] )
        {
          tt = tri_nabe[3*(tt-1)+ee-1];

          if ( tri_vert[3*(tt-1)+0] == a )
          {
            ee = 3;
          }
          else if ( tri_vert[3*(tt-1)+1] == a )
          {
            ee = 1;
          }
          else
          {
            ee = 2;
          }

        }

        tri_nabe[3*(tt-1)+ee-1] = l;

      }

      if ( 0 < r )
      {
        if ( tri_nabe[3*(r-1)+0] == t )
        {
          tri_nabe[3*(r-1)+0] = u;
        }
        else if ( tri_nabe[3*(r-1)+1] == t )
        {
          tri_nabe[3*(r-1)+1] = u;
        }
        else
        {
          tri_nabe[3*(r-1)+2] = u;
        }
      }
      else
      {
        if ( t == *btri && ep1 == *bedg )
        {
          *btri = u;
          *bedg = f;
        }

        l = - ( 3 * u + f - 1 );
        tt = u;
        ee = fm1;

        while ( 0 < tri_nabe[3*(tt-1)+ee-1] )
        {
          tt = tri_nabe[3*(tt-1)+ee-1];

          if ( tri_vert[3*(tt-1)+0] == b )
          {
            ee = 3;
          }
          else if ( tri_vert[3*(tt-1)+1] == b )
          {
            ee = 1;
          }
          else
          {
            ee = 2;
          }
        }
        tri_nabe[3*(tt-1)+ee-1] = l;
      }
    }
  }
  return 0;
}
//**********************************************************************

void timestamp ( void )

//**********************************************************************
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Modified:
//
//    21 August 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 29

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  if ( len != 0 )
  {
    cout << time_buffer << "\n";
  }

  return;
# undef TIME_SIZE
}
//**********************************************************************

char *timestring ( void )

//**********************************************************************
//
//  Purpose:
//
//    TIMESTRING returns the current YMDHMS date as a string.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *TIMESTRING, a string containing the current YMDHMS date.
//
{
# define TIME_SIZE 29

  const struct tm *tm;
  time_t now;
  char *s;

  now = time ( NULL );
  tm = localtime ( &now );

  s = new char[TIME_SIZE];

  strftime ( s, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  return s;
# undef TIME_SIZE
}
//******************************************************************************

double *triangle_circumcenter_2d ( double t[] )

//******************************************************************************
//
//  Purpose:
//
//    TRIANGLE_CIRCUMCENTER_2D computes the circumcenter of a triangle in 2D.
//
//  Discussion:
//
//    The circumcenter of a triangle is the center of the circumcircle, the
//    circle that passes through the three vertices of the triangle.
//
//    The circumcircle contains the triangle, but it is not necessarily the
//    smallest triangle to do so.
//
//    If all angles of the triangle are no greater than 90 degrees, then
//    the center of the circumscribed circle will lie inside the triangle.
//    Otherwise, the center will lie outside the circle.
//
//    The circumcenter is the intersection of the perpendicular bisectors
//    of the sides of the triangle.
//
//    In geometry, the circumcenter of a triangle is often symbolized by "O".
//
//  Modified:
//
//    09 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double *X, *Y, the coordinates of the circumcenter of
//    the triangle.
//
{
# define DIM_NUM 2

  double asq;
  double bot;
  double *center;
  double csq;
  double top1;
  double top2;

  center = new double[DIM_NUM];

  asq = ( t[0+1*2] - t[0+0*2] ) * ( t[0+1*2] - t[0+0*2] )
      + ( t[1+1*2] - t[1+0*2] ) * ( t[1+1*2] - t[1+0*2] );

  csq = ( t[0+2*2] - t[0+0*2] ) * ( t[0+2*2] - t[0+0*2] )
      + ( t[1+2*2] - t[1+0*2] ) * ( t[1+2*2] - t[1+0*2] );

  top1 =  ( t[1+1*2] - t[1+0*2] ) * csq - ( t[1+2*2] - t[1+0*2] ) * asq;
  top2 =  ( t[0+1*2] - t[0+0*2] ) * csq - ( t[0+2*2] - t[0+0*2] ) * asq;

  bot  =  ( t[1+1*2] - t[1+0*2] ) * ( t[0+2*2] - t[0+0*2] )
        - ( t[1+2*2] - t[1+0*2] ) * ( t[0+1*2] - t[0+0*2] );

  center[0] = t[0+0*2] + 0.5 * top1 / bot;
  center[1] = t[1+0*2] + 0.5 * top2 / bot;

  return center;

# undef DIM_NUM
}
//******************************************************************************

bool triangulation_plot_eps ( const char *file_out_name, int g_num,
  double g_xy[], int tri_num, int nod_tri[] )

//******************************************************************************
//
//  Purpose:
//
//    TRIANGULATION_PLOT_EPS plots a triangulation of a pointset.
//
//  Discussion:
//
//    The triangulation is most usually a Delaunay triangulation,
//    but this is not necessary.
//
//    The data can be generated by calling DTRIS2, but this is not
//    necessary.
//
//  Modified:
//
//    08 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, const char *FILE_OUT_NAME, the name of the output file.
//
//    Input, int G_NUM, the number of points.
//
//    Input, double G_XY[G_NUM,2], the coordinates of the points.
//
//    Input, int TRI_NUM, the number of triangles.
//
//    Input, int NOD_TRI[3,TRI_NUM], lists, for each triangle,
//    the indices of the points that form the vertices of the triangle.
//
//    Output, bool TRIANGULATION_PLOT_EPS, is TRUE for success.
//
{
  char *date_time;
  int e;
  ofstream file_out;
  int g;
  int j;
  int k;
  int t;
  double x_max;
  double x_min;
  int x_ps;
  int x_ps_max = 576;
  int x_ps_max_clip = 594;
  int x_ps_min = 36;
  int x_ps_min_clip = 18;
  double y_max;
  double y_min;
  int y_ps;
  int y_ps_max = 666;
  int y_ps_max_clip = 684;
  int y_ps_min = 126;
  int y_ps_min_clip = 108;

  date_time = timestring ( );

  x_max = g_xy[0+0*2];
  x_min = g_xy[0+0*2];
  y_max = g_xy[1+0*2];
  y_min = g_xy[1+0*2];

  for ( g = 0; g < g_num; g++ )
  {
    x_max = d_max ( x_max, g_xy[0+g*2] );
    x_min = d_min ( x_min, g_xy[0+g*2] );
    y_max = d_max ( y_max, g_xy[1+g*2] );
    y_min = d_min ( y_min, g_xy[1+g*2] );
  }
//
//  Plot the Delaunay triangulation.
//
//
//  Open the output file.
//
  file_out.open ( file_out_name );

  if ( !file_out )
  {
    cout << "\n";
    cout << "TRIANGULATION_PLOT_EPS - Fatal error!\n";
    cout << "  Cannot open the output file \"" << file_out_name << "\".\n";
    return false;
  }

  file_out << "%!PS-Adobe-3.0 EPSF-3.0\n";
  file_out << "%%Creator: triangulation_plot_eps.cc\n";
  file_out << "%%Title: " << file_out_name << "\n";
  file_out << "%%CreationDate: " << date_time << "\n";
  file_out << "%%Pages: 1\n";
  file_out << "%%Bounding Box: " << x_ps_min << "  " << y_ps_min << "  "
           << x_ps_max << "  " << y_ps_max << "\n";
  file_out << "%%Document-Fonts: Times-Roman\n";
  file_out << "%%LanguageLevel: 1\n";
  file_out << "%%EndComments\n";
  file_out << "%%BeginProlog\n";
  file_out << "/inch {72 mul} def\n";
  file_out << "%%EndProlog\n";
  file_out << "%%Page: 1 1\n";
  file_out << "save\n";
  file_out << "%\n";
  file_out << "%  Set the RGB line color to very light gray.\n";
  file_out << "%\n";
  file_out << "0.900  0.900  0.900 setrgbcolor\n";
  file_out << "%\n";
  file_out << "%  Draw a gray border around the page.\n";
  file_out << "%\n";
  file_out << "newpath\n";
  file_out << "  " << x_ps_min << "  " << y_ps_min << " moveto\n";
  file_out << "  " << x_ps_max << "  " << y_ps_min << " lineto\n";
  file_out << "  " << x_ps_max << "  " << y_ps_max << " lineto\n";
  file_out << "  " << x_ps_min << "  " << y_ps_max << " lineto\n";
  file_out << "  " << x_ps_min << "  " << y_ps_min << " lineto\n";
  file_out << "stroke\n";
  file_out << "%\n";
  file_out << "%  Set the RGB line color to black.\n";
  file_out << "%\n";
  file_out << "0.000  0.000  0.000 setrgbcolor\n";
  file_out << "%\n";
  file_out << "%  Set the font and its size.\n";
  file_out << "%\n";
  file_out << "/Times-Roman findfont\n";
  file_out << "0.50 inch scalefont\n";
  file_out << "setfont\n";
  file_out << "%\n";
  file_out << "%  Print a title.\n";
  file_out << "%\n";
  file_out << "210  702  moveto\n";
  file_out << "(Triangulation)  show\n";
  file_out << "%\n";
  file_out << "%  Define a clipping polygon.\n";
  file_out << "%\n";
  file_out << "newpath\n";
  file_out << "  " << x_ps_min_clip << "  " << y_ps_min_clip << " moveto\n";
  file_out << "  " << x_ps_max_clip << "  " << y_ps_min_clip << " lineto\n";
  file_out << "  " << x_ps_max_clip << "  " << y_ps_max_clip << " lineto\n";
  file_out << "  " << x_ps_min_clip << "  " << y_ps_max_clip << " lineto\n";
  file_out << "  " << x_ps_min_clip << "  " << y_ps_min_clip << " lineto\n";
  file_out << "clip newpath\n";
  file_out << "%\n";
  file_out << "%  Set the RGB line color to green.\n";
  file_out << "%\n";
  file_out << "0.000  0.750  0.150 setrgbcolor\n";
  file_out << "%\n";
  file_out << "%  Draw the nodes.\n";
  file_out << "%\n";

  for ( g = 0; g < g_num; g++ )
  {
    x_ps = int(
      ( ( x_max - g_xy[0+g*2] ) * double( x_ps_min )
      + ( g_xy[0+g*2] - x_min ) * double( x_ps_max ) )
      / ( x_max - x_min ) );

    y_ps = int(
      ( ( y_max - g_xy[1+g*2] ) * double( y_ps_min )
      + ( g_xy[1+g*2] - y_min ) * double( y_ps_max ) )
      / ( y_max - y_min ) );

    file_out << "newpath " << x_ps << "  "
             << y_ps << " 5 0 360 arc closepath fill\n";
  }

  file_out << "%\n";
  file_out << "%  Set the RGB line color to red.\n";
  file_out << "%\n";
  file_out << "0.900  0.200  0.100 setrgbcolor\n";
  file_out << "%\n";
  file_out << "%  Draw the triangles.\n";
  file_out << "%\n";

  for ( t = 1; t <= tri_num; t++ )
  {
    file_out << "newpath\n";

    for ( j = 1; j <= 4; j++ )
    {
      e = i_wrap ( j, 1, 3 );

      k = nod_tri[3*(t-1)+e-1];

      x_ps = int(
        ( ( x_max - g_xy[0+(k-1)*2] ) * double( x_ps_min )
        + ( g_xy[0+(k-1)*2] - x_min ) * double( x_ps_max ) )
        / ( x_max - x_min ) );

      y_ps = int(
        ( ( y_max - g_xy[1+(k-1)*2] ) * double( y_ps_min )
        + ( g_xy[1+(k-1)*2] - y_min ) * double( y_ps_max ) )
        / ( y_max - y_min ) );

      if ( j == 1 )
      {
        file_out << x_ps << "  " << y_ps << " moveto\n";
      }
      else
      {
        file_out << x_ps << "  " << y_ps << " lineto\n";
      }

    }

    file_out << "stroke\n";

  }

  file_out << "restore  showpage\n";
  file_out << "%\n";
  file_out << "%  End of page.\n";
  file_out << "%\n";
  file_out << "%%Trailer\n";
  file_out << "%%EOF\n";

  file_out.close ( );

  return true;
}
//******************************************************************************

void triangulation_print ( int point_num, double xc[], int tri_num,
  int tri_vert[], int tri_nabe[] )

//******************************************************************************
//
//  Purpose:
//
//    TRIANGULATION_PRINT prints information defining a Delaunay triangulation.
//
//  Discussion:
//
//    Triangulations created by RTRIS include extra information encoded
//    in the negative values of TRI_NABE.
//
//    Because some of the nodes counted in POINT_NUM may not actually be
//    used in the triangulation, I needed to compute the true number
//    of vertices.  I added this calculation on 13 October 2001.
//
//    Ernest Fasse pointed out an error in the indexing of VERTEX_LIST,
//    which was corrected on 19 February 2004.
//
//  Modified:
//
//    19 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double XC[2*POINT_NUM], the point coordinates.
//
//    Input, int TRI_NUM, the number of triangles.
//
//    Input, int TRI_VERT[3*TRI_NUM], the nodes that make up the triangles.
//
//    Input, int TRI_NABE[3*TRI_NUM], the triangle neighbors on each side.
//    If there is no triangle neighbor on a particular side, the value of
//    TRI_NABE should be negative.  If the triangulation data was created by
//    DTRIS2, then there is more information encoded in the negative values.
//
{
# define DIM_NUM 2

  int boundary_num;
  int i;
  int j;
  int k;
  int n1;
  int n2;
  int s;
  int s1;
  int s2;
  bool skip;
  int t;
  int *vertex_list;
  int vertex_num;

  cout << "\n";
  cout << "TRIANGULATION_PRINT\n";
  cout << "  Information defining a triangulation.\n";
  cout << "\n";
  cout << "  The number of points is " << point_num << "\n";

  dmat_transpose_print ( DIM_NUM, point_num, xc, "  Point coordinates" );

  cout << "\n";
  cout << "  The number of triangles is " << tri_num << "\n";
  cout << "\n";
  cout << "  Sets of three points are used as vertices of\n";
  cout << "  the triangles.  For each triangle, the points\n";
  cout << "  are listed in counterclockwise order.\n";

  imat_transpose_print ( 3, tri_num, tri_vert, "  Triangle nodes" );

  cout << "\n";
  cout << "  On each side of a given triangle, there is either\n";
  cout << "  another triangle, or a piece of the convex hull.\n";
  cout << "  For each triangle, we list the indices of the three\n";
  cout << "  neighbors, or (if negative) the codes of the\n";
  cout << "  segments of the convex hull.\n";

  imat_transpose_print ( 3, tri_num, tri_nabe, "  Triangle neighbors" );
//
//  Determine VERTEX_NUM, the number of vertices.  This is not
//  the same as the number of points!
//
  vertex_list = new int[3*tri_num];

  k = 0;
  for ( t = 0; t < tri_num; t++ )
  {
    for ( s = 0; s < 3; s++ )
    {
      vertex_list[k] = tri_vert[s+t*3];
      k = k + 1;
    }
  }

  ivec_sort_heap_a ( 3*tri_num, vertex_list );

  ivec_sorted_unique ( 3*tri_num, vertex_list, &vertex_num );

  delete [] vertex_list;
//
//  Determine the number of boundary points.
//
  boundary_num = 2 * vertex_num - tri_num - 2;

  cout << "\n";
  cout << "  The number of boundary points is " << boundary_num << "\n";
  cout << "\n";
  cout << "  The segments that make up the convex hull can be\n";
  cout << "  determined from the negative entries of the triangle\n";
  cout << "  neighbor list.\n";
  cout << "\n";
  cout << "  # Tri Side  N1  N2\n";
  cout << "\n";

  skip = false;

  k = 0;

  for ( i = 0; i < tri_num; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      if ( tri_nabe[j+i*3] < 0 )
      {
        s = -tri_nabe[j+i*3];
        t = s / 3;

        if ( t < 1 || tri_num < t )
        {
          cout << "\n";
          cout << "  Sorry, this data does not use the DTRIS2\n";
          cout << "  convention for convex hull segments.\n";
          skip = true;
          break;
        }

        s1 = ( s % 3 ) + 1;
        s2 = i_wrap ( s1+1, 1, 3 );
        k = k + 1;
        n1 = tri_vert[s1-1+(t-1)*3];
        n2 = tri_vert[s2-1+(t-1)*3];
        cout << setw(4) << k  << "  "
             << setw(4) << t  << "  "
             << setw(4) << s1 << "  "
             << setw(4) << n1 << "  "
             << setw(4) << n2 << "\n";
      }

    }

    if ( skip )
    {
      break;
    }

  }

  return;
# undef DIM_NUM
}
//******************************************************************************

void vbedg ( double x, double y, int point_num, double point_xy[], int tri_num,
  int tri_vert[], int tri_nabe[], int *ltri, int *ledg, int *rtri, int *redg )

//******************************************************************************
//
//  Purpose:
//
//    VBEDG determines which boundary edges are visible to a point.
//
//  Discussion:
//
//    The point (X,Y) is assumed to be outside the convex hull of the
//    region covered by the 2D triangulation.
//
//  Author:
//
//    Barry Joe,
//    Department of Computing Science,
//    University of Alberta,
//    Edmonton, Alberta, Canada  T6G 2H1
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Modified:
//
//    02 September 2003
//
//  Parameters:
//
//    Input, double X, Y, the coordinates of a point outside the convex hull
//    of the current triangulation.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double POINT_XY[POINT_NUM*2], the coordinates of the vertices.
//
//    Input, int TRI_NUM, the number of triangles.
//
//    Input, int TRI_VERT[TRI_NUM*3], the triangle incidence list.
//
//    Input, int TRI_NABE[TRI_NUM*3], the triangle neighbor list; negative
//    values are used for links of a counter clockwise linked list of boundary
//    edges;
//      LINK = -(3*I + J-1) where I, J = triangle, edge index.
//
//    Input/output, int *LTRI, *LEDG.  If LTRI != 0 then these values are
//    assumed to be already computed and are not changed, else they are updated.
//    On output, LTRI is the index of boundary triangle to the left of the
//    leftmost boundary triangle visible from (X,Y), and LEDG is the boundary
//    edge of triangle LTRI to the left of the leftmost boundary edge visible
//    from (X,Y).  1 <= LEDG <= 3.
//
//    Input/output, int *RTRI.  On input, the index of the boundary triangle
//    to begin the search at.  On output, the index of the rightmost boundary
//    triangle visible from (X,Y).
//
//    Input/output, int *REDG, the edge of triangle RTRI that is visible
//    from (X,Y).  1 <= REDG <= 3.
//
{
  int a;
  double ax;
  double ay;
  int b;
  double bx;
  double by;
  bool done;
  int e;
  int l;
  int lr;
  int t;
//
//  Find the rightmost visible boundary edge using links, then possibly
//  leftmost visible boundary edge using triangle neighbor information.
//
  if ( *ltri == 0 )
  {
    done = false;
    *ltri = *rtri;
    *ledg = *redg;
  }
  else
  {
    done = true;
  }

  for ( ; ; )
  {
    l = -tri_nabe[3*((*rtri)-1)+(*redg)-1];
    t = l / 3;
    e = 1 + l % 3;
    a = tri_vert[3*(t-1)+e-1];

    if ( e <= 2 )
    {
      b = tri_vert[3*(t-1)+e];
    }
    else
    {
      b = tri_vert[3*(t-1)+0];
    }

    ax = point_xy[2*(a-1)+0];
    ay = point_xy[2*(a-1)+1];

    bx = point_xy[2*(b-1)+0];
    by = point_xy[2*(b-1)+1];

    lr = lrline ( x, y, ax, ay, bx, by, 0.0 );

    if ( lr <= 0 )
    {
      break;
    }

    *rtri = t;
    *redg = e;

  }

  if ( done )
  {
    return;
  }

  t = *ltri;
  e = *ledg;

  for ( ; ; )
  {
    b = tri_vert[3*(t-1)+e-1];
    e = i_wrap ( e-1, 1, 3 );

    while ( 0 < tri_nabe[3*(t-1)+e-1] )
    {
      t = tri_nabe[3*(t-1)+e-1];

      if ( tri_vert[3*(t-1)+0] == b )
      {
        e = 3;
      }
      else if ( tri_vert[3*(t-1)+1] == b )
      {
        e = 1;
      }
      else
      {
        e = 2;
      }

    }

    a = tri_vert[3*(t-1)+e-1];
    ax = point_xy[2*(a-1)+0];
    ay = point_xy[2*(a-1)+1];

    bx = point_xy[2*(b-1)+0];
    by = point_xy[2*(b-1)+1];

    lr = lrline ( x, y, ax, ay, bx, by, 0.0 );

    if ( lr <= 0 )
    {
      break;
    }

  }

  *ltri = t;
  *ledg = e;

  return;
}
