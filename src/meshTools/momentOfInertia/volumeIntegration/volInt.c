

	/*******************************************************
        *                                                      *
	*  volInt.c                                            *
	*                                                      *
	*  This code computes volume integrals needed for      *
	*  determining mass properties of polyhedral bodies.   *
	*                                                      *
	*  For more information, see the accompanying README   *
	*  file, and the paper                                 *
	*                                                      *
	*  Brian Mirtich, "Fast and Accurate Computation of    *
	*  Polyhedral Mass Properties," journal of graphics    *
	*  tools, volume 1, number 1, 1996.                    *
	*                                                      *
	*  This source code is public domain, and may be used  *
	*  in any way, shape or form, free of charge.          *
	*                                                      *
	*  Copyright 1995 by Brian Mirtich                     *
	*                                                      *
	*  mirtich@cs.berkeley.edu                             *
	*  http://www.cs.berkeley.edu/~mirtich                 *
        *                                                      *
	*******************************************************/

/*
	Revision history

	26 Jan 1996	Program creation.

	 3 Aug 1996	Corrected bug arising when polyhedron density
			is not 1.0.  Changes confined to function main().
			Thanks to Zoran Popovic for catching this one.

	27 May 1997     Corrected sign error in translation of inertia
	                product terms to center of mass frame.  Changes 
			confined to function main().  Thanks to 
			Chris Hecker.
*/



#include <stdio.h>
#include <string.h>
#include <math.h>

/*
   ============================================================================
   constants
   ============================================================================
*/

#define MAX_VERTS 100     /* maximum number of polyhedral vertices */
#define MAX_FACES 100     /* maximum number of polyhedral faces */
#define MAX_POLYGON_SZ 10 /* maximum number of verts per polygonal face */

#define X 0
#define Y 1
#define Z 2

/*
   ============================================================================
   macros
   ============================================================================
*/

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

/*
   ============================================================================
   data structures
   ============================================================================
*/

typedef struct {
  int numVerts;
  double norm[3];
  double w;
  int verts[MAX_POLYGON_SZ];
  struct polyhedron *poly;
} FACE;

typedef struct polyhedron {
  int numVerts, numFaces;
  double verts[MAX_VERTS][3];
  FACE faces[MAX_FACES];
} POLYHEDRON;


/*
   ============================================================================
   globals
   ============================================================================
*/

static int A;   /* alpha */
static int B;   /* beta */
static int C;   /* gamma */

/* projection integrals */
static double P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;

/* face integrals */
static double Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;

/* volume integrals */
static double T0, T1[3], T2[3], TP[3];


/*
   ============================================================================
   read in a polyhedron
   ============================================================================
*/

void readPolyhedron(char *name, POLYHEDRON *p)
{
  FILE *fp;
  char line[200], *c;
  int i, j, n;
  double dx1, dy1, dz1, dx2, dy2, dz2, nx, ny, nz, len;
  FACE *f;

  
  if (!(fp = fopen(name, "r"))) {
    printf("i/o error\n");
    exit(1);
  }
  
  fscanf(fp, "%d", &p->numVerts);
  printf("Reading in %d vertices\n", p->numVerts);
  for (i = 0; i < p->numVerts; i++)
    fscanf(fp, "%lf %lf %lf", 
	   &p->verts[i][X], &p->verts[i][Y], &p->verts[i][Z]);

  fscanf(fp, "%d", &p->numFaces);
  printf("Reading in %d faces\n", p->numFaces);
  for (i = 0; i < p->numFaces; i++) {
    f = &p->faces[i];
    f->poly = p;
    fscanf(fp, "%d", &f->numVerts);
    for (j = 0; j < f->numVerts; j++) fscanf(fp, "%d", &f->verts[j]);

    /* compute face normal and offset w from first 3 vertices */
    dx1 = p->verts[f->verts[1]][X] - p->verts[f->verts[0]][X];
    dy1 = p->verts[f->verts[1]][Y] - p->verts[f->verts[0]][Y];
    dz1 = p->verts[f->verts[1]][Z] - p->verts[f->verts[0]][Z];
    dx2 = p->verts[f->verts[2]][X] - p->verts[f->verts[1]][X];
    dy2 = p->verts[f->verts[2]][Y] - p->verts[f->verts[1]][Y];
    dz2 = p->verts[f->verts[2]][Z] - p->verts[f->verts[1]][Z];
    nx = dy1 * dz2 - dy2 * dz1;
    ny = dz1 * dx2 - dz2 * dx1;
    nz = dx1 * dy2 - dx2 * dy1;
    len = sqrt(nx * nx + ny * ny + nz * nz);
    f->norm[X] = nx / len;
    f->norm[Y] = ny / len;
    f->norm[Z] = nz / len;
    f->w = - f->norm[X] * p->verts[f->verts[0]][X]
           - f->norm[Y] * p->verts[f->verts[0]][Y]
           - f->norm[Z] * p->verts[f->verts[0]][Z];

  }

  fclose(fp);

}

/*
   ============================================================================
   compute mass properties
   ============================================================================
*/


/* compute various integrations over projection of face */
void compProjectionIntegrals(FACE *f)
{
  double a0, a1, da;
  double b0, b1, db;
  double a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
  double a1_2, a1_3, b1_2, b1_3;
  double C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
  double Cab, Kab, Caab, Kaab, Cabb, Kabb;
  int i;

  P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0;

  for (i = 0; i < f->numVerts; i++) {
    a0 = f->poly->verts[f->verts[i]][A];
    b0 = f->poly->verts[f->verts[i]][B];
    a1 = f->poly->verts[f->verts[(i+1) % f->numVerts]][A];
    b1 = f->poly->verts[f->verts[(i+1) % f->numVerts]][B];
    da = a1 - a0;
    db = b1 - b0;
    a0_2 = a0 * a0; a0_3 = a0_2 * a0; a0_4 = a0_3 * a0;
    b0_2 = b0 * b0; b0_3 = b0_2 * b0; b0_4 = b0_3 * b0;
    a1_2 = a1 * a1; a1_3 = a1_2 * a1; 
    b1_2 = b1 * b1; b1_3 = b1_2 * b1;

    C1 = a1 + a0;
    Ca = a1*C1 + a0_2; Caa = a1*Ca + a0_3; Caaa = a1*Caa + a0_4;
    Cb = b1*(b1 + b0) + b0_2; Cbb = b1*Cb + b0_3; Cbbb = b1*Cbb + b0_4;
    Cab = 3*a1_2 + 2*a1*a0 + a0_2; Kab = a1_2 + 2*a1*a0 + 3*a0_2;
    Caab = a0*Cab + 4*a1_3; Kaab = a1*Kab + 4*a0_3;
    Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
    Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;

    P1 += db*C1;
    Pa += db*Ca;
    Paa += db*Caa;
    Paaa += db*Caaa;
    Pb += da*Cb;
    Pbb += da*Cbb;
    Pbbb += da*Cbbb;
    Pab += db*(b1*Cab + b0*Kab);
    Paab += db*(b1*Caab + b0*Kaab);
    Pabb += da*(a1*Cabb + a0*Kabb);
  }

  P1 /= 2.0;
  Pa /= 6.0;
  Paa /= 12.0;
  Paaa /= 20.0;
  Pb /= -6.0;
  Pbb /= -12.0;
  Pbbb /= -20.0;
  Pab /= 24.0;
  Paab /= 60.0;
  Pabb /= -60.0;
}

compFaceIntegrals(FACE *f)
{
  double *n, w;
  double k1, k2, k3, k4;

  compProjectionIntegrals(f);

  w = f->w;
  n = f->norm;
  k1 = 1 / n[C]; k2 = k1 * k1; k3 = k2 * k1; k4 = k3 * k1;

  Fa = k1 * Pa;
  Fb = k1 * Pb;
  Fc = -k2 * (n[A]*Pa + n[B]*Pb + w*P1);

  Faa = k1 * Paa;
  Fbb = k1 * Pbb;
  Fcc = k3 * (SQR(n[A])*Paa + 2*n[A]*n[B]*Pab + SQR(n[B])*Pbb
	 + w*(2*(n[A]*Pa + n[B]*Pb) + w*P1));

  Faaa = k1 * Paaa;
  Fbbb = k1 * Pbbb;
  Fccc = -k4 * (CUBE(n[A])*Paaa + 3*SQR(n[A])*n[B]*Paab 
	   + 3*n[A]*SQR(n[B])*Pabb + CUBE(n[B])*Pbbb
	   + 3*w*(SQR(n[A])*Paa + 2*n[A]*n[B]*Pab + SQR(n[B])*Pbb)
	   + w*w*(3*(n[A]*Pa + n[B]*Pb) + w*P1));

  Faab = k1 * Paab;
  Fbbc = -k2 * (n[A]*Pabb + n[B]*Pbbb + w*Pbb);
  Fcca = k3 * (SQR(n[A])*Paaa + 2*n[A]*n[B]*Paab + SQR(n[B])*Pabb
	 + w*(2*(n[A]*Paa + n[B]*Pab) + w*Pa));
}

void compVolumeIntegrals(POLYHEDRON *p)
{
  FACE *f;
  double nx, ny, nz;
  int i;

  T0 = T1[X] = T1[Y] = T1[Z] 
     = T2[X] = T2[Y] = T2[Z] 
     = TP[X] = TP[Y] = TP[Z] = 0;

  for (i = 0; i < p->numFaces; i++) {

    f = &p->faces[i];

    nx = fabs(f->norm[X]);
    ny = fabs(f->norm[Y]);
    nz = fabs(f->norm[Z]);
    if (nx > ny && nx > nz) C = X;
    else C = (ny > nz) ? Y : Z;
    A = (C + 1) % 3;
    B = (A + 1) % 3;

    compFaceIntegrals(f);

    T0 += f->norm[X] * ((A == X) ? Fa : ((B == X) ? Fb : Fc));

    T1[A] += f->norm[A] * Faa;
    T1[B] += f->norm[B] * Fbb;
    T1[C] += f->norm[C] * Fcc;
    T2[A] += f->norm[A] * Faaa;
    T2[B] += f->norm[B] * Fbbb;
    T2[C] += f->norm[C] * Fccc;
    TP[A] += f->norm[A] * Faab;
    TP[B] += f->norm[B] * Fbbc;
    TP[C] += f->norm[C] * Fcca;
  }

  T1[X] /= 2; T1[Y] /= 2; T1[Z] /= 2;
  T2[X] /= 3; T2[Y] /= 3; T2[Z] /= 3;
  TP[X] /= 2; TP[Y] /= 2; TP[Z] /= 2;
}


/*
   ============================================================================
   main
   ============================================================================
*/


int main(int argc, char *argv[])
{
  POLYHEDRON p;
  double density, mass;
  double r[3];            /* center of mass */
  double J[3][3];         /* inertia tensor */

  if (argc != 2) {
    printf("usage:  %s <polyhedron geometry filename>\n", argv[0]);
    exit(0);
  }

  readPolyhedron(argv[1], &p);

  compVolumeIntegrals(&p);


  printf("\nT1 =   %+20.6f\n\n", T0);

  printf("Tx =   %+20.6f\n", T1[X]);
  printf("Ty =   %+20.6f\n", T1[Y]);
  printf("Tz =   %+20.6f\n\n", T1[Z]);
  
  printf("Txx =  %+20.6f\n", T2[X]);
  printf("Tyy =  %+20.6f\n", T2[Y]);
  printf("Tzz =  %+20.6f\n\n", T2[Z]);

  printf("Txy =  %+20.6f\n", TP[X]);
  printf("Tyz =  %+20.6f\n", TP[Y]);
  printf("Tzx =  %+20.6f\n\n", TP[Z]);

  density = 1.0;  /* assume unit density */

  mass = density * T0;

  /* compute center of mass */
  r[X] = T1[X] / T0;
  r[Y] = T1[Y] / T0;
  r[Z] = T1[Z] / T0;

  /* compute inertia tensor */
  J[X][X] = density * (T2[Y] + T2[Z]);
  J[Y][Y] = density * (T2[Z] + T2[X]);
  J[Z][Z] = density * (T2[X] + T2[Y]);
  J[X][Y] = J[Y][X] = - density * TP[X];
  J[Y][Z] = J[Z][Y] = - density * TP[Y];
  J[Z][X] = J[X][Z] = - density * TP[Z];

  /* translate inertia tensor to center of mass */
  J[X][X] -= mass * (r[Y]*r[Y] + r[Z]*r[Z]);
  J[Y][Y] -= mass * (r[Z]*r[Z] + r[X]*r[X]);
  J[Z][Z] -= mass * (r[X]*r[X] + r[Y]*r[Y]);
  J[X][Y] = J[Y][X] += mass * r[X] * r[Y]; 
  J[Y][Z] = J[Z][Y] += mass * r[Y] * r[Z]; 
  J[Z][X] = J[X][Z] += mass * r[Z] * r[X]; 

  printf("center of mass:  (%+12.6f,%+12.6f,%+12.6f)\n\n", r[X], r[Y], r[Z]);

  printf("inertia tensor with origin at c.o.m. :\n");
  printf("%+15.6f  %+15.6f  %+15.6f\n", J[X][X], J[X][Y], J[X][Z]);
  printf("%+15.6f  %+15.6f  %+15.6f\n", J[Y][X], J[Y][Y], J[Y][Z]);
  printf("%+15.6f  %+15.6f  %+15.6f\n\n", J[Z][X], J[Z][Y], J[Z][Z]);
  
}
