/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  dev
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// General m4 macros

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'print ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])
define(pi, 3.1415926536)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// User-defined parameters

convertToMeters 1;

// Pipe radius (z-axis)
define(radius, 0.015875)

// Pipe length (x-axis)
define(length, 0.44)

// Center (wedge) angle
define(angle, 1)

// Axis origin
define(x0, 0.0)
define(y0, 0.0)
define(z0, 0.0)

// Number of cells
define(nx, 220)
define(ny, 1)
define(nz, 16)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Derived parameters

define(deg,angle/2)
define(zr,calc(z0 + (radius*cos((pi/180)*deg))))
define(yp,calc(y0 + (radius*sin((pi/180)*deg))))
define(ym,calc(y0 - (radius*sin((pi/180)*deg))))
define(xl,calc(x0 + length))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Parametric description

vertices
(
    (x0 y0 z0)
    (xl y0 z0)
    (xl y0 z0)
    (x0 y0 z0)
    (x0 ym zr)
    (xl ym zr)
    (xl yp zr)
    (x0 yp zr)
);

blocks
(
    hex (0 1 1 0 4 5 6 7) (nx ny nz) simpleGrading (1 1 1)
);

edges
();

boundary
(
    inlet
    {
        type patch;
        faces
        (
            (0 4 7 0)
        );
    }
    outlet
    {
        type patch;
        faces
        (
            (1 1 6 5)
        );
    }
    walls
    {
        type wall;
        faces
        (
            (4 5 6 7)
        );
    }
    front
    {
        type wedge;
        faces
        (
            (0 1 5 4)
        );
    }
    back
    {
        type wedge;
        faces
        (
            (0 7 6 1)
        );
    }
    axis
    {
        type empty;
        faces
        (
            (0 0 1 1)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
