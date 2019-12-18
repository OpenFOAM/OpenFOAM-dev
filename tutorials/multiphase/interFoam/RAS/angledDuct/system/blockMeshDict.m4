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
    `format'      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// block definition for a porosity with an angled inlet/outlet
// the porosity is not aligned with the main axes
//
dnl> -----------------------------------------------------------------
dnl> <STANDARD DEFINITIONS>
dnl>
changecom(//)changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'print ($1)')]) dnl>
define(VCOUNT, 0)  dnl>
define(vlabel, [[// ]pt VCOUNT ($1) define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])  dnl>
dnl>
define(hex2D, hex ($1b $2b $3b $4b $1f $2f $3f $4f)) dnl>
define(quad2D, ($1f $1b $2b $2f))  dnl>
define(frontQuad, ($1f $2f $3f $4f)) dnl>
define(backQuad, ($4b $3b $2b $1b)) dnl>
dnl>
dnl> </STANDARD DEFINITIONS>
dnl> -----------------------------------------------------------------
dnl>
define(ncells, 20) dnl>
define(ninlet, 15) dnl>
define(nporo, 20) dnl>
define(noutlet, 20) dnl>
dnl>
define(x0,0) dnl>
define(y0,0) dnl>
define(y0,0) dnl>
define(Cos,0.7071067812) dnl>   == cos(45)
define(Sin,0.7071067812) dnl>   == sin(45)
dnl>
define(width,50) dnl>
define(zBack,calc(-width/2)) dnl>
define(zFront,calc(width/2)) dnl>
define(leninlet,150)dnl>
define(lenporo,100)dnl>
define(lenoutlet,100)dnl>
dnl>
define(xhyp,calc(Sin*width)) dnl>
define(yhyp,calc(Cos*width)) dnl>
define(xinlet,leninlet)dnl>
define(xporo,calc(Cos*lenporo)) dnl>
define(yporo,calc(Sin*lenporo)) dnl>
define(xoutlet,calc(xporo + Cos*lenoutlet)) dnl>
define(youtlet,calc(yporo + Sin*lenoutlet)) dnl>
dnl>

convertToMeters 0.001;

vertices
(
    // inlet region
    ( -xinlet  y0  zBack )  vlabel(in1b)
    ( -xinlet yhyp  zBack ) vlabel(in2b)
    ( -xinlet  y0  zFront )  vlabel(in1f)
    ( -xinlet yhyp  zFront ) vlabel(in2f)

    // join inlet->outlet
    (  x0 y0  zBack )    vlabel(join1b)
    ( -xhyp   yhyp  zBack ) vlabel(join2b)
    (  x0 y0  zFront )    vlabel(join1f)
    ( -xhyp   yhyp  zFront ) vlabel(join2f)

    // porosity ends ->outlet
    ( xporo yporo  zBack )  vlabel(poro1b)
    ( calc(xporo - xhyp) calc(yporo + yhyp)  zBack )  vlabel(poro2b)
    ( xporo yporo  zFront )  vlabel(poro1f)
    ( calc(xporo - xhyp) calc(yporo + yhyp)  zFront )  vlabel(poro2f)

    // outlet
    ( xoutlet youtlet zBack ) vlabel(out1b)
    ( calc(xoutlet - xhyp) calc(youtlet + yhyp) zBack ) vlabel(out2b)
    ( xoutlet youtlet zFront ) vlabel(out1f)
    ( calc(xoutlet - xhyp) calc(youtlet + yhyp) zFront ) vlabel(out2f)
);

blocks
(
    // inlet block
    hex2D(in1, join1, join2, in2)
    inlet ( ninlet ncells ncells ) simpleGrading (1 1 1)

    // porosity block
    hex2D(join1, poro1, poro2, join2)
    porosity ( nporo ncells ncells ) simpleGrading (1 1 1)

    // outlet block
    hex2D(poro1, out1, out2, poro2)
    outlet ( noutlet ncells ncells )  simpleGrading (1 1 1)
);

defaultPatch
{
    name walls;
    type wall;
}

boundary
(
    porosityWall
    {
        type wall;
        faces
        (
            // porosity block
            frontQuad(join1, poro1, poro2, join2)
            // porosity block
            backQuad(join1, poro1, poro2, join2)
            // porosity block
            quad2D(join1, poro1)
            quad2D(poro2, join2)
        );
    }

    inlet
    {
        type patch;
        faces
        (
            quad2D(in2, in1)
        );
    }

    outlet
    {
        type patch;
        faces
        (
            quad2D(out2, out1)
        );
    }
);

// ************************************************************************* //
