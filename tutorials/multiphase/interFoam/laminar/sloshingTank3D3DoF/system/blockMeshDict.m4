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
// General m4 macros

changecom(//)changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')]) dnl>
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

define(hex2D, hex (b$1 b$2 b$3 b$4 f$1 f$2 f$3 f$4))
define(quad2D, (b$1 b$2 f$2 f$1))
define(frontQuad, (f$1 f$2 f$3 f$4))
define(backQuad, (b$1 b$4 b$3 b$2))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// User-defined parameters

convertToMeters 1;

define(l, 20)       // Length of tank (x-direction)
define(b, 40)       // Breadth of tank (y-direction)
define(h, 30)       // Depth of tank (z-direction)

define(hlc, 5)      // Depth to the top (height) of lower chamfer
define(huc, 10)     // Height of upper chamfer

define(thetalc, 45) // Angle of lower chamfer to the horizontal
define(thetauc, 45) // Angle of upper chamfer to the horizontal

define(CofGy, calc(b/2.0))  // Centre of gravity in y-direction
define(CofGz, 10.0)         // Centre of gravity in z-direction

define(Nl, 19)      // Number of cells in the length (1 for 2D)
define(Nb, 40)      // Number of cells in the breadth
define(Nhlc, 6)     // Number of cells in the height of the lower champfer
define(Nh, 16)      // Number of cells in the height between the chamfers
define(Nhuc, 12)    // Number of cells in the height of the upper champfer

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Derived parameters

define(blc, calc(hlc/tan(deg2rad(thetalc)))) // Breadth to the top (height) of lower chamfer
define(buc, calc(huc/tan(deg2rad(thetauc)))) // Breadth of upper chamfer

define(Yl, -CofGy)
define(Yllc, calc(Yl + blc))
define(Yluc, calc(Yl + buc))

define(Yr, calc(Yl + b))
define(Yrlc, calc(Yr - blc))
define(Yruc, calc(Yr - buc))

define(Zb, -CofGz)
define(Zlc, calc(Zb + hlc))
define(Zt, calc(Zb + h))
define(Zuc, calc(Zt - huc))

define(Xf, calc(l/2.0))
define(Xb, calc(Xf - l))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Parametric description

vertices
(
    (Xb Yllc Zb) vlabel(bllcb)
    (Xb Yl Zlc)  vlabel(bllc)
    (Xb Yl Zuc)  vlabel(bluc)
    (Xb Yluc Zt) vlabel(bluct)
    (Xb Yrlc Zb) vlabel(brlcb)
    (Xb Yr Zlc)  vlabel(brlc)
    (Xb Yr Zuc)  vlabel(bruc)
    (Xb Yruc Zt) vlabel(bruct)

    (Xf Yllc Zb) vlabel(fllcb)
    (Xf Yl Zlc)  vlabel(fllc)
    (Xf Yl Zuc)  vlabel(fluc)
    (Xf Yluc Zt) vlabel(fluct)
    (Xf Yrlc Zb) vlabel(frlcb)
    (Xf Yr Zlc)  vlabel(frlc)
    (Xf Yr Zuc)  vlabel(fruc)
    (Xf Yruc Zt) vlabel(fruct)
);

blocks
(
    // block0
    hex2D(llcb, rlcb, rlc, llc)
    (Nb Nhlc Nl)
    simpleGrading (1 1 1)

    // block1
    hex2D(llc, rlc, ruc, luc)
    (Nb Nh Nl)
    simpleGrading (1 1 1)

    // block2
    hex2D(luc, ruc, ruct, luct)
    (Nb Nhuc Nl)
    simpleGrading (1 1 1)
);

patches
(
    patch walls
    (
        quad2D(llcb, rlcb)
        quad2D(rlcb, rlc)
        quad2D(rlc, ruc)
        quad2D(ruc, ruct)
        quad2D(ruct, luct)
        quad2D(luct, luc)
        quad2D(luc, llc)
        quad2D(llc, llcb)
        frontQuad(llcb, rlcb, rlc, llc)
        frontQuad(llc, rlc, ruc, luc)
        frontQuad(luc, ruc, ruct, luct)
        backQuad(llcb, rlcb, rlc, llc)
        backQuad(llc, rlc, ruc, luc)
        backQuad(luc, ruc, ruct, luct)
    )
);

// ************************************************************************* //
