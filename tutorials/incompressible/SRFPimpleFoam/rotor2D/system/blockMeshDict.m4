/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  dev                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    `format'      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// General macros to create 2D/extruded-2D meshes

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'print ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])
define(pi, 3.14159265)

define(hex2D, hex ($1b $2b $3b $4b $1t $2t $3t $4t))
define(quad2D, ($1b $2b $2t $1t))
define(frontQuad, ($1t $2t $3t $4t))
define(backQuad, ($1b $4b $3b $2b))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.1;

// Hub radius
define(r, 0.05)

// Impeller-tip radius
define(rb, 0.2)

// Baffle-tip radius
define(Rb, 0.7)

// Tank radius
define(R, 1)

// MRF region radius
define(ri, calc(0.5*(rb + Rb)))

// Thickness of 2D slab
define(z, 0.1)

// Base z
define(Zb, 0)

// Top z
define(Zt, calc(Zb + z))

// Number of cells radially between hub and impeller tip
define(Nr, 6)

// Number of cells radially in each of the two regions between
// impeller and baffle tips
define(Ni, 8)

// Number of cells radially between baffle tip and tank
define(NR, 8)

// Number of cells azimuthally in each of the 8 blocks
define(Na, 12)

// Number of cells in the thickness of the slab
define(Nz, 1)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

define(vert, (x$1$2 y$1$2 $3))
define(evert, (ex$1$2 ey$1$2 $3))

define(a0, 0)
define(a1, -45)
define(a2, -90)
define(a3, -135)
define(a4, 180)
define(a5, 135)
define(a6, 90)
define(a7, 45)

define(ea0, -22.5)
define(ea1, -67.5)
define(ea2, -112.5)
define(ea3, -157.5)
define(ea4, 157.5)
define(ea5, 112.5)
define(ea6, 67.5)
define(ea7, 22.5)

define(ca0, calc(cos((pi/180)*a0)))
define(ca1, calc(cos((pi/180)*a1)))
define(ca2, calc(cos((pi/180)*a2)))
define(ca3, calc(cos((pi/180)*a3)))
define(ca4, calc(cos((pi/180)*a4)))
define(ca5, calc(cos((pi/180)*a5)))
define(ca6, calc(cos((pi/180)*a6)))
define(ca7, calc(cos((pi/180)*a7)))

define(sa0, calc(sin((pi/180)*a0)))
define(sa1, calc(sin((pi/180)*a1)))
define(sa2, calc(sin((pi/180)*a2)))
define(sa3, calc(sin((pi/180)*a3)))
define(sa4, calc(sin((pi/180)*a4)))
define(sa5, calc(sin((pi/180)*a5)))
define(sa6, calc(sin((pi/180)*a6)))
define(sa7, calc(sin((pi/180)*a7)))

define(cea0, calc(cos((pi/180)*ea0)))
define(cea1, calc(cos((pi/180)*ea1)))
define(cea2, calc(cos((pi/180)*ea2)))
define(cea3, calc(cos((pi/180)*ea3)))
define(cea4, calc(cos((pi/180)*ea4)))
define(cea5, calc(cos((pi/180)*ea5)))
define(cea6, calc(cos((pi/180)*ea6)))
define(cea7, calc(cos((pi/180)*ea7)))

define(sea0, calc(sin((pi/180)*ea0)))
define(sea1, calc(sin((pi/180)*ea1)))
define(sea2, calc(sin((pi/180)*ea2)))
define(sea3, calc(sin((pi/180)*ea3)))
define(sea4, calc(sin((pi/180)*ea4)))
define(sea5, calc(sin((pi/180)*ea5)))
define(sea6, calc(sin((pi/180)*ea6)))
define(sea7, calc(sin((pi/180)*ea7)))

define(x00, calc(r*ca0))
define(x01, calc(r*ca1))
define(x02, calc(r*ca2))
define(x03, calc(r*ca3))
define(x04, calc(r*ca4))
define(x05, calc(r*ca5))
define(x06, calc(r*ca6))
define(x07, calc(r*ca7))

define(x10, calc(rb*ca0))
define(x11, calc(rb*ca1))
define(x12, calc(rb*ca2))
define(x13, calc(rb*ca3))
define(x14, calc(rb*ca4))
define(x15, calc(rb*ca5))
define(x16, calc(rb*ca6))
define(x17, calc(rb*ca7))

define(x20, calc(ri*ca0))
define(x21, calc(ri*ca1))
define(x22, calc(ri*ca2))
define(x23, calc(ri*ca3))
define(x24, calc(ri*ca4))
define(x25, calc(ri*ca5))
define(x26, calc(ri*ca6))
define(x27, calc(ri*ca7))

define(x30, calc(Rb*ca0))
define(x31, calc(Rb*ca1))
define(x32, calc(Rb*ca2))
define(x33, calc(Rb*ca3))
define(x34, calc(Rb*ca4))
define(x35, calc(Rb*ca5))
define(x36, calc(Rb*ca6))
define(x37, calc(Rb*ca7))

define(x40, calc(R*ca0))
define(x41, calc(R*ca1))
define(x42, calc(R*ca2))
define(x43, calc(R*ca3))
define(x44, calc(R*ca4))
define(x45, calc(R*ca5))
define(x46, calc(R*ca6))
define(x47, calc(R*ca7))

define(y00, calc(r*sa0))
define(y01, calc(r*sa1))
define(y02, calc(r*sa2))
define(y03, calc(r*sa3))
define(y04, calc(r*sa4))
define(y05, calc(r*sa5))
define(y06, calc(r*sa6))
define(y07, calc(r*sa7))

define(y10, calc(rb*sa0))
define(y11, calc(rb*sa1))
define(y12, calc(rb*sa2))
define(y13, calc(rb*sa3))
define(y14, calc(rb*sa4))
define(y15, calc(rb*sa5))
define(y16, calc(rb*sa6))
define(y17, calc(rb*sa7))

define(y20, calc(ri*sa0))
define(y21, calc(ri*sa1))
define(y22, calc(ri*sa2))
define(y23, calc(ri*sa3))
define(y24, calc(ri*sa4))
define(y25, calc(ri*sa5))
define(y26, calc(ri*sa6))
define(y27, calc(ri*sa7))

define(y30, calc(Rb*sa0))
define(y31, calc(Rb*sa1))
define(y32, calc(Rb*sa2))
define(y33, calc(Rb*sa3))
define(y34, calc(Rb*sa4))
define(y35, calc(Rb*sa5))
define(y36, calc(Rb*sa6))
define(y37, calc(Rb*sa7))

define(y40, calc(R*sa0))
define(y41, calc(R*sa1))
define(y42, calc(R*sa2))
define(y43, calc(R*sa3))
define(y44, calc(R*sa4))
define(y45, calc(R*sa5))
define(y46, calc(R*sa6))
define(y47, calc(R*sa7))

define(ex00, calc(r*cea0))
define(ex01, calc(r*cea1))
define(ex02, calc(r*cea2))
define(ex03, calc(r*cea3))
define(ex04, calc(r*cea4))
define(ex05, calc(r*cea5))
define(ex06, calc(r*cea6))
define(ex07, calc(r*cea7))

define(ex10, calc(rb*cea0))
define(ex11, calc(rb*cea1))
define(ex12, calc(rb*cea2))
define(ex13, calc(rb*cea3))
define(ex14, calc(rb*cea4))
define(ex15, calc(rb*cea5))
define(ex16, calc(rb*cea6))
define(ex17, calc(rb*cea7))

define(ex20, calc(ri*cea0))
define(ex21, calc(ri*cea1))
define(ex22, calc(ri*cea2))
define(ex23, calc(ri*cea3))
define(ex24, calc(ri*cea4))
define(ex25, calc(ri*cea5))
define(ex26, calc(ri*cea6))
define(ex27, calc(ri*cea7))

define(ex30, calc(Rb*cea0))
define(ex31, calc(Rb*cea1))
define(ex32, calc(Rb*cea2))
define(ex33, calc(Rb*cea3))
define(ex34, calc(Rb*cea4))
define(ex35, calc(Rb*cea5))
define(ex36, calc(Rb*cea6))
define(ex37, calc(Rb*cea7))

define(ex40, calc(R*cea0))
define(ex41, calc(R*cea1))
define(ex42, calc(R*cea2))
define(ex43, calc(R*cea3))
define(ex44, calc(R*cea4))
define(ex45, calc(R*cea5))
define(ex46, calc(R*cea6))
define(ex47, calc(R*cea7))

define(ey00, calc(r*sea0))
define(ey01, calc(r*sea1))
define(ey02, calc(r*sea2))
define(ey03, calc(r*sea3))
define(ey04, calc(r*sea4))
define(ey05, calc(r*sea5))
define(ey06, calc(r*sea6))
define(ey07, calc(r*sea7))

define(ey10, calc(rb*sea0))
define(ey11, calc(rb*sea1))
define(ey12, calc(rb*sea2))
define(ey13, calc(rb*sea3))
define(ey14, calc(rb*sea4))
define(ey15, calc(rb*sea5))
define(ey16, calc(rb*sea6))
define(ey17, calc(rb*sea7))

define(ey20, calc(ri*sea0))
define(ey21, calc(ri*sea1))
define(ey22, calc(ri*sea2))
define(ey23, calc(ri*sea3))
define(ey24, calc(ri*sea4))
define(ey25, calc(ri*sea5))
define(ey26, calc(ri*sea6))
define(ey27, calc(ri*sea7))

define(ey30, calc(Rb*sea0))
define(ey31, calc(Rb*sea1))
define(ey32, calc(Rb*sea2))
define(ey33, calc(Rb*sea3))
define(ey34, calc(Rb*sea4))
define(ey35, calc(Rb*sea5))
define(ey36, calc(Rb*sea6))
define(ey37, calc(Rb*sea7))

define(ey40, calc(R*sea0))
define(ey41, calc(R*sea1))
define(ey42, calc(R*sea2))
define(ey43, calc(R*sea3))
define(ey44, calc(R*sea4))
define(ey45, calc(R*sea5))
define(ey46, calc(R*sea6))
define(ey47, calc(R*sea7))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vertices
(
    vert(0, 0, Zb) vlabel(r0b)
    vert(0, 0, Zb) vlabel(r0sb)
    vert(0, 1, Zb) vlabel(r1b)
    vert(0, 2, Zb) vlabel(r2b)
    vert(0, 2, Zb) vlabel(r2sb)
    vert(0, 3, Zb) vlabel(r3b)
    vert(0, 4, Zb) vlabel(r4b)
    vert(0, 4, Zb) vlabel(r4sb)
    vert(0, 5, Zb) vlabel(r5b)
    vert(0, 6, Zb) vlabel(r6b)
    vert(0, 6, Zb) vlabel(r6sb)
    vert(0, 7, Zb) vlabel(r7b)

    vert(1, 0, Zb) vlabel(rb0b)
    vert(1, 1, Zb) vlabel(rb1b)
    vert(1, 2, Zb) vlabel(rb2b)
    vert(1, 3, Zb) vlabel(rb3b)
    vert(1, 4, Zb) vlabel(rb4b)
    vert(1, 5, Zb) vlabel(rb5b)
    vert(1, 6, Zb) vlabel(rb6b)
    vert(1, 7, Zb) vlabel(rb7b)

    vert(2, 0, Zb) vlabel(ri0b)
    vert(2, 1, Zb) vlabel(ri1b)
    vert(2, 2, Zb) vlabel(ri2b)
    vert(2, 3, Zb) vlabel(ri3b)
    vert(2, 4, Zb) vlabel(ri4b)
    vert(2, 5, Zb) vlabel(ri5b)
    vert(2, 6, Zb) vlabel(ri6b)
    vert(2, 7, Zb) vlabel(ri7b)

    vert(3, 0, Zb) vlabel(Rb0b)
    vert(3, 1, Zb) vlabel(Rb1b)
    vert(3, 2, Zb) vlabel(Rb2b)
    vert(3, 3, Zb) vlabel(Rb3b)
    vert(3, 4, Zb) vlabel(Rb4b)
    vert(3, 5, Zb) vlabel(Rb5b)
    vert(3, 6, Zb) vlabel(Rb6b)
    vert(3, 7, Zb) vlabel(Rb7b)

    vert(4, 0, Zb) vlabel(R0b)
    vert(4, 1, Zb) vlabel(R1b)
    vert(4, 1, Zb) vlabel(R1b)
    vert(4, 2, Zb) vlabel(R2b)
    vert(4, 3, Zb) vlabel(R3b)
    vert(4, 3, Zb) vlabel(R3b)
    vert(4, 4, Zb) vlabel(R4b)
    vert(4, 5, Zb) vlabel(R5b)
    vert(4, 5, Zb) vlabel(R5b)
    vert(4, 6, Zb) vlabel(R6b)
    vert(4, 7, Zb) vlabel(R7b)
    vert(4, 7, Zb) vlabel(R7b)

    vert(0, 0, Zt) vlabel(r0t)
    vert(0, 0, Zt) vlabel(r0st)
    vert(0, 1, Zt) vlabel(r1t)
    vert(0, 2, Zt) vlabel(r2t)
    vert(0, 2, Zt) vlabel(r2st)
    vert(0, 3, Zt) vlabel(r3t)
    vert(0, 4, Zt) vlabel(r4t)
    vert(0, 4, Zt) vlabel(r4st)
    vert(0, 5, Zt) vlabel(r5t)
    vert(0, 6, Zt) vlabel(r6t)
    vert(0, 6, Zt) vlabel(r6st)
    vert(0, 7, Zt) vlabel(r7t)

    vert(1, 0, Zt) vlabel(rb0t)
    vert(1, 1, Zt) vlabel(rb1t)
    vert(1, 2, Zt) vlabel(rb2t)
    vert(1, 3, Zt) vlabel(rb3t)
    vert(1, 4, Zt) vlabel(rb4t)
    vert(1, 5, Zt) vlabel(rb5t)
    vert(1, 6, Zt) vlabel(rb6t)
    vert(1, 7, Zt) vlabel(rb7t)

    vert(2, 0, Zt) vlabel(ri0t)
    vert(2, 1, Zt) vlabel(ri1t)
    vert(2, 2, Zt) vlabel(ri2t)
    vert(2, 3, Zt) vlabel(ri3t)
    vert(2, 4, Zt) vlabel(ri4t)
    vert(2, 5, Zt) vlabel(ri5t)
    vert(2, 6, Zt) vlabel(ri6t)
    vert(2, 7, Zt) vlabel(ri7t)

    vert(3, 0, Zt) vlabel(Rb0t)
    vert(3, 1, Zt) vlabel(Rb1t)
    vert(3, 2, Zt) vlabel(Rb2t)
    vert(3, 3, Zt) vlabel(Rb3t)
    vert(3, 4, Zt) vlabel(Rb4t)
    vert(3, 5, Zt) vlabel(Rb5t)
    vert(3, 6, Zt) vlabel(Rb6t)
    vert(3, 7, Zt) vlabel(Rb7t)

    vert(4, 0, Zt) vlabel(R0t)
    vert(4, 1, Zt) vlabel(R1t)
    vert(4, 1, Zt) vlabel(R1t)
    vert(4, 2, Zt) vlabel(R2t)
    vert(4, 3, Zt) vlabel(R3t)
    vert(4, 3, Zt) vlabel(R3t)
    vert(4, 4, Zt) vlabel(R4t)
    vert(4, 5, Zt) vlabel(R5t)
    vert(4, 5, Zt) vlabel(R5t)
    vert(4, 6, Zt) vlabel(R6t)
    vert(4, 7, Zt) vlabel(R7t)
    vert(4, 7, Zt) vlabel(R7t)
);

blocks
(
    // block0
    hex2D(r0, r1, rb1, rb0)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block1
    hex2D(r1, r2s, rb2, rb1)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block2
    hex2D(r2, r3, rb3, rb2)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block3
    hex2D(r3, r4s, rb4, rb3)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block4
    hex2D(r4, r5, rb5, rb4)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block5
    hex2D(r5, r6s, rb6, rb5)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block6
    hex2D(r6, r7, rb7, rb6)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block7
    hex2D(r7, r0s, rb0, rb7)
    rotor
    (Na Nr Nz)
    simpleGrading (1 1 1)

    // block0
    hex2D(rb0, rb1, ri1, ri0)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block1
    hex2D(rb1, rb2, ri2, ri1)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block2
    hex2D(rb2, rb3, ri3, ri2)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block3
    hex2D(rb3, rb4, ri4, ri3)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block4
    hex2D(rb4, rb5, ri5, ri4)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block5
    hex2D(rb5, rb6, ri6, ri5)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block6
    hex2D(rb6, rb7, ri7, ri6)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block7
    hex2D(rb7, rb0, ri0, ri7)
    rotor
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block0
    hex2D(ri0, ri1, Rb1, Rb0)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block1
    hex2D(ri1, ri2, Rb2, Rb1)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block2
    hex2D(ri2, ri3, Rb3, Rb2)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block3
    hex2D(ri3, ri4, Rb4, Rb3)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block4
    hex2D(ri4, ri5, Rb5, Rb4)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block5
    hex2D(ri5, ri6, Rb6, Rb5)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block6
    hex2D(ri6, ri7, Rb7, Rb6)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block7
    hex2D(ri7, ri0, Rb0, Rb7)
    (Na Ni Nz)
    simpleGrading (1 1 1)

    // block0
    hex2D(Rb0, Rb1, R1, R0)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block1
    hex2D(Rb1, Rb2, R2, R1)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block2
    hex2D(Rb2, Rb3, R3, R2)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block3
    hex2D(Rb3, Rb4, R4, R3)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block4
    hex2D(Rb4, Rb5, R5, R4)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block5
    hex2D(Rb5, Rb6, R6, R5)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block6
    hex2D(Rb6, Rb7, R7, R6)
    (Na NR Nz)
    simpleGrading (1 1 1)

    // block7
    hex2D(Rb7, Rb0, R0, R7)
    (Na NR Nz)
    simpleGrading (1 1 1)
);

edges
(
    arc r0b r1b evert(0, 0, Zb)
    arc r1b r2sb evert(0, 1, Zb)
    arc r2b r3b evert(0, 2, Zb)
    arc r3b r4sb evert(0, 3, Zb)
    arc r4b r5b evert(0, 4, Zb)
    arc r5b r6sb evert(0, 5, Zb)
    arc r6b r7b evert(0, 6, Zb)
    arc r7b r0sb evert(0, 7, Zb)

    arc rb0b rb1b evert(1, 0, Zb)
    arc rb1b rb2b evert(1, 1, Zb)
    arc rb2b rb3b evert(1, 2, Zb)
    arc rb3b rb4b evert(1, 3, Zb)
    arc rb4b rb5b evert(1, 4, Zb)
    arc rb5b rb6b evert(1, 5, Zb)
    arc rb6b rb7b evert(1, 6, Zb)
    arc rb7b rb0b evert(1, 7, Zb)

    arc ri0b ri1b evert(2, 0, Zb)
    arc ri1b ri2b evert(2, 1, Zb)
    arc ri2b ri3b evert(2, 2, Zb)
    arc ri3b ri4b evert(2, 3, Zb)
    arc ri4b ri5b evert(2, 4, Zb)
    arc ri5b ri6b evert(2, 5, Zb)
    arc ri6b ri7b evert(2, 6, Zb)
    arc ri7b ri0b evert(2, 7, Zb)

    arc Rb0b Rb1b evert(3, 0, Zb)
    arc Rb1b Rb2b evert(3, 1, Zb)
    arc Rb2b Rb3b evert(3, 2, Zb)
    arc Rb3b Rb4b evert(3, 3, Zb)
    arc Rb4b Rb5b evert(3, 4, Zb)
    arc Rb5b Rb6b evert(3, 5, Zb)
    arc Rb6b Rb7b evert(3, 6, Zb)
    arc Rb7b Rb0b evert(3, 7, Zb)

    arc R0b R1b evert(4, 0, Zb)
    arc R1b R2b evert(4, 1, Zb)
    arc R2b R3b evert(4, 2, Zb)
    arc R3b R4b evert(4, 3, Zb)
    arc R4b R5b evert(4, 4, Zb)
    arc R5b R6b evert(4, 5, Zb)
    arc R6b R7b evert(4, 6, Zb)
    arc R7b R0b evert(4, 7, Zb)

    arc r0t r1t evert(0, 0, Zt)
    arc r1t r2st evert(0, 1, Zt)
    arc r2t r3t evert(0, 2, Zt)
    arc r3t r4st evert(0, 3, Zt)
    arc r4t r5t evert(0, 4, Zt)
    arc r5t r6st evert(0, 5, Zt)
    arc r6t r7t evert(0, 6, Zt)
    arc r7t r0st evert(0, 7, Zt)

    arc rb0t rb1t evert(1, 0, Zt)
    arc rb1t rb2t evert(1, 1, Zt)
    arc rb2t rb3t evert(1, 2, Zt)
    arc rb3t rb4t evert(1, 3, Zt)
    arc rb4t rb5t evert(1, 4, Zt)
    arc rb5t rb6t evert(1, 5, Zt)
    arc rb6t rb7t evert(1, 6, Zt)
    arc rb7t rb0t evert(1, 7, Zt)

    arc ri0t ri1t evert(2, 0, Zt)
    arc ri1t ri2t evert(2, 1, Zt)
    arc ri2t ri3t evert(2, 2, Zt)
    arc ri3t ri4t evert(2, 3, Zt)
    arc ri4t ri5t evert(2, 4, Zt)
    arc ri5t ri6t evert(2, 5, Zt)
    arc ri6t ri7t evert(2, 6, Zt)
    arc ri7t ri0t evert(2, 7, Zt)

    arc Rb0t Rb1t evert(3, 0, Zt)
    arc Rb1t Rb2t evert(3, 1, Zt)
    arc Rb2t Rb3t evert(3, 2, Zt)
    arc Rb3t Rb4t evert(3, 3, Zt)
    arc Rb4t Rb5t evert(3, 4, Zt)
    arc Rb5t Rb6t evert(3, 5, Zt)
    arc Rb6t Rb7t evert(3, 6, Zt)
    arc Rb7t Rb0t evert(3, 7, Zt)

    arc R0t R1t evert(4, 0, Zt)
    arc R1t R2t evert(4, 1, Zt)
    arc R2t R3t evert(4, 2, Zt)
    arc R3t R4t evert(4, 3, Zt)
    arc R4t R5t evert(4, 4, Zt)
    arc R5t R6t evert(4, 5, Zt)
    arc R6t R7t evert(4, 6, Zt)
    arc R7t R0t evert(4, 7, Zt)
);

patches
(
    wall rotor
    (
        quad2D(r0, r1)
        quad2D(r1, r2s)
        quad2D(r2, r3)
        quad2D(r3, r4s)
        quad2D(r4, r5)
        quad2D(r5, r6s)
        quad2D(r6, r7)
        quad2D(r7, r0s)

        quad2D(r0, rb0)
        quad2D(r0s, rb0)

        quad2D(r2, rb2)
        quad2D(r2s, rb2)

        quad2D(r4, rb4)
        quad2D(r4s, rb4)

        quad2D(r6, rb6)
        quad2D(r6s, rb6)
    )

    patch freestream
    (
        quad2D(R0, R1)
        quad2D(R1, R2)
        quad2D(R2, R3)
        quad2D(R3, R4)
        quad2D(R4, R5)
        quad2D(R5, R6)
        quad2D(R6, R7)
        quad2D(R7, R0)

        // quad2D(R1, Rb1)
        // quad2D(R1, Rb1)

        // quad2D(R3, Rb3)
        // quad2D(R3, Rb3)

        // quad2D(R5, Rb5)
        // quad2D(R5, Rb5)

        // quad2D(R7, Rb7)
        // quad2D(R7, Rb7)
    )

    empty front
    (
        frontQuad(r0, r1, rb1, rb0)
        frontQuad(r1, r2s, rb2, rb1)
        frontQuad(r2, r3, rb3, rb2)
        frontQuad(r3, r4s, rb4, rb3)
        frontQuad(r4, r5, rb5, rb4)
        frontQuad(r5, r6s, rb6, rb5)
        frontQuad(r6, r7, rb7, rb6)
        frontQuad(r7, r0s, rb0, rb7)
        frontQuad(rb0, rb1, ri1, ri0)
        frontQuad(rb1, rb2, ri2, ri1)
        frontQuad(rb2, rb3, ri3, ri2)
        frontQuad(rb3, rb4, ri4, ri3)
        frontQuad(rb4, rb5, ri5, ri4)
        frontQuad(rb5, rb6, ri6, ri5)
        frontQuad(rb6, rb7, ri7, ri6)
        frontQuad(rb7, rb0, ri0, ri7)
        frontQuad(ri0, ri1, Rb1, Rb0)
        frontQuad(ri1, ri2, Rb2, Rb1)
        frontQuad(ri2, ri3, Rb3, Rb2)
        frontQuad(ri3, ri4, Rb4, Rb3)
        frontQuad(ri4, ri5, Rb5, Rb4)
        frontQuad(ri5, ri6, Rb6, Rb5)
        frontQuad(ri6, ri7, Rb7, Rb6)
        frontQuad(ri7, ri0, Rb0, Rb7)
        frontQuad(Rb0, Rb1, R1, R0)
        frontQuad(Rb1, Rb2, R2, R1)
        frontQuad(Rb2, Rb3, R3, R2)
        frontQuad(Rb3, Rb4, R4, R3)
        frontQuad(Rb4, Rb5, R5, R4)
        frontQuad(Rb5, Rb6, R6, R5)
        frontQuad(Rb6, Rb7, R7, R6)
        frontQuad(Rb7, Rb0, R0, R7)
    )

    empty back
    (
        backQuad(r0, r1, rb1, rb0)
        backQuad(r1, r2s, rb2, rb1)
        backQuad(r2, r3, rb3, rb2)
        backQuad(r3, r4s, rb4, rb3)
        backQuad(r4, r5, rb5, rb4)
        backQuad(r5, r6s, rb6, rb5)
        backQuad(r6, r7, rb7, rb6)
        backQuad(r7, r0s, rb0, rb7)
        backQuad(rb0, rb1, ri1, ri0)
        backQuad(rb1, rb2, ri2, ri1)
        backQuad(rb2, rb3, ri3, ri2)
        backQuad(rb3, rb4, ri4, ri3)
        backQuad(rb4, rb5, ri5, ri4)
        backQuad(rb5, rb6, ri6, ri5)
        backQuad(rb6, rb7, ri7, ri6)
        backQuad(rb7, rb0, ri0, ri7)
        backQuad(ri0, ri1, Rb1, Rb0)
        backQuad(ri1, ri2, Rb2, Rb1)
        backQuad(ri2, ri3, Rb3, Rb2)
        backQuad(ri3, ri4, Rb4, Rb3)
        backQuad(ri4, ri5, Rb5, Rb4)
        backQuad(ri5, ri6, Rb6, Rb5)
        backQuad(ri6, ri7, Rb7, Rb6)
        backQuad(ri7, ri0, Rb0, Rb7)
        backQuad(Rb0, Rb1, R1, R0)
        backQuad(Rb1, Rb2, R2, R1)
        backQuad(Rb2, Rb3, R3, R2)
        backQuad(Rb3, Rb4, R4, R3)
        backQuad(Rb4, Rb5, R5, R4)
        backQuad(Rb5, Rb6, R6, R5)
        backQuad(Rb6, Rb7, R7, R6)
        backQuad(Rb7, Rb0, R0, R7)
    )
);

// ************************************************************************* //
