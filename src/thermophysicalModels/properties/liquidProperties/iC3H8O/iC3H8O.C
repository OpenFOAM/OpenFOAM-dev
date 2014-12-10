/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "iC3H8O.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(iC3H8O, 0);
    addToRunTimeSelectionTable(liquidProperties, iC3H8O,);
    addToRunTimeSelectionTable(liquidProperties, iC3H8O, Istream);
    addToRunTimeSelectionTable(liquidProperties, iC3H8O, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::iC3H8O::iC3H8O()
:
    liquidProperties
    (
        60.096,
        508.31,
        4.7643e+6,
        0.22013,
        0.248,
        185.28,
        3.20e-2,
        355.41,
        5.5372e-30,
        0.6689,
        2.3575e+4
    ),
    rho_(70.91328, 0.26475, 508.31, 0.243),
    pv_(92.935, -8177.1, -10.031, 3.9988e-06, 2.0),
    hl_(508.31, 948149.627263046, 0.087, 0.3007, 0.0, 0.0),
    Cp_
    (
        7760.91586794462,
       -68.3672790202343,
        0.241380457933972,
       -0.000235057241746539,
        0.0,
        0.0
    ),
    h_
    (
       -6227786.27583977,
        7760.91586794462,
       -34.1836395101172,
        0.0804601526446574,
       -5.87643104366347e-05,
        0.0
    ),
    Cpg_(789.73642172524, 3219.8482428115, 1124, 1560.83599574015, 460.0),
    B_
    (
        0.000502529286474973,
       -0.104665867944622,
       -717185.83599574,
        3.3047124600639e+18,
       -1.43270766773163e+21
    ),
    mu_(-8.23, 2282.2, -0.98495, 0.0, 0.0),
    mug_(1.993e-07, 0.7233, 178.0, 0.0),
    K_(0.2029, -0.0002278, 0.0, 0.0, 0.0, 0.0),
    Kg_(-80.642, -1.4549, -604.42, 0.0),
    sigma_(0.03818, -3.818e-05, -6.51e-08, 0.0, 0.0, 0.0),
    D_(4.75e-10, 1.75, 0.0, 0.0, 0.0)
{}


Foam::iC3H8O::iC3H8O
(
    const liquidProperties& l,
    const NSRDSfunc5& density,
    const NSRDSfunc1& vapourPressure,
    const NSRDSfunc6& heatOfVapourisation,
    const NSRDSfunc0& heatCapacity,
    const NSRDSfunc0& enthalpy,
    const NSRDSfunc7& idealGasHeatCapacity,
    const NSRDSfunc4& secondVirialCoeff,
    const NSRDSfunc1& dynamicViscosity,
    const NSRDSfunc2& vapourDynamicViscosity,
    const NSRDSfunc0& thermalConductivity,
    const NSRDSfunc2& vapourThermalConductivity,
    const NSRDSfunc0& surfaceTension,
    const NSRDSfunc1& vapourDiffussivity
)
:
    liquidProperties(l),
    rho_(density),
    pv_(vapourPressure),
    hl_(heatOfVapourisation),
    Cp_(heatCapacity),
    h_(enthalpy),
    Cpg_(idealGasHeatCapacity),
    B_(secondVirialCoeff),
    mu_(dynamicViscosity),
    mug_(vapourDynamicViscosity),
    K_(thermalConductivity),
    Kg_(vapourThermalConductivity),
    sigma_(surfaceTension),
    D_(vapourDiffussivity)
{}


Foam::iC3H8O::iC3H8O(Istream& is)
:
    liquidProperties(is),
    rho_(is),
    pv_(is),
    hl_(is),
    Cp_(is),
    h_(is),
    Cpg_(is),
    B_(is),
    mu_(is),
    mug_(is),
    K_(is),
    Kg_(is),
    sigma_(is),
    D_(is)
{}


Foam::iC3H8O::iC3H8O(const dictionary& dict)
:
    liquidProperties(dict),
    rho_(dict.subDict("rho")),
    pv_(dict.subDict("pv")),
    hl_(dict.subDict("hl")),
    Cp_(dict.subDict("Cp")),
    h_(dict.subDict("h")),
    Cpg_(dict.subDict("Cpg")),
    B_(dict.subDict("B")),
    mu_(dict.subDict("mu")),
    mug_(dict.subDict("mug")),
    K_(dict.subDict("K")),
    Kg_(dict.subDict("Kg")),
    sigma_(dict.subDict("sigma")),
    D_(dict.subDict("D"))
{}


Foam::iC3H8O::iC3H8O(const iC3H8O& liq)
:
    liquidProperties(liq),
    rho_(liq.rho_),
    pv_(liq.pv_),
    hl_(liq.hl_),
    Cp_(liq.Cp_),
    h_(liq.h_),
    Cpg_(liq.Cpg_),
    B_(liq.B_),
    mu_(liq.mu_),
    mug_(liq.mug_),
    K_(liq.K_),
    Kg_(liq.Kg_),
    sigma_(liq.sigma_),
    D_(liq.D_)
{}


// ************************************************************************* //
