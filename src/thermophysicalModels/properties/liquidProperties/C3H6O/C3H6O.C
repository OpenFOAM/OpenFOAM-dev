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

#include "C3H6O.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C3H6O, 0);
    addToRunTimeSelectionTable(liquidProperties, C3H6O,);
    addToRunTimeSelectionTable(liquidProperties, C3H6O, Istream);
    addToRunTimeSelectionTable(liquidProperties, C3H6O, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C3H6O::C3H6O()
:
    liquidProperties
    (
        58.08,
        508.20,
        4.7015e+6,
        0.209,
        0.233,
        178.45,
        2.5938,
        329.44,
        9.6066e-30,
        0.3064,
        1.9774e+4
    ),
    rho_(71.426784, 0.2576, 508.2, 0.29903),
    pv_(70.72, -5685, -7.351, 6.3e-06, 2.0),
    hl_(508.20, 846590.909090909, 1.036, -1.294, 0.672, 0.0),
    Cp_
    (
        2334.71074380165,
       -3.04752066115702,
        0.00488464187327824,
        1.18629476584022e-05,
        0.0,
        0.0
    ),
    h_
    (
        2571201.780143,
        2334.71074380165,
       -1.52376033057851,
        0.00162821395775941,
        2.96573691460055e-06,
        0.0
    ),
    Cpg_(828.512396694215, 2830.57851239669, 1250.0, 1234.50413223141, -524.4),
    B_
    (
        0.00190599173553719,
       -1.70798898071625,
       -525826.446280992,
        1.70282369146006e+17,
       -2.83298898071625e+20
    ),
    mu_(-14.918, 1023.4, 0.5961, 0.0, 0.0),
    mug_(3.1005e-08, 0.9762, 23.139, 0.0),
    K_(0.2502, -0.000298, 0.0, 0.0, 0.0, 0.0),
    Kg_(-26.8, 0.9098, -126500000, 0.0),
    sigma_(508.20, 0.0622, 1.124, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 58.08, 28) // note: Same as nHeptane
{}


Foam::C3H6O::C3H6O
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
    const NSRDSfunc6& surfaceTension,
    const APIdiffCoefFunc& vapourDiffussivity
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


Foam::C3H6O::C3H6O(Istream& is)
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


Foam::C3H6O::C3H6O(const dictionary& dict)
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


Foam::C3H6O::C3H6O(const C3H6O& liq)
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
