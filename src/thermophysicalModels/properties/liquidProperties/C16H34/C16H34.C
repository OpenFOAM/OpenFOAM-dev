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

#include "C16H34.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C16H34, 0);
    addToRunTimeSelectionTable(liquidProperties, C16H34,);
    addToRunTimeSelectionTable(liquidProperties, C16H34, Istream);
    addToRunTimeSelectionTable(liquidProperties, C16H34, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C16H34::C16H34()
:
    liquidProperties
    (
        226.446,
        720.60,
        1.4186e+6,
        0.93,
        0.22,
        291.32,
        8.7467e-2,
        560.01,
        0.0,
        0.7471,
        1.6052e+4
    ),
    rho_(61.94656776, 0.25442, 720.6, 0.3238),
    pv_(233.1, -17346, -32.251, 0.02407, 1.0),
    hl_(720.60, 430654.548987397, 0.4122, 0.0, 0.0, 0.0),
    Cp_
    (
        3769.90540791182,
       -12.5871068599136,
        0.0247211255663602,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
       -2777201.30410301,
        3769.90540791182,
       -6.29355342995681,
        0.00824037518878673,
        0.0,
        0.0
    ),
    Cpg_(1128.74592618108, 3600.8584828171, -1429.7, 2259.69988429913, 679.0),
    B_
    (
        0.0025091191718997,
       -2.46668079807106,
       -1704070.72767901,
       -3.00623548219001e+19,
        7.07320950690231e+21
    ),
    mu_(-18.388, 2056.8, 0.98681, 0.0, 0.0),
    mug_(1.2463e-07, 0.7322, 395.0, 6000.0),
    K_(0.1963, -0.00019, 0.0, 0.0, 0.0, 0.0),
    Kg_(3.075e-06, 1.552, 678.0, 0.0),
    sigma_(720.60, 0.05699, 1.3929, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 226.446, 28.0) // note: Same as nHeptane
{}


Foam::C16H34::C16H34
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


Foam::C16H34::C16H34(Istream& is)
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


Foam::C16H34::C16H34(const dictionary& dict)
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


Foam::C16H34::C16H34(const C16H34& liq)
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
