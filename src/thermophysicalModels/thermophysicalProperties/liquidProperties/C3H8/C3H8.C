/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "C3H8.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C3H8, 0);
    addToRunTimeSelectionTable(liquidProperties, C3H8,);
    addToRunTimeSelectionTable(liquidProperties, C3H8, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C3H8::C3H8()
:
    liquidProperties
    (
        44.096,
        369.83,
        4.248e+6,
        0.2, 0.276,
        85.47,
        1.685e-4,
        231.11,
        0.0,
        0.1523,
        1.31e+4
    ),
    rho_("rho", 60.6628672, 0.27453, 369.83, 0.29359),
    pv_("pv", 59.078, -3492.6, -6.0669, 1.0919e-05, 2.0),
    hl_("hl", 369.83, 662395.682148041, 0.78237, -0.77319, 0.39246, 0.0),
    Cp_
    (
        "Cp",
        369.83,
        9.48470319647089,
        2576.87772133527,
        95.3560311677331,
       -131.535634282099
    ),
    h_("h", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    Cpg_
    (
        "Cpg",
        1177.43105950653,
        4364.34143686502,
        1626.5,
        2648.76632801161,
        723.6
    ),
    B_
    (
        "B",
        0.00255578737300435,
       -2.24963715529753,
       -102276.850507983,
        7.00743831640058e+15,
       -1.59878447024673e+18
    ),
    mu_("mu", -6.9281, 420.76, -0.63276, -1.713e-26, 10.0),
    mug_("mug", 2.4993e-07, 0.68612, 179.34, -8254.6),
    kappa_("kappa", 0.26755, -0.00066457, 2.774e-07, 0.0, 0.0, 0.0),
    kappag_("kappag", -1.12, 0.10972, -9834.6, -7535800),
    sigma_("sigma", 369.83, 0.05092, 1.2197, 0.0, 0.0, 0.0),
    D_("D", 147.18, 20.1, 44.096, 28) // note: Same as nHeptane
{}


Foam::C3H8::C3H8
(
    const liquidProperties& l,
    const Function1s::NSRDS5& density,
    const Function1s::NSRDS1& vapourPressure,
    const Function1s::NSRDS6& heatOfVapourisation,
    const Function1s::NSRDS14& heatCapacity,
    const Function1s::NSRDS0& enthalpy,
    const Function1s::NSRDS7& idealGasHeatCapacity,
    const Function1s::NSRDS4& secondVirialCoeff,
    const Function1s::NSRDS1& dynamicViscosity,
    const Function1s::NSRDS2& vapourDynamicViscosity,
    const Function1s::NSRDS0& thermalConductivity,
    const Function1s::NSRDS2& vapourThermalConductivity,
    const Function1s::NSRDS6& surfaceTension,
    const Function2s::APIdiffCoef& vapourDiffusivity
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
    kappa_(thermalConductivity),
    kappag_(vapourThermalConductivity),
    sigma_(surfaceTension),
    D_(vapourDiffusivity)
{}


Foam::C3H8::C3H8(const dictionary& dict)
:
    C3H8()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C3H8::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
