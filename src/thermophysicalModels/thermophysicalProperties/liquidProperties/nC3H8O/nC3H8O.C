/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "nC3H8O.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nC3H8O, 0);
    addToRunTimeSelectionTable(liquidProperties, nC3H8O,);
    addToRunTimeSelectionTable(liquidProperties, nC3H8O, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nC3H8O::nC3H8O()
:
    liquidProperties
    (
        typeName,
        60.096,
        536.71,
        5.1696e+6,
        0.21853,
        0.253,
        146.95,
        6.5112e-7,
        370.35,
        5.6039e-30,
        0.6279,
        2.4557e+4
    ),
    rho_("rho", 75.300288, 0.272, 536.71, 0.2494),
    pv_("pv", 77.46, -7960, -7.5235, 3e-07, 2.0),
    hl_("hl", 536.71, 1098242.8115016, 0.647, -0.783, 0.613, 0.0),
    Cp_
    (
        "Cp",
        216.320553780618,
        18.5203674121406,
       -0.0751797124600639,
        0.000126464323748669,
        0.0,
        0.0
    ),
    h_
    (
        "h",
       -5533091.96851587,
        216.320553780618,
        9.26018370607029,
       -0.0250599041533546,
        3.16160809371672e-05,
        0.0
    ),
    Cpg_
    (
        "Cpg",
        961.794462193823,
        3467.78487752929,
        1542,
        2046.72523961661,
        649
    ),
    B_
    (
        "B",
        0.000933506389776358,
       -1.09325079872204,
       -531649.361022364,
       -2.32627795527157e+17,
       -3.81888977635783e+20
    ),
    mu_("mu", 0.571, 1521, -2.0894, 0.0, 0.0),
    mug_("mug", 7.942e-07, 0.5491, 415.8, 0.0),
    kappa_("kappa", 0.204, -0.000169, 0.0, 0.0, 0.0, 0.0),
    kappag_("kappag", -613.84, 0.7927, -1157400000.0, 0.0),
    sigma_("sigma", 0.04533, -6.88e-05, -1.6e-08, 0.0, 0.0, 0.0),
    D_("D", 4.75e-10, 1.75, 0.0, 0.0, 0.0), // note: same as iC3H8O,
    Hf_(h_.value(Tstd))
{}


Foam::nC3H8O::nC3H8O
(
    const liquidProperties& l,
    const Function1s::NSRDS5& density,
    const Function1s::NSRDS1& vapourPressure,
    const Function1s::NSRDS6& heatOfVapourisation,
    const Function1s::NSRDS0& heatCapacity,
    const Function1s::NSRDS0& enthalpy,
    const Function1s::NSRDS7& idealGasHeatCapacity,
    const Function1s::NSRDS4& secondVirialCoeff,
    const Function1s::NSRDS1& dynamicViscosity,
    const Function1s::NSRDS2& vapourDynamicViscosity,
    const Function1s::NSRDS0& thermalConductivity,
    const Function1s::NSRDS2& vapourThermalConductivity,
    const Function1s::NSRDS0& surfaceTension,
    const Function1s::NSRDS1& vapourDiffusivity
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
    D_(vapourDiffusivity),
    Hf_(h_.value(Tstd))
{}


Foam::nC3H8O::nC3H8O(const dictionary& dict)
:
    nC3H8O()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nC3H8O::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
