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

#include "C2H6.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C2H6, 0);
    addToRunTimeSelectionTable(liquidProperties, C2H6,);
    addToRunTimeSelectionTable(liquidProperties, C2H6, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C2H6::C2H6()
:
    liquidProperties
    (
        typeName,
        30.070,
        305.32,
        4.872e+6,
        0.14550,
        0.279,
        90.35,
        1.13,
        184.55,
        0.0,
        0.0995,
        1.24e+4
    ),
    rho_("rho", 57.499854, 0.27937, 305.32, 0.29187),
    pv_("pv", 51.857, -2598.7, -5.1283, 1.4913e-05, 2.0),
    hl_("hl", 305.32, 701396.740937812, 0.60646, -0.55492, 0.32799, 0.0),
    Cp_
    (
        "Cp",
        305.32,
        8.02554965861611,
        2983.63817758563,
        167.548325566287,
       -343.93389207094
    ),
    h_("h", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    Cpg_
    (
        "Cpg",
        1341.07083471899,
        4463.58496840705,
        1655.5,
        2435.08480212837,
        752.87
    ),
    B_
    (
        "B",
        0.00269205187894912,
       -2.05221150648487,
       -47721.9820419022,
        2.24808779514466e+15,
       -3.23910874625873e+17
    ),
    mu_("mu", -3.4134, 197.05, -1.2193, -9.2023e-26, 10.0),
    mug_("mug", 2.5906e-07, 0.67988, 98.902, 0.0),
    kappa_("kappa", 0.35758, -0.0011458, 6.1866e-07, 0.0, 0.0, 0.0),
    kappag_("kappag", 7.3869e-05, 1.1689, 500.73, 0.0),
    sigma_("sigma", 305.32, 0.048643, 1.1981, 0.0, 0.0, 0.0),
    D_("D", 147.18, 20.1, 30.070, 28), // note: Same as nHeptane
    Hf_(h_.value(Tstd))
{}


Foam::C2H6::C2H6
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
    D_(vapourDiffusivity),
    Hf_(h_.value(Tstd))
{}


Foam::C2H6::C2H6(const dictionary& dict)
:
    C2H6()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C2H6::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
