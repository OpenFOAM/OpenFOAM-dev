/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "MB.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MB, 0);
    addToRunTimeSelectionTable(liquidProperties, MB,);
    addToRunTimeSelectionTable(liquidProperties, MB, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MB::MB()
:
    liquidProperties
    (
        typeName,
        102.133,
        554.5,
        3.4734e+6,
        0.34,
        0.256,
        187.35,
        1.0102e-1,
        375.90,
        5.7373e-30,
        0.3807,
        1.7713e+4
    ),
    rho_("rho", 76.6099633, 0.257, 554.5, 0.2772),
    pv_("pv", 107.51, -8112.9, -12.77, 9.2919e-06, 2.0),
    hl_("hl", 554.5, 508307.794738233, 0.392, 0.0, 0.0, 0.0),
    Cp_("Cp", 1135.77394182096, 2.89818178257762, 0.0, 0.0, 0.0, 0.0),
    h_
    (
        "h",
        -5255966.14542938,
        1135.77394182096,
        1.44909089128881,
        0.0,
        0.0,
        0.0
    ),
    Cpg_
    (
        "Cpg",
        875.329227575808,
        2849.22600922327,
        1570.0,
        2029.70636327142,
        678.3
    ),
    B_
    (
        "B",
        0.00220496803188,
       -2.42184210783978,
       -401045.695318849,
       -2.85079259397061e+17,
       -3.57377145486767e+19
    ),
    mu_("mu", -12.206, 1141.7, 0.15014, 0.0, 0.0),
    mug_("mug", 3.733e-07, 0.6177, 256.5, 0.0),
    kappa_("kappa", 0.2298, -0.0003002, 0.0, 0.0, 0.0, 0.0),
    kappag_("kappag", 1333.1, 0.9962, 12317000000.0, 0.0),
    sigma_("sigma", 554.5, 0.064084, 1.2418, 0.0, 0.0, 0.0),
    D_("D", 147.18, 20.1, 102.133, 28.0), // note: Same as nHeptane,
    hf_(h_.value(Tstd))
{}


Foam::MB::MB
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
    hf_(h_.value(Tstd))
{}


Foam::MB::MB(const dictionary& dict)
:
    MB()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MB::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
