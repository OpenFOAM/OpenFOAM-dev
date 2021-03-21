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

#include "C4H10O.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C4H10O, 0);
    addToRunTimeSelectionTable(liquidProperties, C4H10O,);
    addToRunTimeSelectionTable(liquidProperties, C4H10O, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C4H10O::C4H10O()
:
    liquidProperties
    (
        typeName,
        74.123,
        466.70,
        3.6376e+6,
        0.28,
        0.262,
        156.85,
        4.0709e-1,
        307.58,
        3.836e-30,
        0.2846,
        1.5532e+4
    ),
    rho_("rho", 75.2793188, 0.27608, 466.7, 0.29358),
    pv_("pv", 101.03, -6311.5, -12.27, 1.377e-05, 2),
    hl_("hl", 466.70, 566355.921913576, 0.40717, 0, 0, 0),
    Cp_
    (
        "Cp",
        599.004357621791,
        17.5519069654493,
       -0.0742009902459426,
        0.00011822241409549,
        0.0,
        0.0
    ),
    h_
    (
        "h",
       -4312350.92187216,
        599.004357621791,
        8.77595348272466,
       -0.0247336634153142,
        2.95556035238725e-05,
        0.0
    ),
    Cpg_
    (
        "Cpg",
        1163.06679438231,
        3441.57683849817,
        1541.3,
        1938.66950878944,
        -688.9
    ),
    B_
    (
        "B",
        0.00215992337061371,
       -1.810504162001,
       -276972.0599544,
       -2.12349742994752e+17,
        3.1016013922804e+19
    ),
    mu_("mu", 10.197, -63.8, -3.226, 0.0, 0.0),
    mug_("mug", 1.948e-06, 0.41, 495.8, 0.0),
    kappa_("kappa", 0.249, -0.0004005, 0.0, 0.0, 0.0, 0.0),
    kappag_("kappag", -0.0044894, 0.6155, -3266.3, 0.0),
    sigma_("sigma", 466.70, 0.057356, 1.288, 0.0, 0.0, 0.0),
    D_("D", 147.18, 20.1, 74.123, 28), // note: Same as nHeptane
    Hf_(h_.value(Tstd))
{}


Foam::C4H10O::C4H10O
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
    Hf_(h_.value(Tstd))
{}


Foam::C4H10O::C4H10O(const dictionary& dict)
:
    C4H10O()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C4H10O::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
