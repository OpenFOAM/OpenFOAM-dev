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

#include "C12H26.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C12H26, 0);
    addToRunTimeSelectionTable(liquidProperties, C12H26,);
    addToRunTimeSelectionTable(liquidProperties, C12H26, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C12H26::C12H26()
:
    liquidProperties
    (
        typeName,
        170.338,
        658.0,
        1.82e+6,
        0.716,
        0.238,
        263.57,
        6.152e-1,
        489.47,
        0.0,
        0.5764,
        1.59e+4
    ),
    rho_("rho", 60.53982858, 0.25511, 658.0, 0.29368),
    pv_("pv", 137.47, -11976.0, -16.698, 8.0906e-06, 2.0),
    hl_("hl", 658.0, 454020.829174935, 0.40681, 0.0, 0.0, 0.0),
    Cp_
    (
        "Cp",
        2983.53861146661,
        -8.0352006011577,
        0.018207916025784,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
        "h",
       -2755166.83820769,
        2983.53861146661,
       -4.01760030057885,
        0.00606930534192801,
        0.0,
        0.0
    ),
    Cpg_
    (
        "Cpg",
        1250.16144371778,
        3894.02247296552,
        1715.5,
        2650.67101879792,
        777.5
    ),
    B_
    (
        "B",
        0.00516619896910848,
       -6.40491258556517,
       -295295.236529723,
       -3.22147729807794e+19,
        8.78195117941974e+21
    ),
    mu_("mu", -20.607, 1943, 1.3205, 0.0, 0.0),
    mug_("mug", 6.344e-08, 0.8287, 219.5, 0.0),
    kappa_("kappa", 0.2047, -0.0002326, 0.0, 0.0, 0.0, 0.0),
    kappag_("kappag", 5.719e-06, 1.4699, 579.4, 0.0),
    sigma_("sigma", 658.0, 0.055493, 1.3262, 0.0, 0.0, 0.0),
    D_("D", 147.18, 20.1, 170.338, 28.0), // note: Same as nHeptane
    hf_(h_.value(Tstd))
{}


Foam::C12H26::C12H26
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


Foam::C12H26::C12H26(const dictionary& dict)
:
    C12H26()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C12H26::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
