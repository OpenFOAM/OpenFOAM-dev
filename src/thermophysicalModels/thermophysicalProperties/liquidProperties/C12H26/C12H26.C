/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    rho_(60.53982858, 0.25511, 658.0, 0.29368),
    pv_(137.47, -11976.0, -16.698, 8.0906e-06, 2.0),
    hl_(658.0, 454020.829174935, 0.40681, 0.0, 0.0, 0.0),
    Cp_(2983.53861146661, -8.0352006011577, 0.018207916025784, 0.0, 0.0, 0.0),
    h_
    (
       -2755166.83820769,
        2983.53861146661,
       -4.01760030057885,
        0.00606930534192801,
        0.0,
        0.0
    ),
    Cpg_(1250.16144371778, 3894.02247296552, 1715.5, 2650.67101879792, 777.5),
    B_
    (
        0.00516619896910848,
       -6.40491258556517,
       -295295.236529723,
       -3.22147729807794e+19,
        8.78195117941974e+21
    ),
    mu_(-20.607, 1943, 1.3205, 0.0, 0.0),
    mug_(6.344e-08, 0.8287, 219.5, 0.0),
    kappa_(0.2047, -0.0002326, 0.0, 0.0, 0.0, 0.0),
    kappag_(5.719e-06, 1.4699, 579.4, 0.0),
    sigma_(658.0, 0.055493, 1.3262, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 170.338, 28.0) // note: Same as nHeptane
{}


Foam::C12H26::C12H26
(
    const liquidProperties& l,
    const thermophysicalFunctions::NSRDS5& density,
    const thermophysicalFunctions::NSRDS1& vapourPressure,
    const thermophysicalFunctions::NSRDS6& heatOfVapourisation,
    const thermophysicalFunctions::NSRDS0& heatCapacity,
    const thermophysicalFunctions::NSRDS0& enthalpy,
    const thermophysicalFunctions::NSRDS7& idealGasHeatCapacity,
    const thermophysicalFunctions::NSRDS4& secondVirialCoeff,
    const thermophysicalFunctions::NSRDS1& dynamicViscosity,
    const thermophysicalFunctions::NSRDS2& vapourDynamicViscosity,
    const thermophysicalFunctions::NSRDS0& thermalConductivity,
    const thermophysicalFunctions::NSRDS2& vapourThermalConductivity,
    const thermophysicalFunctions::NSRDS6& surfaceTension,
    const thermophysicalFunctions::APIdiffCoef& vapourDiffusivity
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
