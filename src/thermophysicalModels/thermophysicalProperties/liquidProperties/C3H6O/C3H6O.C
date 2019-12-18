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

#include "C3H6O.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C3H6O, 0);
    addToRunTimeSelectionTable(liquidProperties, C3H6O,);
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
       -4.905296049462618e06,
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
    kappa_(0.2502, -0.000298, 0.0, 0.0, 0.0, 0.0),
    kappag_(-26.8, 0.9098, -126500000, 0.0),
    sigma_(508.20, 0.0622, 1.124, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 58.08, 28) // note: Same as nHeptane
{}


Foam::C3H6O::C3H6O
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


Foam::C3H6O::C3H6O(const dictionary& dict)
:
    C3H6O()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C3H6O::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
