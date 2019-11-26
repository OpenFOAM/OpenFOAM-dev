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

#include "Ar.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Ar, 0);
    addToRunTimeSelectionTable(liquidProperties, Ar,);
    addToRunTimeSelectionTable(liquidProperties, Ar, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Ar::Ar()
:
    liquidProperties
    (
        39.948,
        150.86,
        4.8981e+6,
        0.07459,
        0.291,
        83.78,
        6.88e+4,
        87.28,
        0.0,
        0.0,
        1.4138e+4
    ),
    rho_(151.922244, 0.286, 150.86, 0.2984),
    pv_(39.233, -1051.7, -3.5895, 5.0444e-05, 2),
    hl_(150.86, 218509.061780314, 0.352, 0.0, 0.0, 0.0),
    Cp_(4562.43116050866, -70.7770101131471, 0.367477721037349, 0.0, 0.0, 0.0),
    h_
    (
       -1460974.49982473,
        4562.43116050866,
       -35.3885050565735,
        0.122492573679116,
        0.0,
        0.0
    ),
    Cpg_(520.326424351657, 0.0, 0.0, 0.0, 0.0, 0.0),
    B_
    (
        0.000952488234705117,
       -0.379993992189847,
       -2022.62941824372,
        4633523580654.85,
        -302893761890458.0
    ),
    mu_(-8.868, 204.3, -0.3831, -1.3e-22, 10.0),
    mug_(8.386e-07, 0.6175, 75.377, -432.5),
    kappa_(0.1819, -0.0003176, -4.11e-06, 0.0, 0.0, 0.0),
    kappag_(0.0001236, 0.8262, -132.8, 16000),
    sigma_(150.86, 0.03823, 1.2927, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 39.948, 28) // note: Same as nHeptane
{}


Foam::Ar::Ar
(
    const liquidProperties& l,
    const thermophysicalFunctions::NSRDS5& density,
    const thermophysicalFunctions::NSRDS1& vapourPressure,
    const thermophysicalFunctions::NSRDS6& heatOfVapourisation,
    const thermophysicalFunctions::NSRDS0& heatCapacity,
    const thermophysicalFunctions::NSRDS0& enthalpy,
    const thermophysicalFunctions::NSRDS0& idealGasHeatCapacity,
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


Foam::Ar::Ar(const dictionary& dict)
:
    Ar()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Ar::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
