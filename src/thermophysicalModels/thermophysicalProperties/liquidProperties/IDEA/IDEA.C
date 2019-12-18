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

#include "IDEA.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IDEA, 0);
    addToRunTimeSelectionTable(liquidProperties, IDEA,);
    addToRunTimeSelectionTable(liquidProperties, IDEA, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IDEA::IDEA()
:
    liquidProperties
    (
        142.26,
        618.074,
        2.11e+6,
        0.523,
        0.247,
        242.67,
        3.4929e-2,
        447.3,
        1.7012e-30,
        0.3478,
        1.57e+4
    ),
    rho_(152.012105, 3.87150382e-1, 618.073893, 4.00790044e-1),
    pv_
    (
        8.4817774623e+01,
       -8.6782398353e+03,
       -9.1277694857,
        4.6153144498e-06,
        2.0
    ),
    hl_
    (
        618.074,
        2.1671983789e+05,
       -4.2413153435e+00,
        1.1656811532e+01,
       -1.1656446689e+01,
        4.3667661492
    ),
    Cp_(1.6604957e+3, -6.250871e-1, 6.1778552e-3, 0.0, 0.0, 0.0),
    h_(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    Cpg_
    (
        1.0457515243e+03,
        3.4410492875e+03,
        1.5976862298e+03,
        2.4697705752e+03,
        7.3699710536e+02
    ),
    B_
    (
        0.00337351091119935,
       -4.13606494008504,
       -534560.916470464,
       -1.13364022911762e+19,
        2.80704220402713e+21
    ),
    mu_(-6.9645853822e+01, 4.4390635942e+03, 8.4680722718e+00, 0.0, 0.0),
    mug_(4.2629382158e-08, 8.8144402122e-01, 9.6918097636e+01, 0.0),
    kappa_(2.03684e-01, -2.3168e-04, 0.0, 0.0, 0.0, 0.0),
    kappag_
    (
       -5.664925956707e+02,
        8.896721676320e-01,
       -2.849783998688e+09,
        6.914935658053e+05
    ),
    sigma_
    (
        618.074,
        8.3846525429e-03,
       -1.0044759047e+01,
        2.7261918781e+01,
       -2.5529134309e+01,
        8.6488806234
    ),
    D_(147.18, 20.1, 142.2, 28.0) // note: Same as nHeptane
{}


Foam::IDEA::IDEA
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


Foam::IDEA::IDEA(const dictionary& dict)
:
    IDEA()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IDEA::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
