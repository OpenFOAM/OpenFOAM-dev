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

#include "IC8H18.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IC8H18, 0);
    addToRunTimeSelectionTable(liquidProperties, IC8H18,);
    addToRunTimeSelectionTable(liquidProperties, IC8H18, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IC8H18::IC8H18()
:
    liquidProperties
    (
        114.231,
        543.96,
        2.5676e+6,
        0.468,
        0.266,
        165.78,
        1.4464e-2,
        372.39,
        0.0,
        0.3031,
        1.4051e+4
    ),
    rho_(67.2363666, 0.27373, 543.96, 0.2846),
    pv_(120.81, -7550, -16.111, 0.017099, 1.0),
    hl_(543.96, 375379.713037617, 0.1549, 0.138, 0.0666, 0.0),
    Cp_
    (
        1219.89652546156,
        1.67205049417409,
        0.00414073237562483,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
       -2743275.10767575,
        1219.89652546156,
        0.836025247087043,
        0.00138024412520828,
        0.0,
        0.0
    ),
    Cpg_(997.10236275617, 4627.4653990598, 1594, 2933.52942721328, 677.94),
    B_
    (
        0.00234936225718063,
       -2.83381919093766,
       -413154.047500241,
       -3.49703670632315e+17,
        3.13750207912038e+19
    ),
    mu_(-15.811, 1282.5, 0.67791, -3.8617e-28, 10.0),
    mug_(1.107e-07, 0.746, 72.4, 0.0),
    kappa_(0.1508, -0.0001712, 0.0, 0.0, 0.0, 0.0),
    kappag_(1.758e-05, 1.3114, 392.9, 0.0),
    sigma_(543.96, 0.047434, 1.1975, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 114.231, 28.0) // note: Same as nHeptane
{}


Foam::IC8H18::IC8H18
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


Foam::IC8H18::IC8H18(const dictionary& dict)
:
    IC8H18()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IC8H18::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
