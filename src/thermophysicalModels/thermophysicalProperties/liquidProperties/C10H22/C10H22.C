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

#include "C10H22.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C10H22, 0);
    addToRunTimeSelectionTable(liquidProperties, C10H22,);
    addToRunTimeSelectionTable(liquidProperties, C10H22, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C10H22::C10H22()
:
    liquidProperties
    (
        142.285,
        617.70,
        2.11e+6,
        0.6,
        0.247,
        243.51,
        1.393,
        447.30,
        0.0,
        0.4923,
        1.57e+4
    ),
    rho_(60.94208835, 0.25745, 617.7, 0.28912),
    pv_(112.73, -9749.6, -13.245, 7.1266e-06, 2.0),
    hl_(617.70, 464743.296904101, 0.39797, 0.0, 0.0, 0.0),
    Cp_
    (
        1958.18252099659,
       -1.39094071757388,
        0.00754612221948905,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
       -2699436.15229142,
        1958.18252099659,
       -0.695470358786942,
        0.00251537407316302,
        0.0,
        0.0
    ),
    Cpg_(1175.10630073444, 3762.16748076045, 1614.1, 2658.04547211582, 742),
    B_
    (
        0.00337351091119935,
       -4.13606494008504,
       -534560.916470464,
       -1.13364022911762e+19,
        2.80704220402713e+21
    ),
    mu_(-16.468, 1533.5, 0.7511, 0.0, 0.0),
    mug_(2.64e-08, 0.9487, 71.0, 0.0),
    kappa_(0.2063, -0.000254, 0.0, 0.0, 0.0, 0.0),
    kappag_(-668.4, 0.9323, -4071000000.0, 0.0),
    sigma_(617.70, 0.055435, 1.3095, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 142.285, 28.0) // note: Same as nHeptane
{}


Foam::C10H22::C10H22
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


Foam::C10H22::C10H22(const dictionary& dict)
:
    C10H22()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C10H22::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
