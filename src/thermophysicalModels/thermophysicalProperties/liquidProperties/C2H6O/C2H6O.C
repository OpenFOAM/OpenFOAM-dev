/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "C2H6O.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C2H6O, 0);
    addToRunTimeSelectionTable(liquidProperties, C2H6O,);
    addToRunTimeSelectionTable(liquidProperties, C2H6O, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C2H6O::C2H6O()
:
    liquidProperties
    (
        46.069,
        400.10,
        5.3702e+6,
        0.17,
        0.274,
        131.65,
        3.0849,
        248.31,
        4.3363e-30,
        0.2036,
        1.7572e+4
    ),
    rho_(69.472052, 0.26325, 400.1, 0.2806),
    pv_(51.566, -3664.4, -4.653, 5.9e-06, 2),
    hl_(400.10, 608435.173326966, 0.2477, -0.089, 0.203, 0),
    Cp_
    (
        1491.24139877141,
        11.3099915344375,
       -0.067273003538171,
        0.000136556035511949,
        0.0,
        0.0
    ),
    h_
    (
       -5024829.22619402,
        1491.24139877141,
        5.65499576721874,
       -0.0224243345127237,
        3.41390088779874e-05,
        0.0
    ),
    Cpg_(950.747791356443, 3160.47667628991, 1284, 1291.5409494454, 520),
    B_
    (
        0.00235082159369641,
       -2.26616596843865,
       -123293.320888233,
       -8.87364605266014e+16,
        1.46389111984198e+19
    ),
    mu_(-10.62, 448.99, 8.3967e-05, 0.0, 0.0),
    mug_(7.27, 0.1091, 440600000, 0.0),
    kappa_(0.31276, -0.0005677, 0.0, 0.0, 0.0, 0.0),
    kappag_(0.2247, 0.1026, 997.06, 1762900),
    sigma_(400.10, 0.06096, 1.2286, 0, 0, 0),
    D_(147.18, 20.1, 46.069, 28) // note: Same as nHeptane
{}


Foam::C2H6O::C2H6O
(
    const liquidProperties& l,
    const NSRDSfunc5& density,
    const NSRDSfunc1& vapourPressure,
    const NSRDSfunc6& heatOfVapourisation,
    const NSRDSfunc0& heatCapacity,
    const NSRDSfunc0& enthalpy,
    const NSRDSfunc7& idealGasHeatCapacity,
    const NSRDSfunc4& secondVirialCoeff,
    const NSRDSfunc1& dynamicViscosity,
    const NSRDSfunc2& vapourDynamicViscosity,
    const NSRDSfunc0& thermalConductivity,
    const NSRDSfunc2& vapourThermalConductivity,
    const NSRDSfunc6& surfaceTension,
    const APIdiffCoefFunc& vapourDiffussivity
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
    D_(vapourDiffussivity)
{}



Foam::C2H6O::C2H6O(const dictionary& dict)
:
    C2H6O()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C2H6O::writeData(Ostream& os) const
{
    liquidProperties::writeData(*this, os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const C2H6O& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
