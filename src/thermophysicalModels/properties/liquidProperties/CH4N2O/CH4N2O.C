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

#include "CH4N2O.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CH4N2O, 0);
    addToRunTimeSelectionTable(liquidProperties, CH4N2O,);
    addToRunTimeSelectionTable(liquidProperties, CH4N2O, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CH4N2O::CH4N2O()
:
    liquidProperties
    (
        60.056,
        705.0,
        9.050e+6,
        0.218,
        0.337,
        405.85,
        9.3131e+1,
        465.0,
        1.52e-29,
        0.3449,
        4.7813e+4
    ),
    rho_(1230.006936, 0.0, 0.0, 0.0, 0.0, 0.0),
    pv_(3015.15611544, -185497.059684, -430.223621983, 0.00017405122622, 2.0),
    hl_(705.0, 2534249.0, 0.5, 0.0, 0.0, 0.0),
    Cp_(2006.46063673904, 0.0, 0.0, 0.0, 0.0, 0.0),
    h_(-6154107.41641135, 2006.46063673904, 0.0, 0.0, 0.0, 0.0),
    Cpg_(811.875582789397, 2099.04089516451, 1627.3, 1603.63660583455, 724.41),
    B_
    (
       -0.000383641934194752,
        0.447249234048222,
       -469062.208605302,
        5.5628080458239e+18,
       -2.3040162514986e+21
    ),
    mu_(-51.964, 3670.6, 5.7331, -5.3495e-29, 10.0),
    mug_(2.6986e-06, 0.498, 1257.7, -19570.0),
    kappa_(-0.4267, 0.0056903, -8.0065e-06, 1.815e-09, 0.0, 0.0),
    kappag_(6.977e-05, 1.1243, 844.9, -148850.0),
    sigma_(705.0, 1.0, 0.0, 0.0, 0.0, 0.0), // note: set to constant
    D_(147.18, 20.1, 60.056, 28.0) // note: Same as nHeptane
{}


Foam::CH4N2O::CH4N2O
(
    const liquidProperties& l,
    const NSRDSfunc0& density,
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


Foam::CH4N2O::CH4N2O(const dictionary& dict)
:
    CH4N2O()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CH4N2O::writeData(Ostream& os) const
{
    liquidProperties::writeData(*this, os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const CH4N2O& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
