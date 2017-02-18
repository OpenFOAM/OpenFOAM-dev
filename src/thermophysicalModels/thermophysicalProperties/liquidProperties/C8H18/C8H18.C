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

#include "C8H18.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C8H18, 0);
    addToRunTimeSelectionTable(liquidProperties, C8H18,);
    addToRunTimeSelectionTable(liquidProperties, C8H18, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C8H18::C8H18()
:
    liquidProperties
    (
        114.231,
        568.70,
        2.49e+6,
        0.486,
        0.256,
        216.38,
        2.1083,
        398.83,
        0.0,
        0.3996,
        1.54e+4
    ),
    rho_(61.37745861, 0.26115, 568.7, 0.28034),
    pv_(96.084, -7900.2, -11.003, 7.1802e-06, 2.0),
    hl_(568.70, 483056.263186, 0.38467, 0.0, 0.0, 0.0),
    Cp_
    (
        1968.20477803749,
       -1.63379467920267,
        0.00839448135795012,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
       -2778787.734126,
        1968.20477803749,
       -0.816897339601334,
        0.00279816045265004,
        0.0,
        0.0
    ),
    Cpg_(1186.54305748877, 3878.9820626625, 1635.6, 2673.52995246474, 746.4),
    B_
    (
        0.00239777293379205,
       -2.81394717721109,
       -585042.589139551,
       -1.11265768486663e+18,
        1.40968738783693e+20
    ),
    mu_(-20.463, 1497.4, 1.379, 0.0, 0.0),
    mug_(3.1191e-08, 0.92925, 55.092, 0.0),
    kappa_(0.2156, -0.00029483, 0.0, 0.0, 0.0, 0.0),
    kappag_(-8758, 0.8448, -27121000000.0, 0.0),
    sigma_(568.70, 0.052789, 1.2323, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 114.231, 28.0) // note: Same as nHeptane
{}


Foam::C8H18::C8H18
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


Foam::C8H18::C8H18(const dictionary& dict)
:
    C8H18()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C8H18::writeData(Ostream& os) const
{
    liquidProperties::writeData(*this, os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const C8H18& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
