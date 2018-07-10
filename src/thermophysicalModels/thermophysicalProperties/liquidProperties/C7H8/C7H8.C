/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "C7H8.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C7H8, 0);
    addToRunTimeSelectionTable(liquidProperties, C7H8,);
    addToRunTimeSelectionTable(liquidProperties, C7H8, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C7H8::C7H8()
:
    liquidProperties
    (
        92.141,
        591.79,
        4.1086e+6,
        0.31579,
        0.264,
        178.18,
        4.1009e-2,
        383.78,
        1.2008e-30,
        0.2641,
        1.8346e+4
    ),
    rho_(81.32088237, 0.27108, 591.79, 0.29889),
    pv_(83.359, -6995, -9.1635, 6.225e-06, 2.0),
    hl_(591.79, 544383.065085033, 0.3834, 0.0, 0.0, 0.0),
    Cp_
    (
        2066.83235476064,
       -8.14664481609707,
        0.0322581695445024,
       -3.01223125427334e-05,
        0.0,
        0.0
    ),
    h_
    (
       -353094.830249075,
        2066.83235476064,
       -4.07332240804853,
        0.0107527231815008,
       -7.53057813568336e-06,
        0.0
    ),
    Cpg_(630.989461803106, 3107.19440856947, 1440.6, 2059.88647833212, -650.43),
    B_
    (
        0.00191120131103418,
       -2.24970425760519,
       -482293.441573241,
       -7.62309938029759e+17,
        1.00986531511488e+20
    ),
    mu_(-13.362, 1183, 0.333, 0.0, 0.0),
    mug_(2.919e-08, 0.9648, 0.0, 0.0),
    kappa_(0.2043, -0.000239, 0.0, 0.0, 0.0, 0.0),
    kappag_(2.392e-05, 1.2694, 537, 0.0),
    sigma_(591.79, 0.06685, 1.2456, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 92.141, 28) // note: Same as nHeptane
{}


Foam::C7H8::C7H8
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


Foam::C7H8::C7H8(const dictionary& dict)
:
    C7H8()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C7H8::writeData(Ostream& os) const
{
    liquidProperties::writeData(*this, os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const C7H8& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
