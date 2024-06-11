/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "reducedUnits.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const  Foam::scalar Foam::reducedUnits::kb = 1.3806504e-23;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::reducedUnits::calcRefValues()
{
    if
    (
        refTime_ < vSmall
     || refLength_ < vSmall
     || refMass_ < vSmall
    )
    {
        FatalErrorInFunction
            << "One of more reference values too small for floating point "
            << "calculation: "
            << "refTime_ = " << refTime_
            << ", refLength = " << refTemp_
            << ", refMass = " << refMass_
            << nl << abort(FatalError);
    }

    refEnergy_ = refLength_*refLength_*refMass_/(refTime_*refTime_);

    refTemp_ = refEnergy_ / kb;

    refForce_ = refEnergy_/refLength_;

    refVelocity_ = Foam::sqrt(refEnergy_/refMass_);

    refVolume_ = Foam::pow(refLength_,3.0);

    refPressure_ = refEnergy_/refVolume_;

    refMassDensity_ = refMass_/refVolume_;

    refNumberDensity_ = 1.0/refVolume_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reducedUnits::reducedUnits()
:
    refLength_(1e-9),
    refTime_(1e-12),
    refMass_(1.660538782e-27)
{
    calcRefValues();
}


Foam::reducedUnits::reducedUnits
(
    scalar refLength,
    scalar refTime,
    scalar refMass
)
:
    refLength_(refLength),
    refTime_(refTime),
    refMass_(refMass)
{
    calcRefValues();
}


Foam::reducedUnits::reducedUnits(const IOdictionary& reducedUnitsDict)
:
    refLength_(),
    refTime_(),
    refMass_()
{
    setRefValues(reducedUnitsDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reducedUnits::~reducedUnits()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reducedUnits::setRefValues
(
    scalar refLength,
    scalar refTime,
    scalar refMass
)
{
    refLength_ = refLength;

    refTime_ = refTime;

    refMass_ = refMass;

    calcRefValues();
}


void Foam::reducedUnits::setRefValues
(
    const IOdictionary& reducedUnitsDict
)
{
    refLength_ = reducedUnitsDict.template lookup<scalar>("refLength");

    refTime_ = reducedUnitsDict.template lookup<scalar>("refTime");

    refMass_  = reducedUnitsDict.template lookup<scalar>("refMass");

    calcRefValues();
}


// ************************************************************************* //
