/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "movingForces.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(movingForces, 0);

    addToRunTimeSelectionTable(functionObject, movingForces, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::vector Foam::functionObjects::movingForces::CofR() const
{
    return CofR_->value(time_.userTimeValue());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::movingForces::movingForces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forcesBase(name, runTime, dict)
{
    read(dict);
}


Foam::functionObjects::movingForces::movingForces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    forcesBase(name, obr, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::movingForces::~movingForces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::movingForces::read(const dictionary& dict)
{
    forcesBase::read(dict);

    // Centre of rotation for moment calculations
    CofR_ = Function1<vector>::New("CofR", dict);

    return true;
}


// ************************************************************************* //
