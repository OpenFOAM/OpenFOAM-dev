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

#include "forces.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forces, 0);

    addToRunTimeSelectionTable(functionObject, forces, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::vector Foam::functionObjects::forces::CofR() const
{
    return CofR_;
}


void Foam::functionObjects::forces::writeCoRValueHeader(Ostream& file)
{
    writeHeaderValue(file, "CofR", CofR());
}


void Foam::functionObjects::forces::writeCoRHeader(Ostream& file)
{}


void Foam::functionObjects::forces::writeCofR(Ostream& file)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forces::forces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forcesBase(name, runTime, dict),
    CofR_(vector::zero)
{
    read(dict);
}


Foam::functionObjects::forces::forces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    forcesBase(name, obr, dict),
    CofR_(vector::zero)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::forces::~forces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forces::read(const dictionary& dict)
{
    forcesBase::read(dict);

    // Centre of rotation for moment calculations
    CofR_ = dict.lookup("CofR");

    return true;
}


// ************************************************************************* //
