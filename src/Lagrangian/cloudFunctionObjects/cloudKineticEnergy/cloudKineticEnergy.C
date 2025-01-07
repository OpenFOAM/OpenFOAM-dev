/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "cloudKineticEnergy.H"
#include "shaped.H"
#include "massive.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudKineticEnergy, 0);
    addToRunTimeSelectionTable(functionObject, cloudKineticEnergy, dictionary);
}
}


const Foam::word Foam::functionObjects::cloudKineticEnergy::KEName_ = "KE";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudKineticEnergy::cloudKineticEnergy
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudLagrangianMeshFunctionObject(name, runTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudKineticEnergy::~cloudKineticEnergy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::cloudKineticEnergy::fields() const
{
    return wordList::null();
}


void Foam::functionObjects::cloudKineticEnergy::preSolve()
{
    objectRegistryFunctionObject::clearObject(KEName_);
}


bool Foam::functionObjects::cloudKineticEnergy::execute()
{
    const Foam::cloud& c = cloud();

    objectRegistryFunctionObject::store
    (
        KEName_,
        isCloud<clouds::massive>()
      ? cloud<clouds::massive>().m(c.mesh())*magSqr(c.U(c.mesh())())/2
      : cloud<clouds::shaped>().v(c.mesh())*magSqr(c.U(c.mesh())())/2
    );

    return true;
}


bool Foam::functionObjects::cloudKineticEnergy::write()
{
    return objectRegistryFunctionObject::writeObject(KEName_);
}


bool Foam::functionObjects::cloudKineticEnergy::clear()
{
    return objectRegistryFunctionObject::clearObject(KEName_);
}


// ************************************************************************* //
