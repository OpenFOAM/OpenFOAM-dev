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

#include "cloudGravitationalPotentialEnergy.H"
#include "lookupUniformDimensionedField.H"
#include "shaped.H"
#include "massive.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudGravitationalPotentialEnergy, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        cloudGravitationalPotentialEnergy,
        dictionary
    );
}
}


const Foam::word
    Foam::functionObjects::cloudGravitationalPotentialEnergy::GPEName_ = "GPE";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudGravitationalPotentialEnergy::
cloudGravitationalPotentialEnergy
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudLagrangianMeshFunctionObject(name, runTime, dict),
    g(lookupUniformDimensionedField<vector>(mesh().time(), "g"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudGravitationalPotentialEnergy::
~cloudGravitationalPotentialEnergy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList
Foam::functionObjects::cloudGravitationalPotentialEnergy::fields() const
{
    return wordList::null();
}


void Foam::functionObjects::cloudGravitationalPotentialEnergy::preSolve()
{
    objectRegistryFunctionObject::clearObject(GPEName_);
}


bool Foam::functionObjects::cloudGravitationalPotentialEnergy::execute()
{
    const Foam::cloud& c = cloud();

    objectRegistryFunctionObject::store
    (
        GPEName_,
        isCloud<clouds::massive>()
      ? cloud<clouds::massive>().m(c.mesh())*(-g & c.mesh().position())
      : cloud<clouds::shaped>().v(c.mesh())*(-g & c.mesh().position())
    );

    return true;
}


bool Foam::functionObjects::cloudGravitationalPotentialEnergy::write()
{
    return objectRegistryFunctionObject::writeObject(GPEName_);
}


bool Foam::functionObjects::cloudGravitationalPotentialEnergy::clear()
{
    return objectRegistryFunctionObject::clearObject(GPEName_);
}


// ************************************************************************* //
