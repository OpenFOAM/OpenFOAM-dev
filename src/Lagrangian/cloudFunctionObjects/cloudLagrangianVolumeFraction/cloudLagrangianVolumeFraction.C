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

#include "cloudLagrangianVolumeFraction.H"
#include "shaped.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudLagrangianVolumeFraction, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        cloudLagrangianVolumeFraction,
        dictionary
    );
}
}


const Foam::word
    Foam::functionObjects::cloudLagrangianVolumeFraction::alphaName_ = "alpha";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudLagrangianVolumeFraction::
cloudLagrangianVolumeFraction
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudLagrangianMeshFunctionObject(name, runTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudLagrangianVolumeFraction::
~cloudLagrangianVolumeFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList
Foam::functionObjects::cloudLagrangianVolumeFraction::fields() const
{
    return wordList::null();
}


void Foam::functionObjects::cloudLagrangianVolumeFraction::preSolve()
{
    objectRegistryFunctionObject::clearObject(alphaName_);
}


bool Foam::functionObjects::cloudLagrangianVolumeFraction::execute()
{
    const Foam::cloud& c = cloud();

    objectRegistryFunctionObject::store
    (
        alphaName_,
        cloud<clouds::shaped>().alpha(c.mesh())
    );

    return true;
}


bool Foam::functionObjects::cloudLagrangianVolumeFraction::write()
{
    return objectRegistryFunctionObject::writeObject(alphaName_);
}


bool Foam::functionObjects::cloudLagrangianVolumeFraction::clear()
{
    return objectRegistryFunctionObject::clearObject(alphaName_);
}


// ************************************************************************* //
