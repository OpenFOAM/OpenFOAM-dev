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

#include "cloudVolumeFraction.H"
#include "grouped.H"
#include "shaped.H"
#include "volFields.H"
#include "LagrangiancAccumulate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudVolumeFraction, 0);
    addToRunTimeSelectionTable(functionObject, cloudVolumeFraction, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudVolumeFraction::cloudVolumeFraction
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudFvMeshFunctionObject(name, runTime, dict),
    alphaName_(cloud().mesh().name() + ":alpha")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudVolumeFraction::~cloudVolumeFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::cloudVolumeFraction::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::cloudVolumeFraction::execute()
{
    const Foam::cloud& c = cloud();

    const LagrangianScalarInternalField v
    (
        "v",
        isCloud<clouds::grouped>()
      ? cloud<clouds::grouped>().number(c.mesh())
       *cloud<clouds::shaped>().v(c.mesh())
      : cloud<clouds::shaped>().v(c.mesh())
    );

    objectRegistryFunctionObject::store
    (
        alphaName_,
        Lagrangianc::accumulate<volMesh>(v)/fvMeshFunctionObject::mesh_.V()
    );

    return true;
}


bool Foam::functionObjects::cloudVolumeFraction::write()
{
    return objectRegistryFunctionObject::writeObject(alphaName_);
}


bool Foam::functionObjects::cloudVolumeFraction::clear()
{
    return objectRegistryFunctionObject::clearObject(alphaName_);
}


// ************************************************************************* //
