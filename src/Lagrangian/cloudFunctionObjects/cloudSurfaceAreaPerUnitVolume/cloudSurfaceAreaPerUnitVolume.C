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

#include "cloudSurfaceAreaPerUnitVolume.H"
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
    defineTypeNameAndDebug(cloudSurfaceAreaPerUnitVolume, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        cloudSurfaceAreaPerUnitVolume,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudSurfaceAreaPerUnitVolume::
cloudSurfaceAreaPerUnitVolume
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudFvMeshFunctionObject(name, runTime, dict),
    AvName_(cloud().mesh().name() + ":Av")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudSurfaceAreaPerUnitVolume::
~cloudSurfaceAreaPerUnitVolume()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList
Foam::functionObjects::cloudSurfaceAreaPerUnitVolume::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::cloudSurfaceAreaPerUnitVolume::execute()
{
    const Foam::cloud& c = cloud();

    const LagrangianScalarInternalField a
    (
        "a",
        isCloud<clouds::grouped>()
      ? cloud<clouds::grouped>().number(c.mesh())
       *cloud<clouds::shaped>().a(c.mesh())
      : cloud<clouds::shaped>().a(c.mesh())
    );

    objectRegistryFunctionObject::store
    (
        AvName_,
        Lagrangianc::accumulate<volMesh>(a)/fvMeshFunctionObject::mesh_.V()
    );

    return true;
}


bool Foam::functionObjects::cloudSurfaceAreaPerUnitVolume::write()
{
    return objectRegistryFunctionObject::writeObject(AvName_);
}


bool Foam::functionObjects::cloudSurfaceAreaPerUnitVolume::clear()
{
    return objectRegistryFunctionObject::clearObject(AvName_);
}


// ************************************************************************* //
