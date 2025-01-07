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

#include "cloudFunctionObject.H"
#include "LagrangianMeshFunctionObject.H"
#include "fvMesh.H"
#include "fvMeshFunctionObject.H"
#include "cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudFunctionObject, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudFunctionObject::cloudFunctionObject
(
    const LagrangianMeshFunctionObject& function
)
:
    function_(function),
    cloud_
    (
        function.mesh().lookupType<Foam::cloud>()
    )
{}


Foam::functionObjects::cloudFunctionObject::cloudFunctionObject
(
    const fvMeshFunctionObject& function,
    const dictionary& dict
)
:
    function_(function),
    cloud_
    (
        function.mesh().lookupObject<LagrangianMesh>
        (
            dict.lookup<word>(cloud::typeName)
        ).lookupType<Foam::cloud>()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudFunctionObject::~cloudFunctionObject()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::cloud& Foam::functionObjects::cloudFunctionObject::cloud() const
{
    return cloud<Foam::cloud>();
}


void Foam::functionObjects::cloudFunctionObject::preSolve()
{}


void Foam::functionObjects::cloudFunctionObject::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{}


void Foam::functionObjects::cloudFunctionObject::preCrossFaces
(
    const LagrangianScalarInternalDynamicField& fraction
)
{}


void Foam::functionObjects::cloudFunctionObject::preCrossFaces
(
    const LagrangianSubScalarSubField& fraction
)
{}


void Foam::functionObjects::cloudFunctionObject::postCrossFaces
(
    const LagrangianSubScalarSubField& fraction
)
{}


void Foam::functionObjects::cloudFunctionObject::postCrossFaces
(
    const LagrangianScalarInternalDynamicField& fraction
)
{}


void Foam::functionObjects::cloudFunctionObject::postSolve()
{}


// ************************************************************************* //
