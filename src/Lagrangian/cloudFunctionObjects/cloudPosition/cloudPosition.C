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

#include "cloudPosition.H"
#include "calculatedLagrangianPatchFields.H"
#include "positionLagrangianVectorFieldSource.H"
#include "zeroLagrangianFieldSources.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudPosition, 0);
    addToRunTimeSelectionTable(functionObject, cloudPosition, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudPosition::cloudPosition
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudLagrangianMeshFunctionObject(name, runTime, dict),
    position_
    (
        IOobject
        (
            "position",
            runTime.name(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh().position()(),
        wordList
        (
            mesh().boundary().size(),
            calculatedLagrangianPatchVectorField::typeName
        ),
        wordList::null(),
        cloud().LagrangianModels().modelTypeFieldSourceTypes
        <
            LagrangianInjection,
            positionLagrangianVectorFieldSource,
            LagrangianSource,
            zeroLagrangianScalarFieldSource
        >()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudPosition::~cloudPosition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::cloudPosition::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::cloudPosition::execute()
{
    return true;
}


void Foam::functionObjects::cloudPosition::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    if (!final) return;

    const LagrangianSubMesh& subMesh = deltaT.mesh();

    LagrangianSubVectorSubField position(subMesh.sub(position_));

    position = mesh().position(subMesh);
}


void Foam::functionObjects::cloudPosition::postCrossFaces
(
    const LagrangianSubScalarSubField& fraction
)
{
    const LagrangianSubMesh& subMesh = fraction.mesh();

    LagrangianSubVectorSubField position(subMesh.sub(position_));

    position = mesh().position(subMesh);
}


bool Foam::functionObjects::cloudPosition::write()
{
    return true;
}


bool Foam::functionObjects::cloudPosition::clear()
{
    return true;
}


// ************************************************************************* //
