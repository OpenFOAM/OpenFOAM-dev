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

#include "HashTable.H"
#include "cloudAge.H"
#include "cloud.H"
#include "calculatedLagrangianPatchFields.H"
#include "zeroLagrangianFieldSources.H"
#include "internalLagrangianFieldSources.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudAge, 0);
    addToRunTimeSelectionTable(functionObject, cloudAge, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudAge::cloudAge
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudLagrangianMeshFunctionObject(name, runTime, dict),
    age_
    (
        IOobject
        (
            "age",
            runTime.name(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar(dimTime, scalar(0)),
        wordList
        (
            mesh().boundary().size(),
            calculatedLagrangianPatchScalarField::typeName
        ),
        wordList::null(),
        cloud().LagrangianModels().modelTypeFieldSourceTypes
        <
            LagrangianInjection,
            zeroLagrangianScalarFieldSource,
            LagrangianSource,
            internalLagrangianScalarFieldSource
        >()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudAge::~cloudAge()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::cloudAge::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::cloudAge::execute()
{
    return true;
}


void Foam::functionObjects::cloudAge::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    if (!final) return;

    const LagrangianSubMesh& subMesh = deltaT.mesh();

    LagrangianSubScalarSubField age(subMesh.sub(age_));

    age = age + deltaT;
}


bool Foam::functionObjects::cloudAge::write()
{
    return true;
}


bool Foam::functionObjects::cloudAge::clear()
{
    return true;
}


// ************************************************************************* //
