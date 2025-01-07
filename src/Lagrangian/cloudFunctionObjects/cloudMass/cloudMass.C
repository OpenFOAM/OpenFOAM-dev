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

#include "cloudMass.H"
#include "calculatedLagrangianPatchFields.H"
#include "massLagrangianScalarFieldSource.H"
#include "massive.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudMass, 0);
    addToRunTimeSelectionTable(functionObject, cloudMass, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudMass::cloudMass
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudLagrangianMeshFunctionObject(name, runTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudMass::~cloudMass()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::cloudMass::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::cloudMass::execute()
{
    if (foundObject<LagrangianScalarDynamicField>(clouds::massive::mName))
    {
        return true;
    }

    store
    (
        LagrangianScalarDynamicField::New
        (
            clouds::massive::mName,
            cloud<clouds::massive>().m(mesh())(),
            wordList
            (
                cloud().mesh().boundary().size(),
                calculatedLagrangianPatchScalarField::typeName
            ),
            wordList::null(),
            cloud().LagrangianModels().modelTypeFieldSourceTypes
            <
                LagrangianInjection,
                massLagrangianScalarFieldSource
            >()
        )
    );

    return true;
}


void Foam::functionObjects::cloudMass::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    LagrangianSubScalarSubField m
    (
        subMesh.sub
        (
            lookupObjectRef<LagrangianScalarDynamicField>
            (
                clouds::massive::mName
            )
        )
    );

    m = cloud<clouds::massive>().m(subMesh);
}


bool Foam::functionObjects::cloudMass::write()
{
    return writeObject(clouds::massive::mName);
}


bool Foam::functionObjects::cloudMass::clear()
{
    return clearObject(clouds::massive::mName);
}


// ************************************************************************* //
