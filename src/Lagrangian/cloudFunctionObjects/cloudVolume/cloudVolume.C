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

#include "cloudVolume.H"
#include "calculatedLagrangianPatchFields.H"
#include "volumeLagrangianScalarFieldSource.H"
#include "shaped.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudVolume, 0);
    addToRunTimeSelectionTable(functionObject, cloudVolume, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudVolume::cloudVolume
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudLagrangianMeshFunctionObject(name, runTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudVolume::~cloudVolume()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::cloudVolume::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::cloudVolume::execute()
{
    if (foundObject<LagrangianScalarDynamicField>(clouds::shaped::vName))
    {
        return true;
    }

    store
    (
        LagrangianScalarDynamicField::New
        (
            clouds::shaped::vName,
            cloud<clouds::shaped>().v(mesh())(),
            wordList
            (
                cloud().mesh().boundary().size(),
                calculatedLagrangianPatchScalarField::typeName
            ),
            wordList::null(),
            cloud().LagrangianModels().modelTypeFieldSourceTypes
            <
                LagrangianInjection,
                volumeLagrangianScalarFieldSource
            >()
        )
    );

    return true;
}


void Foam::functionObjects::cloudVolume::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    LagrangianSubScalarSubField v
    (
        subMesh.sub
        (
            lookupObjectRef<LagrangianScalarDynamicField>
            (
                clouds::shaped::vName
            )
        )
    );

    v = cloud<clouds::shaped>().v(subMesh);
}


bool Foam::functionObjects::cloudVolume::write()
{
    return writeObject(clouds::shaped::vName);
}


bool Foam::functionObjects::cloudVolume::clear()
{
    return clearObject(clouds::shaped::vName);
}


// ************************************************************************* //
