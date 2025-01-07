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

#include "cloudSurfaceArea.H"
#include "calculatedLagrangianPatchFields.H"
#include "surfaceAreaLagrangianScalarFieldSource.H"
#include "shaped.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudSurfaceArea, 0);
    addToRunTimeSelectionTable(functionObject, cloudSurfaceArea, dictionary);
}
}


const Foam::word Foam::functionObjects::cloudSurfaceArea::aName_ = "a";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudSurfaceArea::cloudSurfaceArea
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudLagrangianMeshFunctionObject(name, runTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudSurfaceArea::~cloudSurfaceArea()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::cloudSurfaceArea::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::cloudSurfaceArea::execute()
{
    if (foundObject<LagrangianScalarDynamicField>(aName_))
    {
        return true;
    }

    store
    (
        LagrangianScalarDynamicField::New
        (
            aName_,
            cloud<clouds::shaped>().a(mesh())(),
            wordList
            (
                cloud().mesh().boundary().size(),
                calculatedLagrangianPatchScalarField::typeName
            ),
            wordList::null(),
            cloud().LagrangianModels().modelTypeFieldSourceTypes
            <
                LagrangianInjection,
                surfaceAreaLagrangianScalarFieldSource
            >()
        )
    );

    return true;
}


void Foam::functionObjects::cloudSurfaceArea::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    LagrangianSubScalarSubField a
    (
        subMesh.sub
        (
            lookupObjectRef<LagrangianScalarDynamicField>
            (
                aName_
            )
        )
    );

    a = cloud<clouds::shaped>().a(subMesh);
}


bool Foam::functionObjects::cloudSurfaceArea::write()
{
    return writeObject(aName_);
}


bool Foam::functionObjects::cloudSurfaceArea::clear()
{
    return clearObject(aName_);
}


// ************************************************************************* //
