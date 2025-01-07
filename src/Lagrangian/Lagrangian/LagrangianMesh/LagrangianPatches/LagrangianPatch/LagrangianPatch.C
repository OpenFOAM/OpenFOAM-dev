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

#include "LagrangianPatch.H"
#include "LagrangianMesh.H"
#include "LagrangianFields.H"
#include "tracking.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LagrangianPatch, 0);

    defineRunTimeSelectionTable(LagrangianPatch, polyPatch);

    addToRunTimeSelectionTable(LagrangianPatch, LagrangianPatch, polyPatch);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LagrangianPatch::LagrangianPatch
(
    const polyPatch& patch,
    const LagrangianBoundaryMesh& boundaryMesh
)
:
    patch_(patch),
    boundaryMesh_(boundaryMesh),
    meshPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::LagrangianPatch> Foam::LagrangianPatch::New
(
    const word& type,
    const polyPatch& patch,
    const LagrangianBoundaryMesh& boundaryMesh
)
{
    polyPatchConstructorTable::iterator cstrIter =
        polyPatchConstructorTablePtr_->find(type);

    if (cstrIter == polyPatchConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown LagrangianPatch type " << type
            << nl << nl
            << "Valid LagrangianPatch types are :" << endl
            << polyPatchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<LagrangianPatch>(cstrIter()(patch, boundaryMesh));
}


Foam::autoPtr<Foam::LagrangianPatch> Foam::LagrangianPatch::New
(
    const polyPatch& patch,
    const LagrangianBoundaryMesh& boundaryMesh
)
{
    return New(patch.type(), patch, boundaryMesh);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LagrangianPatch::~LagrangianPatch()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::LagrangianSubMesh& Foam::LagrangianPatch::mesh() const
{
    if (!meshPtr_.valid())
    {
        FatalErrorInFunction
            << "Lagrangian patch mesh requested out of scope of Lagrangian "
            << "mesh changes" << exit(FatalError);
    }

    return meshPtr_;
}


const Foam::objectRegistry& Foam::LagrangianPatch::db() const
{
    return boundaryMesh().mesh();
}


void Foam::LagrangianPatch::initEvaluate
(
    PstreamBuffers&,
    LagrangianMesh&,
    const LagrangianScalarInternalDynamicField& fraction
) const
{}


void Foam::LagrangianPatch::evaluate
(
    PstreamBuffers&,
    LagrangianMesh& mesh,
    const LagrangianScalarInternalDynamicField& fraction
) const
{
    this->mesh().sub(mesh.states()) = LagrangianState::toBeRemoved;
}


void Foam::LagrangianPatch::partition() const
{
    meshPtr_.reset
    (
        boundaryMesh_.mesh().changing()
      ? new LagrangianSubMesh
        (
            boundaryMesh_.mesh().sub
            (
                static_cast<LagrangianGroup>
                (
                    static_cast<label>(LagrangianGroup::onPatchZero)
                  + patch_.index()
                )
            )
        )
      : nullptr
    );
}


// ************************************************************************* //
