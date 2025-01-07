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

#include "cyclicLagrangianPatch.H"
#include "LagrangianMesh.H"
#include "SubField.H"
#include "tracking.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicLagrangianPatch, 0);

    addToRunTimeSelectionTable
    (
        LagrangianPatch,
        cyclicLagrangianPatch,
        polyPatch
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicLagrangianPatch::cyclicLagrangianPatch
(
    const polyPatch& patch,
    const LagrangianBoundaryMesh& boundaryMesh
)
:
    LagrangianPatch(patch, boundaryMesh),
    cyclicPatch_(refCast<const cyclicPolyPatch>(patch)),
    isNbrPatchMesh_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicLagrangianPatch::~cyclicLagrangianPatch()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::LagrangianSubMesh& Foam::cyclicLagrangianPatch::mesh() const
{
    return
        boundaryMesh()
        [
            isNbrPatchMesh_
          ? cyclicPatch_.nbrPatchIndex()
          : patch().index()
        ].LagrangianPatch::mesh();
}


void Foam::cyclicLagrangianPatch::evaluate
(
    PstreamBuffers&,
    LagrangianMesh& mesh,
    const LagrangianScalarInternalDynamicField& fraction
) const
{
    const LagrangianSubMesh& patchMesh = this->mesh();

    // Sub-set the geometry and topology of the elements
    SubField<barycentric> patchCoordinates = patchMesh.sub(mesh.coordinates());
    SubField<label> patchCelli = patchMesh.sub(mesh.celli());
    SubField<label> patchFacei = patchMesh.sub(mesh.facei());
    SubField<label> patchFaceTrii = patchMesh.sub(mesh.faceTrii());

    // Cross between the cyclic patches
    forAll(patchMesh, i)
    {
        tracking::crossCyclic
        (
            cyclicPatch_,
            patchCoordinates[i],
            patchCelli[i],
            patchFacei[i],
            patchFaceTrii[i]
        );
    }

    // Set the elements to be in the adjacent cell
    patchMesh.sub(mesh.states()) = LagrangianState::inCell;

    // The elements are now on the neighbour patch
    isNbrPatchMesh_ = true;
}


void Foam::cyclicLagrangianPatch::partition() const
{
    LagrangianPatch::partition();

    // The elements are back on this patch
    isNbrPatchMesh_ = false;
}


// ************************************************************************* //
