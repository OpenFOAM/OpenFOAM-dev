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

#include "processorLagrangianPatch.H"
#include "LagrangianFields.H"
#include "tracking.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorLagrangianPatch, 0);

    addToRunTimeSelectionTable
    (
        LagrangianPatch,
        processorLagrangianPatch,
        polyPatch
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorLagrangianPatch::processorLagrangianPatch
(
    const polyPatch& patch,
    const LagrangianBoundaryMesh& boundaryMesh
)
:
    LagrangianPatch(patch, boundaryMesh),
    processorPatch_(refCast<const processorPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorLagrangianPatch::~processorLagrangianPatch()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::LagrangianSubMesh& Foam::processorLagrangianPatch::mesh() const
{
    return
        receiveMeshPtr_.valid()
      ? receiveMeshPtr_()
      : LagrangianPatch::mesh();
}


void Foam::processorLagrangianPatch::initEvaluate
(
    PstreamBuffers& pBufs,
    LagrangianMesh& mesh,
    const LagrangianScalarInternalDynamicField& fraction
) const
{
    const LagrangianSubMesh& patchMesh = this->mesh();

    // Sub-set the geometry and topology of the elements
    SubField<barycentric> sendCoordinates = patchMesh.sub(mesh.coordinates());
    labelField sendCelli(patchMesh.sub(mesh.celli()));
    labelField sendFacei(patchMesh.sub(mesh.facei()));
    SubField<label> sendFaceTrii = patchMesh.sub(mesh.faceTrii());
    SubField<scalar> sendFraction = patchMesh.sub(fraction.primitiveField());

    // Pre-communication steps
    forAll(patchMesh, i)
    {
        tracking::inProcessor
        (
            processorPatch_,
            sendCelli[i],
            sendFacei[i]
        );
    }

    // Send
    UOPstream(processorPatch_.neighbProcNo(), pBufs)()
        << sendCoordinates
        << sendCelli
        << sendFacei
        << sendFaceTrii
        << sendFraction;

    // Remove the sent elements
    patchMesh.sub(mesh.states()) = LagrangianState::toBeRemoved;
}


void Foam::processorLagrangianPatch::evaluate
(
    PstreamBuffers& pBufs,
    LagrangianMesh& mesh,
    const LagrangianScalarInternalDynamicField& fraction
) const
{
    // Receive
    UIPstream uips(processorPatch_.neighbProcNo(), pBufs);
    barycentricField receiveCoordinates(uips);
    labelField receiveCelli(uips);
    labelField receiveFacei(uips);
    labelField receiveFaceTrii(uips);
    scalarField receiveFraction(uips);

    // Post-communication steps
    forAll(receiveCoordinates, i)
    {
        tracking::outProcessor
        (
            processorPatch_,
            receiveCoordinates[i],
            receiveCelli[i],
            receiveFacei[i],
            receiveFaceTrii[i]
        );
    }

    // Insert the received elements and store the sub-mesh in which they reside
    receiveMeshPtr_.set
    (
        new LagrangianSubMesh
        (
            mesh.append
            (
                receiveCoordinates,
                receiveCelli,
                receiveFacei,
                receiveFaceTrii
            )
        )
    );

    // Set the step fractions of the elements
    mesh.appendSpecifiedField<scalar, LagrangianInternalDynamicField>
    (
        receiveMeshPtr_(),
        const_cast<LagrangianScalarInternalDynamicField&>(fraction),
        receiveFraction
    );

    // Set the elements to be in the adjacent cell
    receiveMeshPtr_().sub(mesh.states()) = LagrangianState::inCell;
}


void Foam::processorLagrangianPatch::partition() const
{
    LagrangianPatch::partition();

    // The elements are back on this patch, and the receive mesh is now out of
    // date and no longer needed
    receiveMeshPtr_.clear();
}


// ************************************************************************* //
