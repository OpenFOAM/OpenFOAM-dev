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

#include "nonConformalProcessorCyclicLagrangianPatch.H"
#include "LagrangianFields.H"
#include "RemoteData.H"
#include "tracking.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalProcessorCyclicLagrangianPatch, 0);

    addToRunTimeSelectionTable
    (
        LagrangianPatch,
        nonConformalProcessorCyclicLagrangianPatch,
        polyPatch
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonConformalProcessorCyclicLagrangianPatch::
nonConformalProcessorCyclicLagrangianPatch
(
    const polyPatch& patch,
    const LagrangianBoundaryMesh& boundaryMesh
)
:
    processorCyclicLagrangianPatch(patch, boundaryMesh),
    nonConformalProcessorCyclicPatch_
    (
        refCast<const nonConformalProcessorCyclicPolyPatch>(patch)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalProcessorCyclicLagrangianPatch::
~nonConformalProcessorCyclicLagrangianPatch()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::nonConformalProcessorCyclicLagrangianPatch::initEvaluate
(
    PstreamBuffers& pBufs,
    LagrangianMesh& mesh,
    const LagrangianScalarInternalDynamicField& fraction
) const
{
    const LagrangianSubMesh& patchMesh = this->mesh();

    // Sub-set the step fractions
    SubField<scalar> sendFraction = patchMesh.sub(fraction.primitiveField());

    // Sub-set the receiving information
    SubField<label> receivePatchFace =
        patchMesh.sub(mesh.receivePatchFacePtr_());
    SubField<point> receivePosition =
        patchMesh.sub(mesh.receivePositionPtr_());

    // Send
    UOPstream(nonConformalProcessorCyclicPatch_.neighbProcNo(), pBufs)()
        << receivePatchFace
        << receivePosition
        << sendFraction;

    // Remove the sent elements
    patchMesh.sub(mesh.states()) = LagrangianState::toBeRemoved;

    // Invalidate the now used receiving information
    receivePatchFace = -1;
    receivePosition = point::nan;
}


void Foam::nonConformalProcessorCyclicLagrangianPatch::evaluate
(
    PstreamBuffers& pBufs,
    LagrangianMesh& mesh,
    const LagrangianScalarInternalDynamicField& fraction
) const
{
    // Receive
    UIPstream uips(nonConformalProcessorCyclicPatch_.neighbProcNo(), pBufs);
    labelField receivePatchFace(uips);
    pointField receivePosition(uips);
    scalarField receiveFraction(uips);

    // Get a reference to the receiving original patch
    const polyPatch& receivePp =
        nonConformalProcessorCyclicPatch_.referPatch().origPatch();

    // Search for the elements on the receiving side
    barycentricField receiveCoordinates(receivePatchFace.size());
    labelField receiveCelli(receivePatchFace.size());
    labelField receiveFacei(receivePatchFace.size());
    labelField receiveFaceTrii(receivePatchFace.size());
    forAll(receivePatchFace, i)
    {
        receiveCelli[i] =
            mesh.mesh().faceOwner()[receivePatchFace[i] + receivePp.start()];

        if
        (
           !tracking::locate
            (
                mesh.mesh(),
                receivePosition[i],
                receiveCoordinates[i],
                receiveCelli[i],
                receiveFacei[i],
                receiveFaceTrii[i],
                receiveFraction[i]
            )
        )
        {
            // The search hit a boundary. Compute the positional error.
            const scalar positionalErrorSqr =
                magSqr
                (
                    receivePosition[i]
                  - tracking::position
                    (
                        mesh.mesh(),
                        receiveCoordinates[i],
                        receiveCelli[i],
                        receiveFacei[i],
                        receiveFaceTrii[i],
                        receiveFraction[i]
                    )
                );

            // If the positional error is significant, relative to the size of
            // the receiving face, then register as a failure
            if
            (
                sqr(positionalErrorSqr)
              > sqr(small)*magSqr(receivePp.faceAreas()[receivePatchFace[i]])
            )
            {
                referPatch().nPositionalErrors_ ++;

                if (positionalErrorSqr > referPatch().maxPositionalErrorSqr_)
                {
                    referPatch().maxPositionalErrorSqr_ = positionalErrorSqr;
                    referPatch().maxPositionalErrorReceivePosition_ =
                        receivePosition[i];
                }
            }
        }
    }

    // Insert the received elements
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


// ************************************************************************* //
