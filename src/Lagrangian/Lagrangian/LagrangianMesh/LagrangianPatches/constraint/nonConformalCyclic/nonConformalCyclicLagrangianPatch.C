/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "nonConformalCyclicLagrangianPatch.H"
#include "LagrangianFields.H"
#include "meshSearch.H"
#include "RemoteData.H"
#include "tracking.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    struct maxFirstOp
    {
        template<class Type>
        const Type& operator()(const Type& a, const Type& b) const
        {
            return a.first() > b.first() ? a : b;
        }
    };
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalCyclicLagrangianPatch, 0);

    addToRunTimeSelectionTable
    (
        LagrangianPatch,
        nonConformalCyclicLagrangianPatch,
        polyPatch
    );
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonConformalCyclicLagrangianPatch::nonConformalCyclicLagrangianPatch
(
    const polyPatch& poly,
    const LagrangianBoundaryMesh& boundaryMesh
)
:
    LagrangianPatch(poly, boundaryMesh),
    nonConformalCyclicPoly_(refCast<const nonConformalCyclicPolyPatch>(poly)),
    isNbrPatchMesh_(false),
    nPositionalErrors_(0),
    maxPositionalErrorSqr_(-vGreat),
    maxPositionalErrorReceivePosition_(point::uniform(NaN))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalCyclicLagrangianPatch::~nonConformalCyclicLagrangianPatch()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::LagrangianSubMesh&
Foam::nonConformalCyclicLagrangianPatch::mesh() const
{
    return
        boundaryMesh()
        [
            isNbrPatchMesh_
          ? nonConformalCyclicPoly_.nbrPatchIndex()
          : poly().index()
        ].LagrangianPatch::mesh();
}


void Foam::nonConformalCyclicLagrangianPatch::evaluate
(
    PstreamBuffers&,
    LagrangianMesh& mesh,
    const LagrangianInternalScalarDynamicField& fraction
) const
{
    const LagrangianSubMesh& patchMesh = this->mesh();

    const meshSearch& searchEngine = meshSearch::New(mesh.poly());

    // Sub-set the geometry and topology of the elements
    SubField<barycentric> patchCoordinates = patchMesh.sub(mesh.coordinates());
    SubField<label> patchCelli = patchMesh.sub(mesh.celli());
    SubField<label> patchFacei = patchMesh.sub(mesh.facei());
    SubField<label> patchFaceTrii = patchMesh.sub(mesh.faceTrii());

    // Sub-set the receiving information
    SubField<label> receivePatchFace =
        patchMesh.sub(mesh.receivePatchFacePtr_());
    SubField<vector> receivePosition =
        patchMesh.sub(mesh.receivePositionPtr_());

    // Get a reference to the receiving original patch
    const polyPatch& receivePp =
        nonConformalCyclicPoly_.nbrPatch().origPatch();

    // Search for the elements on the receiving side
    forAll(patchMesh, i)
    {
        patchCelli[i] =
            mesh.poly().faceOwner()[receivePatchFace[i] + receivePp.start()];

        if
        (
           !tracking::locate
            (
                searchEngine,
                receivePosition[i],
                patchCoordinates[i],
                patchCelli[i],
                patchFacei[i],
                patchFaceTrii[i],
                fraction[i + patchMesh.start()]
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
                        mesh.poly(),
                        patchCoordinates[i],
                        patchCelli[i],
                        patchFacei[i],
                        patchFaceTrii[i],
                        fraction[i + patchMesh.start()]
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
                nPositionalErrors_ ++;

                if (positionalErrorSqr > maxPositionalErrorSqr_)
                {
                    maxPositionalErrorSqr_ = positionalErrorSqr;
                    maxPositionalErrorReceivePosition_ = receivePosition[i];
                }
            }
        }
    }

    // Set the elements to be in the adjacent cell
    patchMesh.sub(mesh.states()) = LagrangianState::inCell;

    // The elements are now on the neighbour patch
    isNbrPatchMesh_ = true;

    // Invalidate the now used receiving information
    receivePatchFace = -1;
    receivePosition = point::nan;
}


void Foam::nonConformalCyclicLagrangianPatch::partition() const
{
    LagrangianPatch::partition();

    // The elements are back on this patch
    isNbrPatchMesh_ = false;

    // Check for positional errors and report if any are found
    const label nPositionalErrors =
        returnReduce(nPositionalErrors_, sumOp());

    if (nPositionalErrors != 0)
    {
        const label nTransfers =
            returnReduce(this->mesh().size(), sumOp());

        const Tuple2<scalar, vector> maxPositionalErrorInfo =
            returnReduce
            (
                Tuple2<scalar, vector>
                (
                    maxPositionalErrorSqr_,
                    maxPositionalErrorReceivePosition_
                ),
                maxFirstOp()
            );

        maxPositionalErrorSqr_ = maxPositionalErrorInfo.first();
        maxPositionalErrorReceivePosition_ = maxPositionalErrorInfo.second();

        WarningInFunction
            << nPositionalErrors << "/" << nTransfers << " elements "
            << "transferring to patch "
            << nonConformalCyclicPoly_.nbrPatch().name()
            << " were not accurately located. The largest positional error "
            << "was " << sqrt(maxPositionalErrorSqr_)
            << " from " << maxPositionalErrorReceivePosition_ << "."
            << nl;
    }

    nPositionalErrors_ = 0;
    maxPositionalErrorSqr_ = -vGreat;
    maxPositionalErrorReceivePosition_ = point::uniform(NaN);
}


// ************************************************************************* //
