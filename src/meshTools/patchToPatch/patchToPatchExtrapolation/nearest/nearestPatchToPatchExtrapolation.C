/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
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

#include "nearestPatchToPatchExtrapolation.H"
#include "distributionMap.H"
#include "PatchEdgeFacePointData.H"
#include "PatchEdgeFaceWave.H"
#include "SubField.H"
#include "globalIndex.H"
#include "OBJstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatchExtrapolations
{
    defineTypeNameAndDebug(nearest, 0);
    addToRunTimeSelectionTable(patchToPatchExtrapolation, nearest, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::patchToPatchExtrapolations::nearest::extrapolateType
(
    Field<Type>& fld
) const
{
    if (!extrapolation_) return;

    // Communicate remote field values as necessary
    tmp<Field<Type>> srcFld;
    if (Pstream::parRun())
    {
        srcFld = fld.clone();
        extrapolationMapPtr_->distribute(srcFld.ref());
    }
    else
    {
        srcFld = tmp<Field<Type>>(fld);
    }

    // Set the values in the uncoupled faces
    forAll(uncoupledFaces_, uncoupledFacei)
    {
        const label celli = uncoupledFaces_[uncoupledFacei];

        fld[celli] = srcFld()[uncoupledFaceLocalFaces_[uncoupledFacei]];
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatchExtrapolations::nearest::nearest()
:
    patchToPatchExtrapolation(),
    uncoupledFaceLocalFaces_(),
    extrapolationMapPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatchExtrapolations::nearest::~nearest()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::patchToPatchExtrapolations::nearest::update
(
    const polyPatch& patch,
    const PackedBoolList& faceCoupleds
)
{
    patchToPatchExtrapolation::update(faceCoupleds);

    // Quick return if nothing is to be done
    if (!extrapolation_) return;

    // Global patch-face addressing
    const globalIndex globalPatchFaceIndex(patch.size());

    // Construct initial edges. All edges that border a coupled face are added
    // here. The wave will propagate everywhere for just the first iteration.
    // Then most paths will end and subsequent iterations will propagate only
    // through the uncoupled faces. This is a bit odd, but it is easier than
    // doing the necessary synchronisation to determine which edges lie
    // in-between coupled and non-coupled faces.
    typedef PatchEdgeFacePointData<label> info;
    DynamicList<label> initialEdges(patch.nEdges());
    DynamicList<info> initialEdgeInfos(patch.nEdges());
    forAll(patch.edgeFaces(), edgei)
    {
        forAll(patch.edgeFaces()[edgei], edgeFacei)
        {
            const label facei = patch.edgeFaces()[edgei][edgeFacei];

            if (faceCoupleds[facei])
            {
                initialEdges.append(edgei);
                initialEdgeInfos.append
                (
                    info
                    (
                        globalPatchFaceIndex.toGlobal(facei),
                        patch.edges()[edgei].centre(patch.localPoints()),
                        0
                    )
                );
                break;
            }
        }
    }

    // Wave the information about the nearby coupled faces into the un-coupled
    // faces. Base this wave on distance to the cut face. Initialise coupled
    // faces to have a distance of zero, so that we do not waste time waving
    // into coupled regions of the patch.
    List<info> edgeInfos(patch.nEdges()), faceInfos(patch.size());
    forAll(faceCoupleds, facei)
    {
        if (faceCoupleds[facei])
        {
            faceInfos[facei] =
                info
                (
                    globalPatchFaceIndex.toGlobal(facei),
                    patch.faceCentres()[facei],
                    0
                );
        }
    }
    PatchEdgeFaceWave<primitivePatch, info> wave
    (
        patch.mesh(),
        patch,
        initialEdges,
        initialEdgeInfos,
        edgeInfos,
        faceInfos,
        returnReduce(patch.nEdges(), sumOp())
    );

    // Check that the wave connected to all un-mapped faces
    forAll(faceCoupleds, facei)
    {
        if (!faceCoupleds[facei] && !faceInfos[facei].valid(wave.data()))
        {
            FatalErrorInFunction
                << "Un-mapped face " << facei << " of patch " << patch.name()
                << " on processor " << Pstream::myProcNo() << " with centre "
                << "at " << patch.faceCentres()[facei] << " was not connected "
                << "to a mapped cell by the extrapolation wave. This "
                << "indicates that an entire non-contiguous region of patch "
                << "lies outside of the other patch being mapped to. This is "
                << "not recoverable." << exit(FatalError);
        }
    }

    // Construct the cell to local extrapolation cell map
    uncoupledFaceLocalFaces_.resize(uncoupledFaces_.size());
    forAll(uncoupledFaces_, uncoupledFacei)
    {
        const label facei = uncoupledFaces_[uncoupledFacei];

        uncoupledFaceLocalFaces_[uncoupledFacei] = faceInfos[facei].data();
    }

    // Construct the distribution map, if necessary
    if (Pstream::parRun())
    {
        List<Map<label>> compactMap;
        extrapolationMapPtr_.reset
        (
            new distributionMap
            (
                globalPatchFaceIndex,
                uncoupledFaceLocalFaces_,
                compactMap
            )
        );
    }

    // Write out connections
    if (debug)
    {
        OBJstream obj
        (
            typeName + "_" + patch.name()
          + (Pstream::parRun() ? "_proc" + name(Pstream::myProcNo()) : "")
          + "_connections.obj"
        );

        const pointField fcs(patch.faceCentres());
        pointField sfcs(fcs);
        extrapolate(sfcs);

        forAll(fcs, celli)
        {
            const point& c = fcs[celli];
            if (magSqr(c - fcs[celli]) == 0) continue;
            obj.write(linePointRef(fcs[celli], c));
        }
    }
}


#define implementExtrapolateType(Type, nullArg)                                \
    void Foam::patchToPatchExtrapolations::nearest::extrapolate                \
    (                                                                          \
        Field<Type>& fld                                                       \
    ) const                                                                    \
    {                                                                          \
        extrapolateType(fld);                                                  \
    }
FOR_ALL_FIELD_TYPES(implementExtrapolateType);
#undef implementExtrapolateType


// ************************************************************************* //
