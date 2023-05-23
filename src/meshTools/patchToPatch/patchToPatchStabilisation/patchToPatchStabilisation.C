/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "patchToPatchStabilisation.H"
#include "PatchEdgeFacePointData.H"
#include "PatchEdgeFaceWave.H"
#include "SubField.H"
#include "globalIndex.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchToPatchStabilisation, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatchStabilisation::patchToPatchStabilisation()
:
    stabilisation_(false),
    localStabilisationCells_(),
    stabilisationMapPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatchStabilisation::~patchToPatchStabilisation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::patchToPatchStabilisation::update
(
    const polyPatch& patch,
    const PackedBoolList& faceCoupleds
)
{
    // Determine whether or not stabilisation is necessary
    stabilisation_ = false;
    forAll(faceCoupleds, facei)
    {
        if (!faceCoupleds[facei])
        {
            stabilisation_ = true;
            break;
        }
    }
    reduce(stabilisation_, orOp<bool>());

    // Quick return if nothing is to be done
    if (!stabilisation_) return;

    // Construct initial edges. All edges that border a coupled face are added
    // here. The wave will propagate everywhere for just the first iteration.
    // Then most paths will end and subsequent iterations will propagate only
    // through the uncoupled faces. This is a bit odd, but it is easier than
    // doing the necessary synchronisation to determine which edges lie
    // in-between coupled and non-coupled faces.
    typedef PatchEdgeFacePointData<remote> info;
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
                        remote(Pstream::myProcNo(), facei),
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
                    remote(Pstream::myProcNo(), facei),
                    patch.faceCentres()[facei],
                    0
                );
        }
    }
    PatchEdgeFaceWave<primitivePatch, info> wave
    (
        patch.boundaryMesh().mesh(),
        patch,
        initialEdges,
        initialEdgeInfos,
        edgeInfos,
        faceInfos,
        returnReduce(patch.nEdges(), sumOp<label>())
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
                << "to a mapped cell by the stabilisation wave. This "
                << "indicates that an entire non-contiguous region of patch "
                << "lies outside of the other patch being mapped to. This is "
                << "not recoverable." << exit(FatalError);
        }
    }

    // Construct the cell to local stabilisation cell map
    const globalIndex cellGlobalIndex(patch.size());
    localStabilisationCells_.resize(patch.size());
    forAll(faceCoupleds, facei)
    {
        const remote& r = faceInfos[facei].data();
        localStabilisationCells_[facei] =
            faceCoupleds[facei]
          ? cellGlobalIndex.toGlobal(facei)
          : cellGlobalIndex.toGlobal(r.proci, r.elementi);
    }

    // Construct the distribution map, if necessary
    if (Pstream::parRun())
    {
        List<Map<label>> compactMap;
        stabilisationMapPtr_.reset
        (
            new distributionMap
            (
                cellGlobalIndex,
                localStabilisationCells_,
                compactMap
            )
        );
    }

    // Write out stabilisation connections
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
        stabilise(sfcs);

        forAll(fcs, celli)
        {
            const point& c = fcs[celli];
            if (magSqr(c - fcs[celli]) == 0) continue;
            obj.write(linePointRef(fcs[celli], c));
        }
    }
}


// ************************************************************************* //
