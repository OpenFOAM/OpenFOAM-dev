/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "patchRegionSplit.H"
#include "uindirectPrimitivePatch.H"
#include "PatchEdgeFaceWave.H"
#include "globalIndex.H"
#include "patchEdgeFaceRegion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchRegionSplit, 0);
}


// * * * * * * * * * * * * * * Private Constructors  * * * * * * * * * * * * //

template<class Patch>
Foam::patchRegionSplit::patchRegionSplit
(
    const polyMesh& mesh,
    const Patch& patch
)
:
    regionSplitBase(patch.size())
{
    // Create global indices for the edges
    const globalIndex globalPatchEdgeIndex(patch.nEdges());

    // Initialise wave data
    List<patchEdgeFaceRegion> edgeData(patch.nEdges()), faceData(patch.size());

    // Seed all edges with their global index
    labelList seedEdges = identityMap(patch.nEdges());
    List<patchEdgeFaceRegion> seedEdgesData(patch.nEdges());
    forAll(seedEdges, patchEdgei)
    {
        seedEdgesData[patchEdgei] = globalPatchEdgeIndex.toGlobal(patchEdgei);
    }

    // Propagate inwards so that every face data now contains the lowest global
    // edge index in the contiguous region
    PatchEdgeFaceWave<Patch, patchEdgeFaceRegion> deltaCalc
    (
        mesh,
        patch,
        seedEdges,
        seedEdgesData,
        edgeData,
        faceData,
        returnReduce(patch.size(), sumOp()),
        PatchEdgeFaceWave<Patch, patchEdgeFaceRegion>::defaultTrackingData_
    );

    // Unpack indices into the region index list
    forAll(patch, patchFacei)
    {
        this->operator[](patchFacei) = faceData[patchFacei].region();
    }

    // Compact the region indices
    nRegions_ =
        Pstream::parRun()
      ? compactGlobalRegionSplit(globalPatchEdgeIndex, *this)
      : compactLocalRegionSplit(*this);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchRegionSplit::patchRegionSplit
(
    const polyMesh& mesh,
    const labelList& faces
)
:
    patchRegionSplit
    (
        mesh,
        uindirectPrimitivePatch
        (
            UIndirectList<face>(mesh.faces(), faces),
            mesh.points()
        )
    )
{}


Foam::patchRegionSplit::patchRegionSplit(const polyPatch& patch)
:
    patchRegionSplit(patch.mesh(), patch)
{}


// ************************************************************************* //
