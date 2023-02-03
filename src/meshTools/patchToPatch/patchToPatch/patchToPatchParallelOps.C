/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "patchToPatch.H"
#include "uindirectPrimitivePatch.H"
#include "uindirectPrimitiveOldTimePatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelListList Foam::patchToPatch::tgtPatchSendFaces
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch
) const
{
    // Get the bound boxes for the source patch. Just a single box for now.
    List<List<treeBoundBox>> srcPatchProcBbs(Pstream::nProcs());
    if (srcPatch.size())
    {
        srcPatchProcBbs[Pstream::myProcNo()] =
            List<treeBoundBox>
            (
                1,
                srcBox(srcPatch, srcPointNormals, srcPointNormals0)
            );
    }
    else
    {
        srcPatchProcBbs[Pstream::myProcNo()] = List<treeBoundBox>();
    }

    // Distribute the boxes
    Pstream::gatherList(srcPatchProcBbs);
    Pstream::scatterList(srcPatchProcBbs);

    // Send a target face to a process if it overlaps the source patch box
    // for that process
    List<DynamicList<label>> resultDyn(Pstream::nProcs());
    forAll(tgtPatch, tgtFacei)
    {
        const treeBoundBox tgtFaceBb = tgtBox(tgtPatch, tgtFacei);
        forAll(srcPatchProcBbs, proci)
        {
            forAll(srcPatchProcBbs[proci], bbi)
            {
                if (srcPatchProcBbs[proci][bbi].overlaps(tgtFaceBb))
                {
                    resultDyn[proci].append(tgtFacei);
                    break;
                }
            }
        }
    }

    // Transfer to non-dynamic storage
    labelListList result(Pstream::nProcs());
    forAll(result, proci)
    {
        result[proci].transfer(resultDyn[proci]);
    }

    return result;
}


Foam::List<Foam::remote> Foam::patchToPatch::distributePatch
(
    const distributionMap& map,
    const primitiveOldTimePatch& patch,
    autoPtr<PrimitiveOldTimePatch<faceList, pointField>>& localPatchPtr
)
{
    static const label thisProci = Pstream::myProcNo();

    // Exchange per-processor data
    List<labelList> procLocalFaceis(Pstream::nProcs());
    List<faceList> procLocalFaces(Pstream::nProcs());
    List<pointField> procLocalPoints(Pstream::nProcs());
    List<pointField> procLocalPoints0(Pstream::nProcs());
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            const labelList& sendFaceis = map.subMap()[proci];

            if (proci != thisProci && sendFaceis.size())
            {
                uindirectPrimitiveOldTimePatch subPatch
                (
                    UIndirectList<face>(patch, sendFaceis),
                    patch.points(),
                    patch.points0()
                );

                UOPstream(proci, pBufs)()
                    << sendFaceis
                    << subPatch.localFaces()
                    << subPatch.localPoints()
                    << (patch.has0() ? subPatch.localPoints0() : pointField());
            }
        }

        pBufs.finishedSends();

        // Map local data
        {
            const labelList& sendFaceis = map.subMap()[thisProci];

            uindirectPrimitiveOldTimePatch subPatch
            (
                UIndirectList<face>(patch, sendFaceis),
                patch.points(),
                patch.points0()
            );

            procLocalFaceis[thisProci] = sendFaceis;
            procLocalFaces[thisProci] = subPatch.localFaces();
            procLocalPoints[thisProci] = subPatch.localPoints();
            if (patch.has0())
            {
                procLocalPoints0[thisProci] = subPatch.localPoints0();
            }
        }

        // Receive remote data
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            if (proci != thisProci && map.constructMap()[proci].size())
            {
                UIPstream(proci, pBufs)()
                    >> procLocalFaceis[proci]
                    >> procLocalFaces[proci]
                    >> procLocalPoints[proci]
                    >> procLocalPoints0[proci];
            }
        }
    }

    // Allocate
    List<remote> localProcFaces;
    faceList localFaces;
    pointField localPoints;
    pointField localPoints0;
    {
        label nLocalFaces = 0, nLocalPoints = 0;
        forAll(procLocalFaceis, proci)
        {
            nLocalFaces += procLocalFaces[proci].size();
            nLocalPoints += procLocalPoints[proci].size();
        }
        localProcFaces.setSize(nLocalFaces);
        localFaces.setSize(nLocalFaces);
        localPoints.setSize(nLocalPoints);
        if (patch.has0())
        {
            localPoints0.setSize(nLocalPoints);
        }
    }

    // Construct the result
    label localFacei = 0, localPointi = 0;
    forAll(procLocalFaces, proci)
    {
        const labelList& fis = procLocalFaceis[proci];
        const faceList& fs = procLocalFaces[proci];
        forAll(fis, i)
        {
            localProcFaces[localFacei] = {proci, fis[i]};
            localFaces[localFacei] = face(fs[i] + localPointi);
            localFacei ++;
        }

        const pointField& ps = procLocalPoints[proci];
        const pointField& ps0 = procLocalPoints0[proci];
        forAll(ps, i)
        {
            localPoints[localPointi] = ps[i];
            if (patch.has0())
            {
                localPoints0[localPointi] = ps0[i];
            }
            localPointi ++;
        }
    }

    // Construct the local patch
    localPatchPtr.reset
    (
        patch.has0()
      ? new PrimitiveOldTimePatch<faceList, pointField>
        (
            localFaces,
            localPoints,
            localPoints0
        )
      : new PrimitiveOldTimePatch<faceList, pointField>
        (
            localFaces,
            localPoints
        )
    );

    return localProcFaces;
}


// ************************************************************************* //
