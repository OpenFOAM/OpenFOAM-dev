/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
#include "treeBoundBoxList.H"
#include "uindirectPrimitivePatch.H"
#include "uindirectPrimitiveOldTimePatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchToPatch::calcSingleProcess
(
    const primitiveOldTimePatch& srcPatch,
    const primitiveOldTimePatch& tgtPatch
)
{
    singleProcess_ = 0;

    if (Pstream::parRun())
    {
        boolList procHasFaces(Pstream::nProcs(), false);

        if ((srcPatch.size() > 0) || (tgtPatch.size() > 0))
        {
            procHasFaces[Pstream::myProcNo()] = true;
        }

        Pstream::gatherList(procHasFaces);
        Pstream::scatterList(procHasFaces);

        const label nProcsHaveFaces = count(procHasFaces, true);

        if (nProcsHaveFaces == 0)
        {
            singleProcess_ = 0;
        }
        if (nProcsHaveFaces == 1)
        {
            singleProcess_ = findIndex(procHasFaces, true);
        }
        if (nProcsHaveFaces > 1)
        {
            singleProcess_ = -1;
        }
    }
}


Foam::labelListList Foam::patchToPatch::tgtPatchSendFaces
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch
) const
{
    // Get the bound boxes for the source patch. Just a single box for now.
    List<treeBoundBoxList> srcPatchProcBbs(Pstream::nProcs());
    if (srcPatch.size())
    {
        srcPatchProcBbs[Pstream::myProcNo()] =
            treeBoundBoxList
            (
                1,
                srcBox(srcPatch, srcPointNormals, srcPointNormals0)
            );
    }
    else
    {
        srcPatchProcBbs[Pstream::myProcNo()] = treeBoundBoxList();
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


Foam::labelListList Foam::patchToPatch::srcPatchSendFaces() const
{
    // Send a source face to a proc if target face on that proc intersects it
    List<labelHashSet> resultDyn(Pstream::nProcs());
    forAll(tgtLocalSrcFaces_, tgtFacei)
    {
        const label tgtProci = localTgtProcFacesPtr_()[tgtFacei].proci;

        forAll(tgtLocalSrcFaces_[tgtFacei], i)
        {
            const label srcFacei = tgtLocalSrcFaces_[tgtFacei][i];

            resultDyn[tgtProci].insert(srcFacei);
        }
    }

    // Transfer to non-dynamic storage
    labelListList result(Pstream::nProcs());
    forAll(result, proci)
    {
        result[proci] = resultDyn[proci].toc();
    }

    return result;
}


Foam::distributionMap Foam::patchToPatch::patchDistributionMap
(
    labelListList&& sendFaces
) const
{
    // Figure out how many target faces are to be received
    labelList nReceiveFaces(Pstream::nProcs());
    {
        labelListList nSendFaces(Pstream::nProcs());

        nSendFaces[Pstream::myProcNo()].setSize(Pstream::nProcs());
        forAll(sendFaces, proci)
        {
            nSendFaces[Pstream::myProcNo()][proci] =
                sendFaces[proci].size();
        }

        Pstream::gatherList(nSendFaces);
        Pstream::scatterList(nSendFaces);

        forAll(sendFaces, proci)
        {
            nReceiveFaces[proci] =
                nSendFaces[proci][Pstream::myProcNo()];
        }
    }

    // Determine order of receiving
    labelListList receiveFaces(Pstream::nProcs());

    // Local faces first
    receiveFaces[Pstream::myProcNo()] =
        identity(sendFaces[Pstream::myProcNo()].size());

    // Remote faces next
    label localFacei = receiveFaces[Pstream::myProcNo()].size();
    forAll(receiveFaces, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            const label n = nReceiveFaces[proci];
            receiveFaces[proci].setSize(n);
            for (label i = 0; i < n; i ++)
            {
                receiveFaces[proci][i] = localFacei ++;
            }
        }
    }

    // Construct and return the map
    return
        distributionMap
        (
            localFacei,
            move(sendFaces),
            move(receiveFaces)
        );
}


void Foam::patchToPatch::distributePatch
(
    const distributionMap& map,
    List<procFace>& localProcFaces
) const
{
    static const label thisProci = Pstream::myProcNo();

    // Exchange per-processor data
    List<labelList> procLocalFaceis(Pstream::nProcs());
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            const labelList& sendFaces = map.subMap()[proci];

            if (proci != thisProci && sendFaces.size())
            {
                UOPstream(proci, pBufs)() << sendFaces;
            }
        }

        pBufs.finishedSends();

        // Map local data
        {
            const labelList& sendFaces = map.subMap()[thisProci];

            procLocalFaceis[thisProci] = sendFaces;
        }

        // Receive remote data
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            const labelList& receiveNewFaces = map.constructMap()[proci];

            if (proci != thisProci && receiveNewFaces.size())
            {
                UIPstream(proci, pBufs)() >> procLocalFaceis[proci];
            }
        }
    }

    // Allocate
    {
        label nLocalFaces = 0;
        forAll(procLocalFaceis, proci)
        {
            nLocalFaces += procLocalFaceis[proci].size();
        }
        localProcFaces.setSize(nLocalFaces);
    }

    // Renumber and flatten
    label localTgtFacei = 0;

    // Local data first
    {
        const labelList& fis = procLocalFaceis[thisProci];
        forAll(fis, i)
        {
            localProcFaces[localTgtFacei] = {thisProci, fis[i]};
            localTgtFacei ++;
        }
    }

    // Remote data after
    forAll(procLocalFaceis, proci)
    {
        if (proci != thisProci)
        {
            const labelList& fis = procLocalFaceis[proci];
            forAll(fis, i)
            {
                localProcFaces[localTgtFacei] = {proci, fis[i]};
                localTgtFacei ++;
            }
        }
    }
}


Foam::PrimitiveOldTimePatch<Foam::faceList, Foam::pointField>
Foam::patchToPatch::distributePatch
(
    const distributionMap& map,
    const primitiveOldTimePatch& patch,
    List<procFace>& localProcFaces
) const
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
            const labelList& sendTgtFaces = map.subMap()[proci];

            if (proci != thisProci && sendTgtFaces.size())
            {
                uindirectPrimitiveOldTimePatch subPatch
                (
                    UIndirectList<face>(patch, sendTgtFaces),
                    patch.points(),
                    patch.points0()
                );

                UOPstream(proci, pBufs)()
                    << sendTgtFaces
                    << subPatch.localFaces()
                    << subPatch.localPoints()
                    << (patch.has0() ? subPatch.localPoints0() : pointField());
            }
        }

        pBufs.finishedSends();

        // Map local data
        {
            const labelList& sendTgtFaces = map.subMap()[thisProci];

            uindirectPrimitiveOldTimePatch subPatch
            (
                UIndirectList<face>(patch, sendTgtFaces),
                patch.points(),
                patch.points0()
            );

            procLocalFaceis[thisProci] = sendTgtFaces;
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
            const labelList& receiveNewTgtFaces = map.constructMap()[proci];

            if (proci != thisProci && receiveNewTgtFaces.size())
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
    faceList localTgtFaces;
    pointField localTgtPoints;
    pointField localTgtPoints0;
    {
        label nLocalFaces = 0, nLocalPoints = 0;
        forAll(procLocalFaceis, proci)
        {
            nLocalFaces += procLocalFaces[proci].size();
            nLocalPoints += procLocalPoints[proci].size();
        }
        localProcFaces.setSize(nLocalFaces);
        localTgtFaces.setSize(nLocalFaces);
        localTgtPoints.setSize(nLocalPoints);
        if (patch.has0())
        {
            localTgtPoints0.setSize(nLocalPoints);
        }
    }

    // Renumber and flatten
    label localTgtFacei = 0, localTgtPointi = 0;

    // Local data first
    {
        const labelList& fis = procLocalFaceis[thisProci];
        const faceList& fs = procLocalFaces[thisProci];
        forAll(fis, i)
        {
            localProcFaces[localTgtFacei] = {thisProci, fis[i]};
            localTgtFaces[localTgtFacei] = face(fs[i] + localTgtPointi);
            localTgtFacei ++;
        }

        const pointField& ps = procLocalPoints[thisProci];
        const pointField& ps0 = procLocalPoints0[thisProci];
        forAll(ps, i)
        {
            localTgtPoints[localTgtPointi] = ps[i];
            if (patch.has0())
            {
                localTgtPoints0[localTgtPointi] = ps0[i];
            }
            localTgtPointi ++;
        }
    }

    // Remote data after
    forAll(procLocalFaces, proci)
    {
        if (proci != thisProci)
        {
            const labelList& fis = procLocalFaceis[proci];
            const faceList& fs = procLocalFaces[proci];
            forAll(fis, i)
            {
                localProcFaces[localTgtFacei] = {proci, fis[i]};
                localTgtFaces[localTgtFacei] = face(fs[i] + localTgtPointi);
                localTgtFacei ++;
            }

            const pointField& ps = procLocalPoints[proci];
            const pointField& ps0 = procLocalPoints0[proci];
            forAll(ps, i)
            {
                localTgtPoints[localTgtPointi] = ps[i];
                if (patch.has0())
                {
                    localTgtPoints0[localTgtPointi] = ps0[i];
                }
                localTgtPointi ++;
            }
        }
    }

    return
        patch.has0()
      ? PrimitiveOldTimePatch<faceList, pointField>
        (
            localTgtFaces,
            localTgtPoints,
            localTgtPoints0
        )
      : PrimitiveOldTimePatch<faceList, pointField>
        (
            localTgtFaces,
            localTgtPoints
        );
}


// ************************************************************************* //
