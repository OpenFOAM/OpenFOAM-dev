/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "patchToPatchTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::patchToPatchTools::singleProcess
(
    const label srcSize,
    const label tgtSize
)
{
    label result = 0;

    if (Pstream::parRun())
    {
        boolList procHasElements(Pstream::nProcs(), false);

        if (srcSize > 0 || tgtSize > 0)
        {
            procHasElements[Pstream::myProcNo()] = true;
        }

        Pstream::gatherList(procHasElements);
        Pstream::scatterList(procHasElements);

        const label nProcsHaveElements = count(procHasElements, true);

        if (nProcsHaveElements == 0)
        {
            result = 0;
        }

        if (nProcsHaveElements == 1)
        {
            result = findIndex(procHasElements, true);
        }

        if (nProcsHaveElements > 1)
        {
            result = -1;
        }
    }

    return result;
}


Foam::autoPtr<Foam::distributionMap>
Foam::patchToPatchTools::constructDistributionMap
(
    const labelListList& procSendIndices
)
{
    // Communicate the size of all the transfers
    labelListList transferSizes(Pstream::nProcs());
    transferSizes[Pstream::myProcNo()].setSize(Pstream::nProcs());
    forAll(procSendIndices, proci)
    {
        transferSizes[Pstream::myProcNo()][proci] =
            procSendIndices[proci].size();
    }
    Pstream::gatherList(transferSizes);
    Pstream::scatterList(transferSizes);

    // Determine order of receiving
    labelListList constructMap(Pstream::nProcs());
    label index = constructMap[Pstream::myProcNo()].size();
    forAll(constructMap, proci)
    {
        const label n = transferSizes[proci][Pstream::myProcNo()];
        constructMap[proci].setSize(n);

        for (label i = 0; i < n; i ++)
        {
            constructMap[proci][i] = index ++;
        }
    }

    // Construct and return the map
    return
        autoPtr<distributionMap>
        (
            new distributionMap
            (
                index,
                labelListList(procSendIndices),
                move(constructMap)
            )
        );
}


Foam::List<Foam::remote> Foam::patchToPatchTools::distributeAddressing
(
    const distributionMap& map
)
{
    static const label thisProci = Pstream::myProcNo();

    // Exchange per-processor data
    List<labelList> procLocalFaceis(Pstream::nProcs());
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            const labelList& sendFaceis = map.subMap()[proci];

            if (proci != thisProci && sendFaceis.size())
            {
                UOPstream(proci, pBufs)() << sendFaceis;
            }
        }

        pBufs.finishedSends();

        // Map local data
        {
            const labelList& sendFaceis = map.subMap()[thisProci];

            procLocalFaceis[thisProci] = sendFaceis;
        }

        // Receive remote data
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            if (proci != thisProci && map.constructMap()[proci].size())
            {
                UIPstream(proci, pBufs)() >> procLocalFaceis[proci];
            }
        }
    }

    // Allocate
    List<remote> localProcFaces;
    {
        label nLocalFaces = 0;
        forAll(procLocalFaceis, proci)
        {
            nLocalFaces += procLocalFaceis[proci].size();
        }
        localProcFaces.setSize(nLocalFaces);
    }

    // Construct the result
    label localFacei = 0;
    forAll(procLocalFaceis, proci)
    {
        const labelList& fis = procLocalFaceis[proci];
        forAll(fis, i)
        {
            localProcFaces[localFacei] = {proci, fis[i]};
            localFacei ++;
        }
    }

    return localProcFaces;
}


Foam::labelListList Foam::patchToPatchTools::procSendIndices
(
    const labelListList& tgtLocalSrcFaces,
    const List<remote>& localTgtProcFaces
)
{
    // Send a source face to a proc if target face on that proc intersects it
    List<labelHashSet> resultDyn(Pstream::nProcs());
    forAll(tgtLocalSrcFaces, tgtFacei)
    {
        const label tgtProci = localTgtProcFaces[tgtFacei].proci;

        forAll(tgtLocalSrcFaces[tgtFacei], i)
        {
            const label srcFacei = tgtLocalSrcFaces[tgtFacei][i];

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


Foam::labelListList Foam::patchToPatchTools::procSendIndices
(
    const List<DynamicList<label>>& tgtLocalSrcFaces,
    const List<remote>& localTgtProcFaces
)
{
    return
        procSendIndices
        (
            labelListList(tgtLocalSrcFaces),
            localTgtProcFaces
        );
}


void Foam::patchToPatchTools::trimDistributionMap
(
    const boolList& oldIsUsed,
    distributionMap& map,
    labelList& oldToNew,
    labelList& newToOld
)
{
    // Create re-indexing
    oldToNew.resize(oldIsUsed.size());
    newToOld.resize(count(oldIsUsed, true));

    oldToNew = -1;
    newToOld = -1;

    label newi = 0;
    forAll(oldIsUsed, oldi)
    {
        if (oldIsUsed[oldi])
        {
            oldToNew[oldi] = newi;
            newToOld[newi] = oldi;
            ++ newi;
        }
    }

    // Per-processor used list for construction
    List<boolList> allOldIsUsed(Pstream::nProcs());
    forAll(map.constructMap(), proci)
    {
        allOldIsUsed[proci] =
            UIndirectList<bool>(oldIsUsed, map.constructMap()[proci]);
    }

    // Communicate to form a per-processor used list for subsetting
    List<boolList> allSubOldIsUsed(Pstream::nProcs());
    Pstream::exchange<boolList, bool>(allOldIsUsed, allSubOldIsUsed);

    // Subset the sub map
    forAll(map.subMap(), proci)
    {
        label newi = 0;
        forAll(map.subMap()[proci], oldi)
        {
            if (allSubOldIsUsed[proci][oldi])
            {
                map.subMap()[proci][newi ++] =
                    map.subMap()[proci][oldi];
            }
        }
        map.subMap()[proci].resize(newi);
    }

    // Subset and renumber the construct map
    forAll(map.constructMap(), proci)
    {
        label newi = 0;
        forAll(map.constructMap()[proci], oldi)
        {
            if (allOldIsUsed[proci][oldi])
            {
                map.constructMap()[proci][newi ++] =
                    oldToNew[map.constructMap()[proci][oldi]];
            }
        }
        map.constructMap()[proci].resize(newi);
    }
}


Foam::List<Foam::List<Foam::remote>> Foam::patchToPatchTools::localToRemote
(
    const labelListList& indices,
    const List<remote>& indexToProcIndex
)
{
    List<List<remote>> result(indices.size());

    forAll(indices, thisFacei)
    {
        result[thisFacei].resize(indices[thisFacei].size());

        forAll(indices[thisFacei], i)
        {
            result[thisFacei][i] =
                isNull(indexToProcIndex)
              ? remote(Pstream::myProcNo(), indices[thisFacei][i])
              : indexToProcIndex[indices[thisFacei][i]];
        }
    }

    return result;
}


Foam::List<Foam::List<Foam::remote>> Foam::patchToPatchTools::localToRemote
(
    const List<DynamicList<label>>& indices,
    const List<remote>& indexToProcIndex
)
{
    return
        localToRemote
        (
            labelListList(indices),
            indexToProcIndex
        );
}


void Foam::patchToPatchTools::rDistributeTgtAddressing
(
    const label tgtSize,
    const distributionMap& tgtMap,
    const List<remote>& localSrcProcFaces,
    labelListList& tgtLocalSrcFaces
)
{
    // Create a map from source procFace to local source face
    HashTable<label, remote, Hash<remote>> srcProcFaceToLocal;
    forAll(localSrcProcFaces, localSrcFacei)
    {
        srcProcFaceToLocal.insert
        (
            localSrcProcFaces[localSrcFacei],
            localSrcFacei
        );
    }

    // Collect the source procFaces on the target and convert to local
    // source face addressing
    List<List<remote>> tgtSrcProcFaces = localToRemote(tgtLocalSrcFaces);

    rDistributeListList(tgtSize, tgtMap, tgtSrcProcFaces);

    tgtLocalSrcFaces.resize(tgtSize);
    forAll(tgtSrcProcFaces, tgtFacei)
    {
        tgtLocalSrcFaces[tgtFacei].resize(tgtSrcProcFaces[tgtFacei].size());

        forAll(tgtSrcProcFaces[tgtFacei], i)
        {
            tgtLocalSrcFaces[tgtFacei][i] =
                srcProcFaceToLocal[tgtSrcProcFaces[tgtFacei][i]];
        }
    }
}


void Foam::patchToPatchTools::rDistributeTgtAddressing
(
    const label tgtSize,
    const distributionMap& tgtMap,
    const List<remote>& localSrcProcFaces,
    List<DynamicList<label>>& tgtLocalSrcFaces
)
{
    labelListList tTgtLocalSrcFaces;
    transferListList(tTgtLocalSrcFaces, tgtLocalSrcFaces);
    rDistributeTgtAddressing
    (
        tgtSize,
        tgtMap,
        localSrcProcFaces,
        tTgtLocalSrcFaces
    );
    transferListList(tgtLocalSrcFaces, tTgtLocalSrcFaces);
}


// ************************************************************************* //
