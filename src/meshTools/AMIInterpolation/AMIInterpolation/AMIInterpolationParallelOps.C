/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "AMIInterpolation.H"
#include "mergePoints.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::label Foam::AMIInterpolation<SourcePatch, TargetPatch>::calcDistribution
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch
) const
{
    label procI = 0;

    if (Pstream::parRun())
    {
        List<label> facesPresentOnProc(Pstream::nProcs(), 0);
        if ((srcPatch.size() > 0) || (tgtPatch.size() > 0))
        {
            facesPresentOnProc[Pstream::myProcNo()] = 1;
        }
        else
        {
            facesPresentOnProc[Pstream::myProcNo()] = 0;
        }

        Pstream::gatherList(facesPresentOnProc);
        Pstream::scatterList(facesPresentOnProc);

        label nHaveFaces = sum(facesPresentOnProc);

        if (nHaveFaces > 1)
        {
            procI = -1;
            if (debug)
            {
                Info<< "AMIInterpolation::calcDistribution: "
                    << "AMI split across multiple processors" << endl;
            }
        }
        else if (nHaveFaces == 1)
        {
            procI = findIndex(facesPresentOnProc, 1);
            if (debug)
            {
                Info<< "AMIInterpolation::calcDistribution: "
                    << "AMI local to processor" << procI << endl;
            }
        }
    }


    // Either not parallel or no faces on any processor
    return procI;
}


template<class SourcePatch, class TargetPatch>
Foam::label
Foam::AMIInterpolation<SourcePatch, TargetPatch>::calcOverlappingProcs
(
    const List<treeBoundBoxList>& procBb,
    const treeBoundBox& bb,
    boolList& overlaps
) const
{
    overlaps.setSize(procBb.size());
    overlaps = false;

    label nOverlaps = 0;

    forAll(procBb, procI)
    {
        const List<treeBoundBox>& bbs = procBb[procI];

        forAll(bbs, bbI)
        {
            if (bbs[bbI].overlaps(bb))
            {
                overlaps[procI] = true;
                nOverlaps++;
                break;
            }
        }
    }
    return nOverlaps;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::distributePatches
(
    const mapDistribute& map,
    const TargetPatch& pp,
    const globalIndex& gi,
    List<faceList>& faces,
    List<pointField>& points,
    List<labelList>& faceIDs
) const
{
    PstreamBuffers pBufs(Pstream::nonBlocking);

    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& sendElems = map.subMap()[domain];

        if (domain != Pstream::myProcNo() && sendElems.size())
        {
            labelList globalElems(sendElems.size());
            forAll(sendElems, i)
            {
                globalElems[i] = gi.toGlobal(sendElems[i]);
            }

            faceList subFaces(UIndirectList<face>(pp, sendElems));
            primitivePatch subPatch
            (
                SubList<face>(subFaces, subFaces.size()),
                pp.points()
            );

            if (debug & 2)
            {
                Pout<< "distributePatches: to processor " << domain
                    << " sending faces " << subPatch.faceCentres() << endl;
            }

            UOPstream toDomain(domain, pBufs);
            toDomain
                << subPatch.localFaces() << subPatch.localPoints()
                << globalElems;
        }
    }

    // Start receiving
    pBufs.finishedSends();

    faces.setSize(Pstream::nProcs());
    points.setSize(Pstream::nProcs());
    faceIDs.setSize(Pstream::nProcs());

    {
        // Set up 'send' to myself
        const labelList& sendElems = map.subMap()[Pstream::myProcNo()];
        faceList subFaces(UIndirectList<face>(pp, sendElems));
        primitivePatch subPatch
        (
            SubList<face>(subFaces, subFaces.size()),
            pp.points()
        );

        // Receive
        if (debug & 2)
        {
            Pout<< "distributePatches: to processor " << Pstream::myProcNo()
                << " sending faces " << subPatch.faceCentres() << endl;
        }

        faces[Pstream::myProcNo()] = subPatch.localFaces();
        points[Pstream::myProcNo()] = subPatch.localPoints();

        faceIDs[Pstream::myProcNo()].setSize(sendElems.size());
        forAll(sendElems, i)
        {
            faceIDs[Pstream::myProcNo()][i] = gi.toGlobal(sendElems[i]);
        }
    }

    // Consume
    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& recvElems = map.constructMap()[domain];

        if (domain != Pstream::myProcNo() && recvElems.size())
        {
            UIPstream str(domain, pBufs);

            str >> faces[domain]
                >> points[domain]
                >> faceIDs[domain];
        }
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::
distributeAndMergePatches
(
    const mapDistribute& map,
    const TargetPatch& tgtPatch,
    const globalIndex& gi,
    faceList& tgtFaces,
    pointField& tgtPoints,
    labelList& tgtFaceIDs
) const
{
    // Exchange per-processor data
    List<faceList> allFaces;
    List<pointField> allPoints;
    List<labelList> allTgtFaceIDs;
    distributePatches(map, tgtPatch, gi, allFaces, allPoints, allTgtFaceIDs);

    // Renumber and flatten
    label nFaces = 0;
    label nPoints = 0;
    forAll(allFaces, procI)
    {
        nFaces += allFaces[procI].size();
        nPoints += allPoints[procI].size();
    }

    tgtFaces.setSize(nFaces);
    tgtPoints.setSize(nPoints);
    tgtFaceIDs.setSize(nFaces);

    nFaces = 0;
    nPoints = 0;

    // My own data first
    {
        const labelList& faceIDs = allTgtFaceIDs[Pstream::myProcNo()];
        SubList<label>(tgtFaceIDs, faceIDs.size()).assign(faceIDs);

        const faceList& fcs = allFaces[Pstream::myProcNo()];
        forAll(fcs, i)
        {
            const face& f = fcs[i];
            face& newF = tgtFaces[nFaces++];
            newF.setSize(f.size());
            forAll(f, fp)
            {
                newF[fp] = f[fp] + nPoints;
            }
        }

        const pointField& pts = allPoints[Pstream::myProcNo()];
        forAll(pts, i)
        {
            tgtPoints[nPoints++] = pts[i];
        }
    }


    // Other proc data follows
    forAll(allFaces, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            const labelList& faceIDs = allTgtFaceIDs[procI];
            SubList<label>(tgtFaceIDs, faceIDs.size(), nFaces).assign(faceIDs);

            const faceList& fcs = allFaces[procI];
            forAll(fcs, i)
            {
                const face& f = fcs[i];
                face& newF = tgtFaces[nFaces++];
                newF.setSize(f.size());
                forAll(f, fp)
                {
                    newF[fp] = f[fp] + nPoints;
                }
            }

            const pointField& pts = allPoints[procI];
            forAll(pts, i)
            {
                tgtPoints[nPoints++] = pts[i];
            }
        }
    }

    // Merge
    labelList oldToNew;
    pointField newTgtPoints;
    bool hasMerged = mergePoints
    (
        tgtPoints,
        SMALL,
        false,
        oldToNew,
        newTgtPoints
    );

    if (hasMerged)
    {
        if (debug)
        {
            Pout<< "Merged from " << tgtPoints.size()
                << " down to " << newTgtPoints.size() << " points" << endl;
        }

        tgtPoints.transfer(newTgtPoints);
        forAll(tgtFaces, i)
        {
            inplaceRenumber(oldToNew, tgtFaces[i]);
        }
    }
}


template<class SourcePatch, class TargetPatch>
Foam::autoPtr<Foam::mapDistribute>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::calcProcMap
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch
) const
{
    // Get decomposition of patch
    List<treeBoundBoxList> procBb(Pstream::nProcs());

    if (srcPatch.size())
    {
        procBb[Pstream::myProcNo()] = treeBoundBoxList
        (
            1,  // For now single bounding box per proc
            treeBoundBox
            (
                srcPatch.points(),
                srcPatch.meshPoints()
            )
        );
    }
    else
    {
        procBb[Pstream::myProcNo()] = treeBoundBoxList();
    }

    // slightly increase size of bounding boxes to allow for cases where
    // bounding boxes are perfectly alligned
    forAll(procBb[Pstream::myProcNo()], bbI)
    {
        treeBoundBox& bb = procBb[Pstream::myProcNo()][bbI];
        bb.inflate(0.01);
    }

    Pstream::gatherList(procBb);
    Pstream::scatterList(procBb);


    if (debug)
    {
        Info<< "Determining extent of srcPatch per processor:" << nl
            << "\tproc\tbb" << endl;
        forAll(procBb, procI)
        {
            Info<< '\t' << procI << '\t' << procBb[procI] << endl;
        }
    }


    // Determine which faces of tgtPatch overlaps srcPatch per proc
    const faceList& faces = tgtPatch.localFaces();
    const pointField& points = tgtPatch.localPoints();

    labelListList sendMap;

    {
        // Per processor indices into all segments to send
        List<DynamicList<label> > dynSendMap(Pstream::nProcs());

        // Work array - whether processor bb overlaps the face bounds
        boolList procBbOverlaps(Pstream::nProcs());

        forAll(faces, faceI)
        {
            if (faces[faceI].size())
            {
                treeBoundBox faceBb(points, faces[faceI]);

                // Find the processor this face overlaps
                calcOverlappingProcs(procBb, faceBb, procBbOverlaps);

                forAll(procBbOverlaps, procI)
                {
                    if (procBbOverlaps[procI])
                    {
                        dynSendMap[procI].append(faceI);
                    }
                }
            }
        }

        // Convert dynamicList to labelList
        sendMap.setSize(Pstream::nProcs());
        forAll(sendMap, procI)
        {
            sendMap[procI].transfer(dynSendMap[procI]);
        }
    }

    // Debug printing
    if (debug)
    {
        Pout<< "Of my " << faces.size() << " I need to send to:" << nl
            << "\tproc\tfaces" << endl;
        forAll(sendMap, procI)
        {
            Pout<< '\t' << procI << '\t' << sendMap[procI].size() << endl;
        }
    }


    // Send over how many faces I need to receive
    labelListList sendSizes(Pstream::nProcs());
    sendSizes[Pstream::myProcNo()].setSize(Pstream::nProcs());
    forAll(sendMap, procI)
    {
        sendSizes[Pstream::myProcNo()][procI] = sendMap[procI].size();
    }
    Pstream::gatherList(sendSizes);
    Pstream::scatterList(sendSizes);


    // Determine order of receiving
    labelListList constructMap(Pstream::nProcs());

    // My local segment first
    constructMap[Pstream::myProcNo()] = identity
    (
        sendMap[Pstream::myProcNo()].size()
    );

    label segmentI = constructMap[Pstream::myProcNo()].size();
    forAll(constructMap, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            // What I need to receive is what other processor is sending to me
            label nRecv = sendSizes[procI][Pstream::myProcNo()];
            constructMap[procI].setSize(nRecv);

            for (label i = 0; i < nRecv; i++)
            {
                constructMap[procI][i] = segmentI++;
            }
        }
    }

    autoPtr<mapDistribute> mapPtr
    (
        new mapDistribute
        (
            segmentI,       // size after construction
            sendMap.xfer(),
            constructMap.xfer()
        )
    );

    return mapPtr;
}


// ************************************************************************* //
