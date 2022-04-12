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
#include "cpuTime.H"
#include "distributionMap.H"
#include "globalIndex.H"
#include "indexedOctree.H"
#include "treeDataPrimitivePatch.H"
#include "vtkWritePolyData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    labelList first(const UList<labelPair>& p)
    {
        labelList f(p.size());
        forAll(p, i)
        {
            f[i] = p[i].first();
        }
        return f;
    }

    labelList second(const UList<labelPair>& p)
    {
        labelList s(p.size());
        forAll(p, i)
        {
            s[i] = p[i].second();
        }
        return s;
    }

    template<class ListType>
    label findNotIndex
    (
        const ListType& l,
        typename ListType::const_reference t,
        const label start=0
    )
    {
        for (label i = start; i < l.size(); i++)
        {
            if (l[i] != t)
            {
                return i;
            }
        }

        return -1;
    }

    treeBoundBox combine(const treeBoundBox& a, const treeBoundBox& b)
    {
        return treeBoundBox(min(a.min(), b.min()), max(a.max(), b.max()));
    }
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchToPatch, 0);
    defineRunTimeSelectionTable(patchToPatch, bool);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::List<Foam::List<Foam::patchToPatch::procFace>>
Foam::patchToPatch::localFacesToProcFaces
(
    const List<DynamicList<label>>& localFaces,
    const List<procFace>& map
)
{
    List<List<procFace>> result(localFaces.size());

    forAll(localFaces, thisFacei)
    {
        result[thisFacei].resize(localFaces[thisFacei].size());

        forAll(localFaces[thisFacei], i)
        {
            result[thisFacei][i] =
                isNull(map)
              ? procFace({Pstream::myProcNo(), localFaces[thisFacei][i]})
              : map[localFaces[thisFacei][i]];
        }
    }

    return result;
}


Foam::List<Foam::DynamicList<Foam::label>>
Foam::patchToPatch::procFacesToLocalFaces
(
    const List<List<procFace>>& procFaces,
    const HashTable<label, procFace, Hash<procFace>>& map
)
{
    List<DynamicList<label>> result(procFaces.size());

    forAll(procFaces, tgtFacei)
    {
        result[tgtFacei].resize(procFaces[tgtFacei].size());

        forAll(procFaces[tgtFacei], i)
        {
            result[tgtFacei][i] =
                isNull(map)
              ? procFaces[tgtFacei][i].facei
              : map[procFaces[tgtFacei][i]];
        }
    }

    return result;
}


Foam::treeBoundBox Foam::patchToPatch::srcBox
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const label srcFacei
) const
{
    const face& srcFace = srcPatch.localFaces()[srcFacei];
    const pointField& srcPoints = srcPatch.localPoints();
    treeBoundBox box = srcBox(srcFace, srcPoints, srcPointNormals);

    if (srcPatch.has0())
    {
        const pointField& srcPoints0 = srcPatch.localPoints();
        box = combine(box, srcBox(srcFace, srcPoints0, srcPointNormals0));
    }

    return box;
}


Foam::treeBoundBox Foam::patchToPatch::srcBox
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0
) const
{
    treeBoundBox box = treeBoundBox::invertedBox;

    forAll(srcPatch, srcFacei)
    {
        box =
            combine
            (
                box,
                srcBox(srcPatch, srcPointNormals, srcPointNormals0, srcFacei)
            );
    }

    box.inflate(rootSmall);

    return box;
}


Foam::treeBoundBox Foam::patchToPatch::tgtBox
(
    const primitiveOldTimePatch& tgtPatch,
    const label tgtFacei
) const
{
    const face& tgtFace = tgtPatch.localFaces()[tgtFacei];
    const pointField& tgtPoints = tgtPatch.localPoints();
    treeBoundBox box(tgtPoints, tgtFace);

    if (tgtPatch.has0())
    {
        const pointField& tgtPoints0 = tgtPatch.localPoints();
        box = combine(box, treeBoundBox(tgtPoints0, tgtFace));
    }

    return box;
}


Foam::treeBoundBox Foam::patchToPatch::tgtBox
(
    const primitiveOldTimePatch& tgtPatch
) const
{
    treeBoundBox box = treeBoundBox::invertedBox;

    forAll(tgtPatch, tgtFacei)
    {
        box = combine(box, tgtBox(tgtPatch, tgtFacei));
    }

    box.inflate(rootSmall);

    return box;
}


bool Foam::patchToPatch::findOrIntersectFaces
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch,
    const label srcFacei,
    const label tgtFacei
)
{
    // Return if these faces have already been intersected
    auto intersected = []
    (
        const List<DynamicList<label>>& faceLocalFaces,
        const label facei,
        const label otherFacei
    )
    {
        forAll(faceLocalFaces[facei], i)
        {
            if (faceLocalFaces[facei][i] == otherFacei)
            {
                return true;
            }
        }
        return false;
    };
    if
    (
        intersected(srcLocalTgtFaces_, srcFacei, tgtFacei)
     || intersected(tgtLocalSrcFaces_, tgtFacei, srcFacei)
    )
    {
        return true;
    }

    // Try to intersect these faces
    return
        intersectFaces
        (
            srcPatch,
            srcPointNormals,
            srcPointNormals0,
            tgtPatch,
            srcFacei,
            tgtFacei
        );
}


Foam::label Foam::patchToPatch::intersectPatchQueue
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch,
    const bool isSrc,
    const DynamicList<labelPair>& queue,
    labelList& faceComplete,
    DynamicList<labelPair>& otherQueue,
    const labelList& otherFaceComplete,
    boolList& otherFaceQueued,
    boolList& otherFaceVisited
)
{
    const primitivePatch& otherPatch = isSrc ? tgtPatch : srcPatch;

    const faceList& otherLocalFaces = otherPatch.localFaces();
    //const labelListList& otherFaceFaces = otherPatch.faceFaces();
    const labelListList& otherPointFaces = otherPatch.pointFaces();

    DynamicList<label> otherQueuedFaces;
    label nFaceComplete = 0;
    forAll(queue, queuei)
    {
        const label facei = queue[queuei].first();
        const label otherFacei = queue[queuei].second();

        // Wave out from the initial target face until all intersections fail
        DynamicList<label> otherCurrentFaces(1, otherFacei);
        DynamicList<label> otherVisitedFaces(1, otherFacei);
        otherFaceVisited[otherFacei] = true;
        label otherPerimeterReached = false;
        while (otherCurrentFaces.size())
        {
            labelList otherNextFaces;
            forAll(otherCurrentFaces, otheri)
            {
                const label otherFacej = otherCurrentFaces[otheri];

                if
                (
                    findOrIntersectFaces
                    (
                        srcPatch,
                        srcPointNormals,
                        srcPointNormals0,
                        tgtPatch,
                        isSrc ? facei : otherFacej,
                        isSrc ? otherFacej : facei
                    )
                )
                {
                    /*
                    // Propagate to edge neighbours
                    forAll(otherFaceFaces[otherFacej], otherFaceFacej)
                    {
                        const label otherFacek =
                            otherFaceFaces[otherFacej][otherFaceFacej];
                        if (!otherFaceVisited[otherFacek])
                        {
                            otherFaceVisited[otherFacek] = true;
                            otherVisitedFaces.append(otherFacek);
                            otherNextFaces.append(otherFacek);
                        }
                    }
                    */

                    // Propagate to point neighbours
                    forAll(otherPatch[otherFacej], otherFacePointj)
                    {
                        const label otherPointj =
                            otherLocalFaces[otherFacej][otherFacePointj];
                        forAll(otherPointFaces[otherPointj], otherPointFacej)
                        {
                            const label otherFacek =
                                otherPointFaces[otherPointj][otherPointFacej];
                            if (!otherFaceVisited[otherFacek])
                            {
                                otherFaceVisited[otherFacek] = true;
                                otherVisitedFaces.append(otherFacek);
                                otherNextFaces.append(otherFacek);
                            }
                        }
                    }

                    // Check if this face is connected to the perimeter
                    forAll(otherPatch[otherFacej], otherFaceEdgej)
                    {
                        const label otherEdgej =
                            otherPatch.faceEdges()[otherFacej][otherFaceEdgej];

                        // !!! Two face-edges is not a sufficient condition for
                        // manifoldness. The edges also need to be numbered in
                        // opposite directions in the two faces.
                        if (otherPatch.edgeFaces()[otherEdgej].size() != 2)
                        {
                            otherPerimeterReached = true;
                        }
                    }

                    // If this face is not complete, then add it to the next
                    // iteration's queue
                    if
                    (
                        otherFaceComplete[otherFacej] == 0
                     && !otherFaceQueued[otherFacej]
                    )
                    {
                        otherFaceQueued[otherFacej] = true;
                        otherQueuedFaces.append(otherFacej);
                        otherQueue.append(labelPair(otherFacej, facei));
                    }
                }
            }

            otherCurrentFaces.transfer(otherNextFaces);
        }

        // Reset the visited array
        forAll(otherVisitedFaces, otherVisitedFacei)
        {
            otherFaceVisited[otherVisitedFaces[otherVisitedFacei]] = false;
        }

        // Set this face as complete
        nFaceComplete -= faceComplete[facei] == 2;
        faceComplete[facei] =
            max(faceComplete[facei], otherPerimeterReached ? 1 : 2);
        nFaceComplete += faceComplete[facei] == 2;
    }

    // Reset the queued array
    forAll(otherQueuedFaces, otherQueuedFacei)
    {
        otherFaceQueued[otherQueuedFaces[otherQueuedFacei]] = false;
    }

    return nFaceComplete;
}


void Foam::patchToPatch::intersectPatches
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch
)
{
    if (srcPatch.empty() || tgtPatch.empty()) return;

    if (debug)
    {
        Info<< indent << "Writing patches" << incrIndent << endl;

        const fileName srcFileName = typeName + "_srcPatch.vtk";
        Info<< indent << "Writing patch to " << srcFileName << endl;
        vtkWritePolyData::write
        (
            srcFileName,
            "source",
            false,
            srcPatch.localPoints(),
            labelList(),
            labelListList(),
            srcPatch.localFaces()
        );

        const fileName tgtFileName = typeName + "_tgtPatch.vtk";
        Info<< indent << "Writing patch to " << tgtFileName << endl;
        vtkWritePolyData::write
        (
            tgtFileName,
            "target",
            false,
            tgtPatch.localPoints(),
            labelList(),
            labelListList(),
            tgtPatch.localFaces()
        );

        Info<< decrIndent;
    }

    auto writeQueue = [&]
    (
        const label restarti,
        const label iterationi,
        const bool isSrc,
        const DynamicList<labelPair>& queue
    )
    {
        const primitivePatch& patch = isSrc ? srcPatch : tgtPatch;

        const fileName queueFileNamePart0 =
            typeName + "_restarti=" + name(restarti) + "_iterationi=";
        const fileName queueFileNamePart1 =
            word(isSrc ? "src" : "tgt") + "Queue";

        Info<< incrIndent;

        if (iterationi == 0 || iterationi % 2 == !isSrc)
        {
            const fileName queueFileName =
                queueFileNamePart0 + name(iterationi) + '_'
              + queueFileNamePart1 + ".vtk";
            Info<< indent << "Writing queue to " << queueFileName << endl;
            vtkWritePolyData::write
            (
                queueFileName,
                queueFileNamePart1,
                false,
                patch.localPoints(),
                labelList(),
                labelListList(),
                IndirectList<face>(patch.localFaces(), first(queue))
            );
        }
        else
        {
            ln
            (
                queueFileNamePart0 + name(iterationi - 1) + '_'
              + queueFileNamePart1 + ".vtk",
                queueFileNamePart0 + name(iterationi) + '_'
              + queueFileNamePart1 + ".vtk"
            );
        }

        Info<< decrIndent;
    };

    auto writeQueues = [&]
    (
        const label restarti,
        const label iterationi,
        const DynamicList<labelPair>& srcQueue,
        const DynamicList<labelPair>& tgtQueue
    )
    {
        if (debug)
        {
            Info<< indent << "Iteration #" << iterationi
                << ": processing " << srcQueue.size() << '/'
                << srcPatch.size() << " source and "
                << tgtQueue.size() << '/' << tgtPatch.size()
                << " target faces " << endl;
        }

        if (debug > 1)
        {
            writeQueue(restarti, iterationi, true, srcQueue);
            writeQueue(restarti, iterationi, false, tgtQueue);
        }
    };

    // Build a search tree for the target patch
    typedef treeDataPrimitivePatch<primitivePatch> treeType;
    const treeBoundBox tgtTreeBox =
        treeBoundBox(tgtPatch.localPoints()).extend(1e-4);
    indexedOctree<treeType> tgtTree
    (
        treeType
        (
            false,
            tgtPatch,
            indexedOctree<treeType>::perturbTol()
        ),
        tgtTreeBox,
        8,
        10,
        3
    );

    // Set up complete arrays and loop until they are full. Note that the
    // *FaceComplete lists can take three values; 0 is incomplete, 1 is
    // complete for this iteration only, 2 is fully complete.
    if (debug)
    {
        Info<< indent << "Calculating coupling" << incrIndent << endl;
    }

    label nSrcFaceComplete = 0, nTgtFaceComplete = 0;
    labelList srcFaceComplete(srcPatch.size(), 0);
    labelList tgtFaceComplete(tgtPatch.size(), 0);
    boolList srcFaceQueued(srcPatch.size(), false);
    boolList tgtFaceQueued(tgtPatch.size(), false);
    boolList srcFaceVisited(srcPatch.size(), false);
    boolList tgtFaceVisited(tgtPatch.size(), false);
    label srcFacei = 0;
    label restarti = 0;
    while (srcFacei < srcPatch.size() && srcFacei != -1)
    {
        // Consider this face only once
        srcFaceComplete[srcFacei] = 2;
        nSrcFaceComplete ++;

        // Find target faces that overlap this source face's bound box
        const labelList seedTgtFaces =
            tgtTree.findBox
            (
                srcBox(srcPatch, srcPointNormals, srcPointNormals0, srcFacei)
            );

        if (!seedTgtFaces.empty())
        {
            if (debug)
            {
                Info<< indent << "Restart #" << restarti
                    << " from at source face at "
                    << srcPatch.faceCentres()[srcFacei]
                    << incrIndent << endl;
            }

            // Initialise queues with the target faces identified
            DynamicList<labelPair> srcQueue, tgtQueue;
            forAll(seedTgtFaces, seedTgtFacei)
            {
                const label tgtFacei = seedTgtFaces[seedTgtFacei];
                srcQueue.append(labelPair(srcFacei, tgtFacei));
                tgtQueue.append(labelPair(tgtFacei, srcFacei));
            }

            if (debug)
            {
                writeQueues(restarti, 0, srcQueue, tgtQueue);
            }

            // Do intersections until queues are empty
            label iterationi = 0;
            while (true)
            {
                tgtQueue.clear();

                nSrcFaceComplete +=
                    intersectPatchQueue
                    (
                        srcPatch,
                        srcPointNormals,
                        srcPointNormals0,
                        tgtPatch,
                        true,
                        srcQueue,
                        srcFaceComplete,
                        tgtQueue,
                        tgtFaceComplete,
                        tgtFaceQueued,
                        tgtFaceVisited
                    );

                if (debug)
                {
                    writeQueues(restarti, 2*iterationi + 1, srcQueue, tgtQueue);
                }

                if (!tgtQueue.size()) break;

                srcQueue.clear();

                nTgtFaceComplete +=
                    intersectPatchQueue
                    (
                        srcPatch,
                        srcPointNormals,
                        srcPointNormals0,
                        tgtPatch,
                        false,
                        tgtQueue,
                        tgtFaceComplete,
                        srcQueue,
                        srcFaceComplete,
                        srcFaceQueued,
                        srcFaceVisited
                    );

                if (debug)
                {
                    writeQueues(restarti, 2*iterationi + 2, srcQueue, tgtQueue);
                }

                if (!srcQueue.size()) break;

                ++ iterationi;
            }

            if (debug)
            {
                Info<< indent << "Completed " << nSrcFaceComplete << '/'
                    << srcPatch.size() << " source faces " << decrIndent
                    << endl;
            }
        }

        // Find the next incomplete face
        srcFacei = findNotIndex(srcFaceComplete, 2, srcFacei);

        ++ restarti;
    }

    if (debug)
    {
        Info<< decrIndent;
    }
}


void Foam::patchToPatch::initialise
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch
)
{
    srcLocalTgtFaces_.resize(srcPatch.size());
    forAll(srcLocalTgtFaces_, i)
    {
        srcLocalTgtFaces_[i].clear();
    }

    tgtLocalSrcFaces_.resize(tgtPatch.size());
    forAll(tgtLocalSrcFaces_, i)
    {
        tgtLocalSrcFaces_[i].clear();
    }
}


Foam::tmpNrc<Foam::PrimitiveOldTimePatch<Foam::faceList, Foam::pointField>>
Foam::patchToPatch::distributeTgt
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch,
    distributionMap& tgtMap
)
{
    tgtMap =
        patchDistributionMap
        (
            tgtPatchSendFaces
            (
                srcPatch,
                srcPointNormals,
                srcPointNormals0,
                tgtPatch
            )
        );

    if (localTgtProcFacesPtr_.empty())
    {
        localTgtProcFacesPtr_.set(new List<procFace>());
    }

    return
        tmpNrc<PrimitiveOldTimePatch<faceList, pointField>>
        (
            new PrimitiveOldTimePatch<faceList, pointField>
            (
                distributePatch(tgtMap, tgtPatch, localTgtProcFacesPtr_())
            )
        );
}


Foam::tmpNrc<Foam::PrimitiveOldTimePatch<Foam::faceList, Foam::pointField>>
Foam::patchToPatch::distributeSrc
(
    const primitiveOldTimePatch& srcPatch,
    distributionMap& srcMap
)
{
    srcMap = patchDistributionMap(srcPatchSendFaces());

    if (localSrcProcFacesPtr_.empty())
    {
        localSrcProcFacesPtr_.set(new List<procFace>());
    }

    return
        tmpNrc<PrimitiveOldTimePatch<faceList, pointField>>
        (
            new PrimitiveOldTimePatch<faceList, pointField>
            (
                distributePatch(srcMap, srcPatch, localSrcProcFacesPtr_())
            )
        );
}


void Foam::patchToPatch::rDistributeTgt
(
    const primitiveOldTimePatch& tgtPatch,
    const distributionMap& tgtMap
)
{
    // Create a map from source procFace to local source face
    HashTable<label, procFace, Hash<procFace>> srcProcFaceToLocal;
    forAll(localSrcProcFacesPtr_(), localSrcFacei)
    {
        srcProcFaceToLocal.insert
        (
            localSrcProcFacesPtr_()[localSrcFacei],
            localSrcFacei
        );
    }

    // Collect the source procFaces on the target and convert to local
    // source face addressing
    List<List<procFace>> tgtSrcProcFaces =
        localFacesToProcFaces(tgtLocalSrcFaces_);

    rDistributeListList(tgtPatch.size(), tgtMap, tgtSrcProcFaces);

    tgtLocalSrcFaces_ =
        procFacesToLocalFaces(tgtSrcProcFaces, srcProcFaceToLocal);
}


Foam::label Foam::patchToPatch::finalise
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch,
    const transformer& tgtToSrc
)
{
    label nCouples = 0;

    forAll(srcLocalTgtFaces_, srcFacei)
    {
        nCouples += srcLocalTgtFaces_[srcFacei].size();
    }
    forAll(tgtLocalSrcFaces_, tgtFacei)
    {
        nCouples += tgtLocalSrcFaces_[tgtFacei].size();
    }

    reduce(nCouples, sumOp<label>());

    return nCouples;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatch::patchToPatch(const bool reverse)
:
    reverse_(reverse),
    singleProcess_(-labelMax),
    localSrcProcFacesPtr_(nullptr),
    localTgtProcFacesPtr_(nullptr),
    srcLocalTgtFaces_(),
    tgtLocalSrcFaces_()
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::patchToPatch::~patchToPatch()
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::patchToPatch> Foam::patchToPatch::New
(
    const word& patchToPatchType,
    const bool reverse
)
{
    boolConstructorTable::iterator cstrIter =
        boolConstructorTablePtr_->find(patchToPatchType);

    if (cstrIter == boolConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type "
            << patchToPatchType << endl << endl
            << "Valid " << typeName << " types are : " << endl
            << boolConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(reverse);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchToPatch::update
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch,
    const transformer& tgtToSrc
)
{
    cpuTime time;

    // Determine numbers of faces on both sides, report, and quit if either
    // side is empty
    const label srcTotalSize = returnReduce(srcPatch.size(), sumOp<label>());
    const label tgtTotalSize = returnReduce(tgtPatch.size(), sumOp<label>());
    if (srcTotalSize == 0 || tgtTotalSize == 0)
    {
        return;
    }

    // If a transformation is given then transform the target to the source
    tmpNrc<primitiveOldTimePatch> tTgtPatchPtr(tgtPatch);
    tmpNrc<pointField> tTgtPointsPtr(tgtPatch.localPoints());
    tmpNrc<pointField> tTgtPoints0Ptr
    (
        tgtPatch.has0() ? tgtPatch.localPoints0() : NullObjectRef<pointField>()
    );
    if (!isNull(tgtToSrc))
    {
        tTgtPointsPtr = new pointField(tgtPatch.localPoints());
        tgtToSrc.transformPosition
        (
            tTgtPointsPtr.ref(),
            tTgtPointsPtr.ref()
        );

        if (tgtPatch.has0())
        {
            tTgtPoints0Ptr = new pointField(tgtPatch.localPoints0());
            tgtToSrc.transformPosition
            (
                tTgtPoints0Ptr.ref(),
                tTgtPoints0Ptr.ref()
            );
        }

        tTgtPatchPtr =
            tgtPatch.has0()
          ? new primitiveOldTimePatch
            (
                SubList<face>(tgtPatch.localFaces(), tgtPatch.size()),
                tTgtPointsPtr(),
                tTgtPoints0Ptr()
            )
          : new primitiveOldTimePatch
            (
                SubList<face>(tgtPatch.localFaces(), tgtPatch.size()),
                tTgtPointsPtr()
            );
    }
    const primitiveOldTimePatch& tTgtPatch = tTgtPatchPtr();

    Info<< indent << typeName << ": Calculating couplings between "
        << srcTotalSize << " source faces and " << tgtTotalSize
        << " target faces" << incrIndent << endl;

    // Determine if patches are present on multiple processors
    calcSingleProcess(srcPatch, tTgtPatch);

    // Do intersection in serial or parallel as appropriate
    if (isSingleProcess())
    {
        // Initialise the workspace
        initialise(srcPatch, srcPointNormals, srcPointNormals0, tTgtPatch);

        // Intersect the patches
        const treeBoundBox srcPatchBox =
            srcBox(srcPatch, srcPointNormals, srcPointNormals0);
        const treeBoundBox tgtPatchBox = tgtBox(tTgtPatch);
        if (srcPatchBox.overlaps(tgtPatchBox))
        {
            intersectPatches
            (
                srcPatch,
                srcPointNormals,
                srcPointNormals0,
                tTgtPatch
            );
        }
    }
    else
    {
        // Distribute the target
        distributionMap tgtMap;
        tmpNrc<PrimitiveOldTimePatch<faceList, pointField>> localTTgtPatchPtr =
            distributeTgt
            (
                srcPatch,
                srcPointNormals,
                srcPointNormals0,
                tTgtPatch,
                tgtMap
            );

        // Massage target patch into form that can be used by the serial
        // intersection interface
        const primitiveOldTimePatch localTTgtPatch =
            tgtPatch.has0()
          ? primitiveOldTimePatch
            (
                SubList<face>(localTTgtPatchPtr(), localTTgtPatchPtr().size()),
                localTTgtPatchPtr().points(),
                localTTgtPatchPtr().points0()
            )
          : primitiveOldTimePatch
            (
                SubList<face>(localTTgtPatchPtr(), localTTgtPatchPtr().size()),
                localTTgtPatchPtr().points()
            );

        // Initialise the workspace
        initialise(srcPatch, srcPointNormals, srcPointNormals, localTTgtPatch);

        // Intersect the patches
        if (localTTgtPatch.size())
        {
            intersectPatches
            (
                srcPatch,
                srcPointNormals,
                srcPointNormals,
                localTTgtPatch
            );
        }

        // Distribute the source
        distributionMap srcMap;
        tmpNrc<PrimitiveOldTimePatch<faceList, pointField>> localSrcPatchPtr =
            distributeSrc(srcPatch, srcMap);

        // Reverse distribute coupling data back to the target
        rDistributeTgt(tgtPatch, tgtMap);
    }

    // Finalise the intersection
    const label nCouples =
        finalise
        (
            srcPatch,
            srcPointNormals,
            srcPointNormals,
            tgtPatch,
            tgtToSrc
        );

    if (nCouples != 0)
    {
        Info<< indent << nCouples << " couplings calculated in "
            << time.cpuTimeIncrement() << 's' << endl;
    }
    else
    {
        Info<< indent << "No couplings found" << endl;
    }

    Info<< decrIndent;
}


void Foam::patchToPatch::update
(
    const primitivePatch& srcPatch,
    const vectorField& srcPointNormals,
    const primitivePatch& tgtPatch,
    const transformer& tgtToSrc
)
{
    update
    (
        primitiveOldTimePatch(srcPatch),
        srcPointNormals,
        NullObjectRef<vectorField>(),
        primitiveOldTimePatch(tgtPatch),
        tgtToSrc
    );
}


// ************************************************************************* //
