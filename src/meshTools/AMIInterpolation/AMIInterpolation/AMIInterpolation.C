/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
#include "AMIMethod.H"
#include "meshTools.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::word
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolationMethodToWord
(
    const interpolationMethod& im
)
{
    word method = "unknown-interpolationMethod";

    switch (im)
    {
        case imDirect:
        {
            method = "directAMI";
            break;
        }
        case imMapNearest:
        {
            method = "mapNearestAMI";
            break;
        }
        case imFaceAreaWeight:
        {
            method = "faceAreaWeightAMI";
            break;
        }
        case imPartialFaceAreaWeight:
        {
            method = "partialFaceAreaWeightAMI";
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "const Foam::word"
                "Foam::AMIInterpolation<SourcePatch, TargetPatch>::"
                "interpolationMethodToWord"
                "("
                    "const interpolationMethod&"
                ")"
            )
                << "Unhandled interpolationMethod enumeration " << method
                << abort(FatalError);
        }
    }

    return method;
}


template<class SourcePatch, class TargetPatch>
typename Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolationMethod
Foam::AMIInterpolation<SourcePatch, TargetPatch>::wordTointerpolationMethod
(
    const word& im
)
{
    interpolationMethod method = imDirect;

    wordList methods
    (
        IStringStream
        (
            "("
                "directAMI "
                "mapNearestAMI "
                "faceAreaWeightAMI "
                "partialFaceAreaWeightAMI"
            ")"
        )()
    );

    if (im == "directAMI")
    {
        method = imDirect;
    }
    else if (im == "mapNearestAMI")
    {
        method = imMapNearest;
    }
    else if (im == "faceAreaWeightAMI")
    {
        method = imFaceAreaWeight;
    }
    else if (im == "partialFaceAreaWeightAMI")
    {
        method = imPartialFaceAreaWeight;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::AMIInterpolation<SourcePatch, TargetPatch>::"
            "interpolationMethod"
            "Foam::AMIInterpolation<SourcePatch, TargetPatch>::"
            "wordTointerpolationMethod"
            "("
                "const word&"
            ")"
        )
            << "Invalid interpolationMethod " << im
            << ".  Valid methods are:" << methods
            << exit(FatalError);
    }

    return method;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::projectPointsToSurface
(
    const searchableSurface& surf,
    pointField& pts
) const
{
    if (debug)
    {
        Info<< "AMI: projecting points to surface" << endl;
    }

    List<pointIndexHit> nearInfo;

    surf.findNearest(pts, scalarField(pts.size(), GREAT), nearInfo);

    label nMiss = 0;
    forAll(nearInfo, i)
    {
        const pointIndexHit& pi = nearInfo[i];

        if (pi.hit())
        {
            pts[i] = pi.hitPoint();
        }
        else
        {
            pts[i] = pts[i];
            nMiss++;
        }
    }

    if (nMiss > 0)
    {
        FatalErrorIn
        (
            "void Foam::AMIInterpolation<SourcePatch, TargetPatch>::"
            "projectPointsToSurface"
            "("
                "const searchableSurface&, "
                "pointField&"
            ") const"
        )
            << "Error projecting points to surface: "
            << nMiss << " faces could not be determined"
            << abort(FatalError);
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::normaliseWeights
(
    const scalarField& patchAreas,
    const word& patchName,
    const labelListList& addr,
    scalarListList& wght,
    scalarField& wghtSum,
    const bool conformal,
    const bool output,
    const scalar lowWeightTol
)
{
    // Normalise the weights
    wghtSum.setSize(wght.size(), 0.0);
    label nLowWeight = 0;

    forAll(wght, faceI)
    {
        scalarList& w = wght[faceI];

        if (w.size())
        {
            scalar denom = patchAreas[faceI];

            scalar s = sum(w);
            scalar t = s/denom;

            if (conformal)
            {
                denom = s;
            }

            forAll(w, i)
            {
                w[i] /= denom;
            }

            wghtSum[faceI] = t;

            if (t < lowWeightTol)
            {
                nLowWeight++;
            }
        }
        else
        {
            wghtSum[faceI] = 0;
        }
    }


    if (output)
    {
        const label nFace = returnReduce(wght.size(), sumOp<label>());

        if (nFace)
        {
            Info<< indent
                << "AMI: Patch " << patchName
                << " sum(weights) min/max/average = "
                << gMin(wghtSum) << ", "
                << gMax(wghtSum) << ", "
                << gAverage(wghtSum) << endl;

            const label nLow = returnReduce(nLowWeight, sumOp<label>());

            if (nLow)
            {
                Info<< indent
                    << "AMI: Patch " << patchName
                    << " identified " << nLow
                    << " faces with weights less than " << lowWeightTol
                    << endl;
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::agglomerate
(
    const autoPtr<mapDistribute>& targetMapPtr,
    const scalarField& fineSrcMagSf,
    const labelListList& fineSrcAddress,
    const scalarListList& fineSrcWeights,

    const labelList& sourceRestrictAddressing,
    const labelList& targetRestrictAddressing,

    scalarField& srcMagSf,
    labelListList& srcAddress,
    scalarListList& srcWeights,
    scalarField& srcWeightsSum,
    autoPtr<mapDistribute>& tgtMap
)
{
    label sourceCoarseSize =
    (
        sourceRestrictAddressing.size()
      ? max(sourceRestrictAddressing)+1
      : 0
    );

    label targetCoarseSize =
    (
        targetRestrictAddressing.size()
      ? max(targetRestrictAddressing)+1
      : 0
    );

    // Agglomerate face areas
    {
        srcMagSf.setSize(sourceRestrictAddressing.size(), 0.0);
        forAll(sourceRestrictAddressing, faceI)
        {
            label coarseFaceI = sourceRestrictAddressing[faceI];
            srcMagSf[coarseFaceI] += fineSrcMagSf[faceI];
        }
    }


    // Agglomerate weights and indices
    if (targetMapPtr.valid())
    {
        const mapDistribute& map = targetMapPtr();

        // Get all restriction addressing.
        labelList allRestrict(targetRestrictAddressing);
        map.distribute(allRestrict);

        // So now we have agglomeration of the target side in
        // allRestrict:
        //  0..size-1 : local agglomeration (= targetRestrictAddressing)
        //  size..    : agglomeration data from other processors

        labelListList tgtSubMap(Pstream::nProcs());

        // Local subMap is just identity
        {
            tgtSubMap[Pstream::myProcNo()] = identity(targetCoarseSize);
        }

        forAll(map.subMap(), procI)
        {
            if (procI != Pstream::myProcNo())
            {
                // Combine entries that point to the same coarse element. All
                // the elements refer to local data so index into
                // targetRestrictAddressing or allRestrict (since the same
                // for local data).
                const labelList& elems = map.subMap()[procI];
                labelList& newSubMap = tgtSubMap[procI];
                newSubMap.setSize(elems.size());

                labelList oldToNew(targetCoarseSize, -1);
                label newI = 0;

                forAll(elems, i)
                {
                    label fineElem = elems[i];
                    label coarseElem = allRestrict[fineElem];
                    if (oldToNew[coarseElem] == -1)
                    {
                        oldToNew[coarseElem] = newI;
                        newSubMap[newI] = coarseElem;
                        newI++;
                    }
                }
                newSubMap.setSize(newI);
            }
        }

        // Reconstruct constructMap by combining entries. Note that order
        // of handing out indices should be the same as loop above to compact
        // the sending map

        labelListList tgtConstructMap(Pstream::nProcs());
        labelList tgtCompactMap;

        // Local constructMap is just identity
        {
            tgtConstructMap[Pstream::myProcNo()] =
                identity(targetCoarseSize);

            tgtCompactMap = targetRestrictAddressing;
        }
        tgtCompactMap.setSize(map.constructSize());
        label compactI = targetCoarseSize;

        // Compact data from other processors
        forAll(map.constructMap(), procI)
        {
            if (procI != Pstream::myProcNo())
            {
                // Combine entries that point to the same coarse element. All
                // elements now are remote data so we cannot use any local
                // data here - use allRestrict instead.
                const labelList& elems = map.constructMap()[procI];

                labelList& newConstructMap = tgtConstructMap[procI];
                newConstructMap.setSize(elems.size());

                if (elems.size())
                {
                    // Get the maximum target coarse size for this set of
                    // received data.
                    label remoteTargetCoarseSize = labelMin;
                    forAll(elems, i)
                    {
                        remoteTargetCoarseSize = max
                        (
                            remoteTargetCoarseSize,
                            allRestrict[elems[i]]
                        );
                    }
                    remoteTargetCoarseSize += 1;

                    // Combine locally data coming from procI
                    labelList oldToNew(remoteTargetCoarseSize, -1);
                    label newI = 0;

                    forAll(elems, i)
                    {
                        label fineElem = elems[i];
                        // fineElem now points to section from procI
                        label coarseElem = allRestrict[fineElem];
                        if (oldToNew[coarseElem] == -1)
                        {
                            oldToNew[coarseElem] = newI;
                            tgtCompactMap[fineElem] = compactI;
                            newConstructMap[newI] = compactI++;
                            newI++;
                        }
                        else
                        {
                            // Get compact index
                            label compactI = oldToNew[coarseElem];
                            tgtCompactMap[fineElem] = newConstructMap[compactI];
                        }
                    }
                    newConstructMap.setSize(newI);
                }
            }
        }


        srcAddress.setSize(sourceCoarseSize);
        srcWeights.setSize(sourceCoarseSize);

        forAll(fineSrcAddress, faceI)
        {
            // All the elements contributing to faceI. Are slots in
            // mapDistribute'd data.
            const labelList& elems = fineSrcAddress[faceI];
            const scalarList& weights = fineSrcWeights[faceI];
            const scalar fineArea = fineSrcMagSf[faceI];

            label coarseFaceI = sourceRestrictAddressing[faceI];

            labelList& newElems = srcAddress[coarseFaceI];
            scalarList& newWeights = srcWeights[coarseFaceI];

            forAll(elems, i)
            {
                label elemI = elems[i];
                label coarseElemI = tgtCompactMap[elemI];

                label index = findIndex(newElems, coarseElemI);
                if (index == -1)
                {
                    newElems.append(coarseElemI);
                    newWeights.append(fineArea*weights[i]);
                }
                else
                {
                    newWeights[index] += fineArea*weights[i];
                }
            }
        }

        tgtMap.reset
        (
            new mapDistribute
            (
                compactI,
                tgtSubMap.xfer(),
                tgtConstructMap.xfer()
            )
        );
    }
    else
    {
        srcAddress.setSize(sourceCoarseSize);
        srcWeights.setSize(sourceCoarseSize);

        forAll(fineSrcAddress, faceI)
        {
            // All the elements contributing to faceI. Are slots in
            // mapDistribute'd data.
            const labelList& elems = fineSrcAddress[faceI];
            const scalarList& weights = fineSrcWeights[faceI];
            const scalar fineArea = fineSrcMagSf[faceI];

            label coarseFaceI = sourceRestrictAddressing[faceI];

            labelList& newElems = srcAddress[coarseFaceI];
            scalarList& newWeights = srcWeights[coarseFaceI];

            forAll(elems, i)
            {
                label elemI = elems[i];
                label coarseElemI = targetRestrictAddressing[elemI];

                label index = findIndex(newElems, coarseElemI);
                if (index == -1)
                {
                    newElems.append(coarseElemI);
                    newWeights.append(fineArea*weights[i]);
                }
                else
                {
                    newWeights[index] += fineArea*weights[i];
                }
            }
        }
    }

    // weights normalisation
    normaliseWeights
    (
        srcMagSf,
        "source",
        srcAddress,
        srcWeights,
        srcWeightsSum,
        true,
        false,
        -1
    );
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::constructFromSurface
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr
)
{
    if (surfPtr.valid())
    {
        // create new patches for source and target
        pointField srcPoints = srcPatch.points();
        SourcePatch srcPatch0
        (
            SubList<face>
            (
                srcPatch,
                srcPatch.size(),
                0
            ),
            srcPoints
        );

        if (debug)
        {
            OFstream os("amiSrcPoints.obj");
            forAll(srcPoints, i)
            {
                meshTools::writeOBJ(os, srcPoints[i]);
            }
        }

        pointField tgtPoints = tgtPatch.points();
        TargetPatch tgtPatch0
        (
            SubList<face>
            (
                tgtPatch,
                tgtPatch.size(),
                0
            ),
            tgtPoints
        );

        if (debug)
        {
            OFstream os("amiTgtPoints.obj");
            forAll(tgtPoints, i)
            {
                meshTools::writeOBJ(os, tgtPoints[i]);
            }
        }


        // map source and target patches onto projection surface
        projectPointsToSurface(surfPtr(), srcPoints);
        projectPointsToSurface(surfPtr(), tgtPoints);


        // calculate AMI interpolation
        update(srcPatch0, tgtPatch0);
    }
    else
    {
        update(srcPatch, tgtPatch);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const interpolationMethod& method,
    const scalar lowWeightCorrection,
    const bool reverseTarget
)
:
    methodName_(interpolationMethodToWord(method)),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    singlePatchProc_(-999),
    lowWeightCorrection_(lowWeightCorrection),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    triMode_(triMode),
    srcMapPtr_(NULL),
    tgtMapPtr_(NULL)
{
    update(srcPatch, tgtPatch);
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const word& methodName,
    const scalar lowWeightCorrection,
    const bool reverseTarget
)
:
    methodName_(methodName),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    singlePatchProc_(-999),
    lowWeightCorrection_(lowWeightCorrection),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    triMode_(triMode),
    srcMapPtr_(NULL),
    tgtMapPtr_(NULL)
{
    update(srcPatch, tgtPatch);
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const interpolationMethod& method,
    const scalar lowWeightCorrection,
    const bool reverseTarget
)
:
    methodName_(interpolationMethodToWord(method)),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    singlePatchProc_(-999),
    lowWeightCorrection_(lowWeightCorrection),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    triMode_(triMode),
    srcMapPtr_(NULL),
    tgtMapPtr_(NULL)
{
    constructFromSurface(srcPatch, tgtPatch, surfPtr);
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const word& methodName,
    const scalar lowWeightCorrection,
    const bool reverseTarget
)
:
    methodName_(methodName),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    singlePatchProc_(-999),
    lowWeightCorrection_(lowWeightCorrection),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    triMode_(triMode),
    srcMapPtr_(NULL),
    tgtMapPtr_(NULL)
{
    constructFromSurface(srcPatch, tgtPatch, surfPtr);
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const AMIInterpolation<SourcePatch, TargetPatch>& fineAMI,
    const labelList& sourceRestrictAddressing,
    const labelList& targetRestrictAddressing
)
:
    methodName_(fineAMI.methodName_),
    reverseTarget_(fineAMI.reverseTarget_),
    requireMatch_(fineAMI.requireMatch_),
    singlePatchProc_(fineAMI.singlePatchProc_),
    lowWeightCorrection_(-1.0),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    triMode_(fineAMI.triMode_),
    srcMapPtr_(NULL),
    tgtMapPtr_(NULL)
{
    label sourceCoarseSize =
    (
        sourceRestrictAddressing.size()
      ? max(sourceRestrictAddressing)+1
      : 0
    );

    label neighbourCoarseSize =
    (
        targetRestrictAddressing.size()
      ? max(targetRestrictAddressing)+1
      : 0
    );

    if (debug & 2)
    {
        Pout<< "AMI: Creating addressing and weights as agglomeration of AMI :"
            << " source:" << fineAMI.srcAddress().size()
            << " target:" << fineAMI.tgtAddress().size()
            << " coarse source size:" << sourceCoarseSize
            << " neighbour source size:" << neighbourCoarseSize
            << endl;
    }

    if
    (
        fineAMI.srcAddress().size() != sourceRestrictAddressing.size()
     || fineAMI.tgtAddress().size() != targetRestrictAddressing.size()
    )
    {
        FatalErrorIn
        (
            "AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation"
            "("
                "const AMIInterpolation<SourcePatch, TargetPatch>&, "
                "const labelList&, "
                "const labelList&"
            ")"
        )   << "Size mismatch." << nl
            << "Source patch size:" << fineAMI.srcAddress().size() << nl
            << "Source agglomeration size:"
            << sourceRestrictAddressing.size() << nl
            << "Target patch size:" << fineAMI.tgtAddress().size() << nl
            << "Target agglomeration size:"
            << targetRestrictAddressing.size()
            << exit(FatalError);
    }


    // Agglomerate addresses and weights

    agglomerate
    (
        fineAMI.tgtMapPtr_,
        fineAMI.srcMagSf(),
        fineAMI.srcAddress(),
        fineAMI.srcWeights(),

        sourceRestrictAddressing,
        targetRestrictAddressing,

        srcMagSf_,
        srcAddress_,
        srcWeights_,
        srcWeightsSum_,
        tgtMapPtr_
    );

    //if (tgtMapPtr_.valid())
    //{
    //    Pout<< "tgtMap:" << endl;
    //    string oldPrefix = Pout.prefix();
    //    Pout.prefix() = oldPrefix + "  ";
    //    tgtMapPtr_().printLayout(Pout);
    //    Pout.prefix() = oldPrefix;
    //}

    agglomerate
    (
        fineAMI.srcMapPtr_,
        fineAMI.tgtMagSf(),
        fineAMI.tgtAddress(),
        fineAMI.tgtWeights(),

        targetRestrictAddressing,
        sourceRestrictAddressing,

        tgtMagSf_,
        tgtAddress_,
        tgtWeights_,
        tgtWeightsSum_,
        srcMapPtr_
    );

    //if (srcMapPtr_.valid())
    //{
    //    Pout<< "srcMap:" << endl;
    //    string oldPrefix = Pout.prefix();
    //    Pout.prefix() = oldPrefix + "  ";
    //    srcMapPtr_().printLayout(Pout);
    //    Pout.prefix() = oldPrefix;
    //}
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::~AMIInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::update
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch
)
{
    label srcTotalSize = returnReduce(srcPatch.size(), sumOp<label>());
    label tgtTotalSize = returnReduce(tgtPatch.size(), sumOp<label>());

    if (srcTotalSize == 0)
    {
        if (debug)
        {
            Info<< "AMI: no source faces present - no addressing constructed"
                << endl;
        }

        return;
    }

    Info<< indent
        << "AMI: Creating addressing and weights between "
        << srcTotalSize << " source faces and "
        << tgtTotalSize << " target faces"
        << endl;

    // Calculate face areas
    srcMagSf_.setSize(srcPatch.size());
    forAll(srcMagSf_, faceI)
    {
        srcMagSf_[faceI] = srcPatch[faceI].mag(srcPatch.points());
    }
    tgtMagSf_.setSize(tgtPatch.size());
    forAll(tgtMagSf_, faceI)
    {
        tgtMagSf_[faceI] = tgtPatch[faceI].mag(tgtPatch.points());
    }

    // Calculate if patches present on multiple processors
    singlePatchProc_ = calcDistribution(srcPatch, tgtPatch);

    if (singlePatchProc_ == -1)
    {
        // convert local addressing to global addressing
        globalIndex globalSrcFaces(srcPatch.size());
        globalIndex globalTgtFaces(tgtPatch.size());

        // Create processor map of overlapping faces. This map gets
        // (possibly remote) faces from the tgtPatch such that they (together)
        // cover all of the srcPatch
        autoPtr<mapDistribute> mapPtr = calcProcMap(srcPatch, tgtPatch);
        const mapDistribute& map = mapPtr();

        // create new target patch that fully encompasses source patch

        // faces and points
        faceList newTgtFaces;
        pointField newTgtPoints;
        // original faces from tgtPatch (in globalIndexing since might be
        // remote)
        labelList tgtFaceIDs;
        distributeAndMergePatches
        (
            map,
            tgtPatch,
            globalTgtFaces,
            newTgtFaces,
            newTgtPoints,
            tgtFaceIDs
        );

        TargetPatch
            newTgtPatch
            (
                SubList<face>
                (
                    newTgtFaces,
                    newTgtFaces.size()
                ),
                newTgtPoints
            );

        // calculate AMI interpolation
        autoPtr<AMIMethod<SourcePatch, TargetPatch> > AMIPtr
        (
            AMIMethod<SourcePatch, TargetPatch>::New
            (
                methodName_,
                srcPatch,
                newTgtPatch,
                srcMagSf_,
                tgtMagSf_,
                triMode_,
                reverseTarget_,
                requireMatch_
            )
        );

        AMIPtr->calculate
        (
            srcAddress_,
            srcWeights_,
            tgtAddress_,
            tgtWeights_
        );

        // Now
        // ~~~
        //  srcAddress_ :   per srcPatch face a list of the newTgtPatch (not
        //                  tgtPatch) faces it overlaps
        //  tgtAddress_ :   per newTgtPatch (not tgtPatch) face a list of the
        //                  srcPatch faces it overlaps


        // Rework newTgtPatch indices into globalIndices of tgtPatch
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        forAll(srcAddress_, i)
        {
            labelList& addressing = srcAddress_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = tgtFaceIDs[addressing[addrI]];
            }
        }

        forAll(tgtAddress_, i)
        {
            labelList& addressing = tgtAddress_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = globalSrcFaces.toGlobal(addressing[addrI]);
            }
        }

        // send data back to originating procs. Note that contributions
        // from different processors get added (ListAppendEqOp)

        mapDistribute::distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(),
            tgtPatch.size(),
            map.constructMap(),
            map.subMap(),
            tgtAddress_,
            ListAppendEqOp<label>(),
            labelList()
        );

        mapDistribute::distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(),
            tgtPatch.size(),
            map.constructMap(),
            map.subMap(),
            tgtWeights_,
            ListAppendEqOp<scalar>(),
            scalarList()
        );

        // weights normalisation
        normaliseWeights
        (
            srcMagSf_,
            "source",
            srcAddress_,
            srcWeights_,
            srcWeightsSum_,
            AMIPtr->conformal(),
            true,
            lowWeightCorrection_
        );
        normaliseWeights
        (
            tgtMagSf_,
            "target",
            tgtAddress_,
            tgtWeights_,
            tgtWeightsSum_,
            AMIPtr->conformal(),
            true,
            lowWeightCorrection_
        );

        // cache maps and reset addresses
        List<Map<label> > cMap;
        srcMapPtr_.reset(new mapDistribute(globalSrcFaces, tgtAddress_, cMap));
        tgtMapPtr_.reset(new mapDistribute(globalTgtFaces, srcAddress_, cMap));

        if (debug)
        {
            writeFaceConnectivity(srcPatch, newTgtPatch, srcAddress_);
        }
    }
    else
    {
        // calculate AMI interpolation
        autoPtr<AMIMethod<SourcePatch, TargetPatch> > AMIPtr
        (
            AMIMethod<SourcePatch, TargetPatch>::New
            (
                methodName_,
                srcPatch,
                tgtPatch,
                srcMagSf_,
                tgtMagSf_,
                triMode_,
                reverseTarget_,
                requireMatch_
            )
        );

        AMIPtr->calculate
        (
            srcAddress_,
            srcWeights_,
            tgtAddress_,
            tgtWeights_
        );

        normaliseWeights
        (
            srcMagSf_,
            "source",
            srcAddress_,
            srcWeights_,
            srcWeightsSum_,
            AMIPtr->conformal(),
            true,
            lowWeightCorrection_
        );
        normaliseWeights
        (
            tgtMagSf_,
            "target",
            tgtAddress_,
            tgtWeights_,
            tgtWeightsSum_,
            AMIPtr->conformal(),
            true,
            lowWeightCorrection_
        );
    }

    if (debug)
    {
        Info<< "AMIInterpolation : Constructed addressing and weights" << nl
            << "    triMode        :"
            << faceAreaIntersect::triangulationModeNames_[triMode_] << nl
            << "    singlePatchProc:" << singlePatchProc_ << nl
            << "    srcMagSf       :" << gSum(srcMagSf_) << nl
            << "    tgtMagSf       :" << gSum(tgtMagSf_) << nl
            << endl;
    }
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    if (fld.size() != srcAddress_.size())
    {
        FatalErrorIn
        (
            "AMIInterpolation::interpolateToTarget"
            "("
                "const UList<Type>&, "
                "const CombineOp&, "
                "List<Type>&, "
                "const UList<Type>&"
            ") const"
        )   << "Supplied field size is not equal to source patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }

    if (lowWeightCorrection_ > 0)
    {
        if (defaultValues.size() != tgtAddress_.size())
        {
            FatalErrorIn
            (
                "AMIInterpolation::interpolateToTarget"
                "("
                    "const UList<Type>&, "
                    "const CombineOp&, "
                    "List<Type>&, "
                    "const UList<Type>&"
                ") const"
            )   << "Employing default values when sum of weights falls below "
                << lowWeightCorrection_
                << " but supplied default field size is not equal to target "
                << "patch size" << nl
                << "    default values = " << defaultValues.size() << nl
                << "    target patch   = " << tgtAddress_.size() << nl
                << abort(FatalError);
        }
    }

    result.setSize(tgtAddress_.size());

    if (singlePatchProc_ == -1)
    {
        const mapDistribute& map = srcMapPtr_();

        List<Type> work(fld);
        map.distribute(work);

        forAll(result, faceI)
        {
            if (tgtWeightsSum_[faceI] < lowWeightCorrection_)
            {
                result[faceI] = defaultValues[faceI];
            }
            else
            {
                const labelList& faces = tgtAddress_[faceI];
                const scalarList& weights = tgtWeights_[faceI];

                forAll(faces, i)
                {
                    cop(result[faceI], faceI, work[faces[i]], weights[i]);
                }
            }
        }
    }
    else
    {
        forAll(result, faceI)
        {
            if (tgtWeightsSum_[faceI] < lowWeightCorrection_)
            {
                result[faceI] = defaultValues[faceI];
            }
            else
            {
                const labelList& faces = tgtAddress_[faceI];
                const scalarList& weights = tgtWeights_[faceI];

                forAll(faces, i)
                {
                    cop(result[faceI], faceI, fld[faces[i]], weights[i]);
                }
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    if (fld.size() != tgtAddress_.size())
    {
        FatalErrorIn
        (
            "AMIInterpolation::interpolateToSource"
            "("
                "const UList<Type>&, "
                "const CombineOp&, "
                "List<Type>&, "
                "const UList<Type>&"
            ") const"
        )   << "Supplied field size is not equal to target patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }

    if (lowWeightCorrection_ > 0)
    {
        if (defaultValues.size() != srcAddress_.size())
        {
            FatalErrorIn
            (
                "AMIInterpolation::interpolateToSource"
                "("
                    "const UList<Type>&, "
                    "const CombineOp&, "
                    "List<Type>&, "
                    "const UList<Type>&"
                ") const"
            )   << "Employing default values when sum of weights falls below "
                << lowWeightCorrection_
                << " but supplied default field size is not equal to target "
                << "patch size" << nl
                << "    default values = " << defaultValues.size() << nl
                << "    source patch   = " << srcAddress_.size() << nl
                << abort(FatalError);
        }
    }

    result.setSize(srcAddress_.size());

    if (singlePatchProc_ == -1)
    {
        const mapDistribute& map = tgtMapPtr_();

        List<Type> work(fld);
        map.distribute(work);

        forAll(result, faceI)
        {
            if (srcWeightsSum_[faceI] < lowWeightCorrection_)
            {
                result[faceI] = defaultValues[faceI];
            }
            else
            {
                const labelList& faces = srcAddress_[faceI];
                const scalarList& weights = srcWeights_[faceI];

                forAll(faces, i)
                {
                    cop(result[faceI], faceI, work[faces[i]], weights[i]);
                }
            }
        }
    }
    else
    {
        forAll(result, faceI)
        {
            if (srcWeightsSum_[faceI] < lowWeightCorrection_)
            {
                result[faceI] = defaultValues[faceI];
            }
            else
            {
                const labelList& faces = srcAddress_[faceI];
                const scalarList& weights = srcWeights_[faceI];

                forAll(faces, i)
                {
                    cop(result[faceI], faceI, fld[faces[i]], weights[i]);
                }
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            srcAddress_.size(),
            pTraits<Type>::zero
        )
    );

    interpolateToSource
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult(),
        defaultValues
    );

    return tresult;
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const tmp<Field<Type> >& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), cop, defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            tgtAddress_.size(),
            pTraits<Type>::zero
        )
    );

    interpolateToTarget
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult(),
        defaultValues
    );

    return tresult;
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const tmp<Field<Type> >& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), cop, defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(fld, plusEqOp<Type>(), defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const tmp<Field<Type> >& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), plusEqOp<Type>(), defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(fld, plusEqOp<Type>(), defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const tmp<Field<Type> >& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), plusEqOp<Type>(), defaultValues);
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcPointFace
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const vector& n,
    const label tgtFaceI,
    point& tgtPoint
)
const
{
    const pointField& srcPoints = srcPatch.points();

    // source face addresses that intersect target face tgtFaceI
    const labelList& addr = tgtAddress_[tgtFaceI];

    forAll(addr, i)
    {
        label srcFaceI = addr[i];
        const face& f = srcPatch[srcFaceI];

        pointHit ray = f.ray(tgtPoint, n, srcPoints);

        if (ray.hit())
        {
            tgtPoint = ray.rawPoint();

            return srcFaceI;
        }
    }

    // no hit registered - try with face normal instead of input normal
    forAll(addr, i)
    {
        label srcFaceI = addr[i];
        const face& f = srcPatch[srcFaceI];

        vector nFace(-srcPatch.faceNormals()[srcFaceI]);
        nFace += tgtPatch.faceNormals()[tgtFaceI];

        pointHit ray = f.ray(tgtPoint, nFace, srcPoints);

        if (ray.hit())
        {
            tgtPoint = ray.rawPoint();

            return srcFaceI;
        }
    }

    return -1;
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtPointFace
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const vector& n,
    const label srcFaceI,
    point& srcPoint
)
const
{
    const pointField& tgtPoints = tgtPatch.points();

    // target face addresses that intersect source face srcFaceI
    const labelList& addr = srcAddress_[srcFaceI];

    forAll(addr, i)
    {
        label tgtFaceI = addr[i];
        const face& f = tgtPatch[tgtFaceI];

        pointHit ray = f.ray(srcPoint, n, tgtPoints);

        if (ray.hit())
        {
            srcPoint = ray.rawPoint();

            return tgtFaceI;
        }
    }

    // no hit registered - try with face normal instead of input normal
    forAll(addr, i)
    {
        label tgtFaceI = addr[i];
        const face& f = tgtPatch[tgtFaceI];

        vector nFace(-srcPatch.faceNormals()[srcFaceI]);
        nFace += tgtPatch.faceNormals()[tgtFaceI];

        pointHit ray = f.ray(srcPoint, n, tgtPoints);

        if (ray.hit())
        {
            srcPoint = ray.rawPoint();

            return tgtFaceI;
        }
    }

    return -1;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::writeFaceConnectivity
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const labelListList& srcAddress
)
const
{
    OFstream os("faceConnectivity" + Foam::name(Pstream::myProcNo()) + ".obj");

    label ptI = 1;

    forAll(srcAddress, i)
    {
        const labelList& addr = srcAddress[i];
        const point& srcPt = srcPatch.faceCentres()[i];
        forAll(addr, j)
        {
            label tgtPtI = addr[j];
            const point& tgtPt = tgtPatch.faceCentres()[tgtPtI];

            meshTools::writeOBJ(os, srcPt);
            meshTools::writeOBJ(os, tgtPt);

            os  << "l " << ptI << " " << ptI + 1 << endl;

            ptI += 2;
        }
    }
}


// ************************************************************************* //
