/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "mapNearestAMI.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::mapNearestAMI<SourcePatch, TargetPatch>::findNearestFace
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const label& srcFaceI,
    label& tgtFaceI
) const
{
    const vectorField& srcCf = srcPatch.faceCentres();
    const vectorField& tgtCf = tgtPatch.faceCentres();

    const vector srcP = srcCf[srcFaceI];

    DynamicList<label> tgtFaces(10);
    tgtFaces.append(tgtFaceI);

    DynamicList<label> visitedFaces(10);

    scalar d = GREAT;

    do
    {
        label tgtI = tgtFaces.remove();
        visitedFaces.append(tgtI);

        scalar dTest = magSqr(tgtCf[tgtI] - srcP);
        if (dTest < d)
        {
            tgtFaceI = tgtI;
            d = dTest;

            this->appendNbrFaces
            (
                tgtFaceI,
                tgtPatch,
                visitedFaces,
                tgtFaces
            );
        }

    } while (tgtFaces.size() > 0);
}


template<class SourcePatch, class TargetPatch>
void Foam::mapNearestAMI<SourcePatch, TargetPatch>::setNextNearestFaces
(
    boolList& mapFlag,
    label& startSeedI,
    label& srcFaceI,
    label& tgtFaceI
) const
{
    const labelList& srcNbr = this->srcPatch_.faceFaces()[srcFaceI];

    srcFaceI = -1;

    forAll(srcNbr, i)
    {
        label faceI = srcNbr[i];
        if (mapFlag[faceI])
        {
            srcFaceI = faceI;
            startSeedI = faceI + 1;

            return;
        }
    }

    forAll(mapFlag, faceI)
    {
        if (mapFlag[faceI])
        {
            srcFaceI = faceI;
            tgtFaceI = this->findTargetFace(faceI);

            if (tgtFaceI == -1)
            {
                const vectorField& srcCf = this->srcPatch_.faceCentres();

                FatalErrorIn
                (
                    "void Foam::mapNearestAMI<SourcePatch, TargetPatch>::"
                    "setNextNearestFaces"
                    "("
                        "boolList&, "
                        "label&, "
                        "label&, "
                        "label&"
                    ") const"
                )
                    << "Unable to find target face for source face "
                    << srcFaceI << " with face centre " << srcCf[srcFaceI]
                    << abort(FatalError);
            }

            break;
        }
    }
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::mapNearestAMI<SourcePatch, TargetPatch>::findMappedSrcFace
(
    const label tgtFaceI,
    const List<DynamicList<label> >& tgtToSrc
) const
{
    DynamicList<label> testFaces(10);
    DynamicList<label> visitedFaces(10);

    testFaces.append(tgtFaceI);

    do
    {
        // search target tgtFaceI neighbours for match with source face
        label tgtI = testFaces.remove();

        if (findIndex(visitedFaces, tgtI) == -1)
        {
            visitedFaces.append(tgtI);

            if (tgtToSrc[tgtI].size())
            {
                return tgtToSrc[tgtI][0];
            }
            else
            {
                const labelList& nbrFaces = this->tgtPatch_.faceFaces()[tgtI];

                forAll(nbrFaces, i)
                {
                    if (findIndex(visitedFaces, nbrFaces[i]) == -1)
                    {
                        testFaces.append(nbrFaces[i]);
                    }
                }
            }
        }
    } while (testFaces.size());

    // did not find any match - should not be possible to get here!
    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::mapNearestAMI<SourcePatch, TargetPatch>::mapNearestAMI
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const scalarField& srcMagSf,
    const scalarField& tgtMagSf,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget,
    const bool requireMatch
)
:
    AMIMethod<SourcePatch, TargetPatch>
    (
        srcPatch,
        tgtPatch,
        srcMagSf,
        tgtMagSf,
        triMode,
        reverseTarget,
        requireMatch
    )
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::mapNearestAMI<SourcePatch, TargetPatch>::~mapNearestAMI()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::mapNearestAMI<SourcePatch, TargetPatch>::calculate
(
    labelListList& srcAddress,
    scalarListList& srcWeights,
    labelListList& tgtAddress,
    scalarListList& tgtWeights,
    label srcFaceI,
    label tgtFaceI
)
{
    bool ok =
        this->initialise
        (
            srcAddress,
            srcWeights,
            tgtAddress,
            tgtWeights,
            srcFaceI,
            tgtFaceI
        );

    if (!ok)
    {
        return;
    }


    // temporary storage for addressing and weights
    List<DynamicList<label> > srcAddr(this->srcPatch_.size());
    List<DynamicList<label> > tgtAddr(this->tgtPatch_.size());


    // construct weights and addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // list to keep track of whether src face can be mapped
    boolList mapFlag(srcAddr.size(), true);

    // reset starting seed
    label startSeedI = 0;

    DynamicList<label> nonOverlapFaces;
    do
    {
        findNearestFace(this->srcPatch_, this->tgtPatch_, srcFaceI, tgtFaceI);

        srcAddr[srcFaceI].append(tgtFaceI);
        tgtAddr[tgtFaceI].append(srcFaceI);

        mapFlag[srcFaceI] = false;

        // Do advancing front starting from srcFaceI, tgtFaceI
        setNextNearestFaces
        (
            mapFlag,
            startSeedI,
            srcFaceI,
            tgtFaceI
        );
    } while (srcFaceI >= 0);


    // for the case of multiple source faces per target face, select the
    // nearest source face only and discard the others
    const vectorField& srcCf = this->srcPatch_.faceCentres();
    const vectorField& tgtCf = this->tgtPatch_.faceCentres();

    forAll(tgtAddr, targetFaceI)
    {
        if (tgtAddr[targetFaceI].size() > 1)
        {
            const vector& tgtC = tgtCf[tgtFaceI];

            DynamicList<label>& srcFaces = tgtAddr[targetFaceI];

            label srcFaceI = srcFaces[0];
            scalar d = magSqr(tgtC - srcCf[srcFaceI]);

            for (label i = 1; i < srcFaces.size(); i++)
            {
                label srcI = srcFaces[i];
                scalar dNew = magSqr(tgtC - srcCf[srcI]);
                if (dNew < d)
                {
                    d = dNew;
                    srcFaceI = srcI;
                }
            }

            srcFaces.clear();
            srcFaces.append(srcFaceI);
        }
    }

    // If there are more target faces than source faces, some target faces
    // might not yet be mapped
    forAll(tgtAddr, tgtFaceI)
    {
        if (tgtAddr[tgtFaceI].empty())
        {
            label srcFaceI = findMappedSrcFace(tgtFaceI, tgtAddr);

            if (srcFaceI >= 0)
            {
                // note - reversed search from src->tgt to tgt->src
                findNearestFace
                (
                    this->tgtPatch_,
                    this->srcPatch_,
                    tgtFaceI,
                    srcFaceI
                );

                tgtAddr[tgtFaceI].append(srcFaceI);
            }
        }
    }


    // transfer data to persistent storage
    forAll(srcAddr, i)
    {
        scalar magSf = this->srcMagSf_[i];
        srcAddress[i].transfer(srcAddr[i]);
        srcWeights[i] = scalarList(1, magSf);
    }
    forAll(tgtAddr, i)
    {
        scalar magSf = this->tgtMagSf_[i];
        tgtAddress[i].transfer(tgtAddr[i]);
        tgtWeights[i] = scalarList(1, magSf);
    }
}


// ************************************************************************* //
