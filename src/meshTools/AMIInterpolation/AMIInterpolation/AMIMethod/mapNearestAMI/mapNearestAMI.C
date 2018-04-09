/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    const label& srcFacei,
    label& tgtFacei
) const
{
    const vectorField& srcCf = srcPatch.faceCentres();
    const vectorField& tgtCf = tgtPatch.faceCentres();

    const vector srcP = srcCf[srcFacei];

    DynamicList<label> tgtFaces(10);
    tgtFaces.append(tgtFacei);

    DynamicList<label> visitedFaces(10);

    scalar d = great;

    do
    {
        label tgtI = tgtFaces.remove();
        visitedFaces.append(tgtI);

        scalar dTest = magSqr(tgtCf[tgtI] - srcP);
        if (dTest < d)
        {
            tgtFacei = tgtI;
            d = dTest;

            this->appendNbrFaces
            (
                tgtFacei,
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
    label& srcFacei,
    label& tgtFacei
) const
{
    const labelList& srcNbr = this->srcPatch_.faceFaces()[srcFacei];

    srcFacei = -1;

    forAll(srcNbr, i)
    {
        label facei = srcNbr[i];
        if (mapFlag[facei])
        {
            srcFacei = facei;
            startSeedI = facei + 1;

            return;
        }
    }

    forAll(mapFlag, facei)
    {
        if (mapFlag[facei])
        {
            srcFacei = facei;
            tgtFacei = this->findTargetFace(facei);

            if (tgtFacei == -1)
            {
                const vectorField& srcCf = this->srcPatch_.faceCentres();

                FatalErrorInFunction
                    << "Unable to find target face for source face "
                    << srcFacei << " with face centre " << srcCf[srcFacei]
                    << abort(FatalError);
            }

            break;
        }
    }
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::mapNearestAMI<SourcePatch, TargetPatch>::findMappedSrcFace
(
    const label tgtFacei,
    const List<DynamicList<label>>& tgtToSrc
) const
{
    DynamicList<label> testFaces(10);
    DynamicList<label> visitedFaces(10);

    testFaces.append(tgtFacei);

    do
    {
        // search target tgtFacei neighbours for match with source face
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
    label srcFacei,
    label tgtFacei
)
{
    bool ok =
        this->initialise
        (
            srcAddress,
            srcWeights,
            tgtAddress,
            tgtWeights,
            srcFacei,
            tgtFacei
        );

    if (!ok)
    {
        return;
    }


    // temporary storage for addressing and weights
    List<DynamicList<label>> srcAddr(this->srcPatch_.size());
    List<DynamicList<label>> tgtAddr(this->tgtPatch_.size());


    // construct weights and addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // list to keep track of whether src face can be mapped
    boolList mapFlag(srcAddr.size(), true);

    // reset starting seed
    label startSeedI = 0;

    DynamicList<label> nonOverlapFaces;
    do
    {
        findNearestFace(this->srcPatch_, this->tgtPatch_, srcFacei, tgtFacei);

        srcAddr[srcFacei].append(tgtFacei);
        tgtAddr[tgtFacei].append(srcFacei);

        mapFlag[srcFacei] = false;

        // Do advancing front starting from srcFacei, tgtFacei
        setNextNearestFaces
        (
            mapFlag,
            startSeedI,
            srcFacei,
            tgtFacei
        );
    } while (srcFacei >= 0);


    // for the case of multiple source faces per target face, select the
    // nearest source face only and discard the others
    const vectorField& srcCf = this->srcPatch_.faceCentres();
    const vectorField& tgtCf = this->tgtPatch_.faceCentres();

    forAll(tgtAddr, targetFacei)
    {
        if (tgtAddr[targetFacei].size() > 1)
        {
            const vector& tgtC = tgtCf[tgtFacei];

            DynamicList<label>& srcFaces = tgtAddr[targetFacei];

            label srcFacei = srcFaces[0];
            scalar d = magSqr(tgtC - srcCf[srcFacei]);

            for (label i = 1; i < srcFaces.size(); i++)
            {
                label srcI = srcFaces[i];
                scalar dNew = magSqr(tgtC - srcCf[srcI]);
                if (dNew < d)
                {
                    d = dNew;
                    srcFacei = srcI;
                }
            }

            srcFaces.clear();
            srcFaces.append(srcFacei);
        }
    }

    // If there are more target faces than source faces, some target faces
    // might not yet be mapped
    forAll(tgtAddr, tgtFacei)
    {
        if (tgtAddr[tgtFacei].empty())
        {
            label srcFacei = findMappedSrcFace(tgtFacei, tgtAddr);

            if (srcFacei >= 0)
            {
                // note - reversed search from src->tgt to tgt->src
                findNearestFace
                (
                    this->tgtPatch_,
                    this->srcPatch_,
                    tgtFacei,
                    srcFacei
                );

                tgtAddr[tgtFacei].append(srcFacei);
            }
        }
    }


    // transfer data to persistent storage
    forAll(srcAddr, i)
    {
        srcAddress[i].transfer(srcAddr[i]);
        srcWeights[i] = scalarList(1, 1.0);
    }
    forAll(tgtAddr, i)
    {
        tgtAddress[i].transfer(tgtAddr[i]);
        tgtWeights[i] = scalarList(1, 1.0);
    }
}


// ************************************************************************* //
