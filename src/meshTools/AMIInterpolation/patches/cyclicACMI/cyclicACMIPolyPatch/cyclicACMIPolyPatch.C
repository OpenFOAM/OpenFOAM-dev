/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "cyclicACMIPolyPatch.H"
#include "SubField.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicACMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicACMIPolyPatch, dictionary);
}

const Foam::scalar Foam::cyclicACMIPolyPatch::tolerance_ = 1e-10;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::cyclicACMIPolyPatch::resetAMI() const
{
    if (owner())
    {
        const polyPatch& nonOverlapPatch = this->nonOverlapPatch();

        if (debug)
        {
            Pout<< "cyclicACMIPolyPatch::resetAMI : recalculating weights"
                << " for " << name() << " and " << nonOverlapPatch.name()
                << endl;
        }

        if (boundaryMesh().mesh().hasCellCentres())
        {
            if (debug)
            {
                Pout<< "cyclicACMIPolyPatch::resetAMI : clearing cellCentres"
                    << " for " << name() << " and " << nonOverlapPatch.name()
                    << endl;
            }

            const_cast<polyMesh&>
            (
                boundaryMesh().mesh()
            ).primitiveMesh::clearGeom();
        }

        // Trigger re-building of faceAreas
        (void)boundaryMesh().mesh().faceAreas();


        // Calculate the AMI using partial face-area-weighted. This leaves
        // the weights as fractions of local areas (sum(weights) = 1 means
        // face is fully covered)
        cyclicAMIPolyPatch::resetAMI();

        AMIInterpolation& AMI = this->AMIs_[0];

        srcMask_ =
            min(scalar(1) - tolerance_, max(tolerance_, AMI.srcWeightsSum()));

        tgtMask_ =
            min(scalar(1) - tolerance_, max(tolerance_, AMI.tgtWeightsSum()));


        // Adapt owner side areas. Note that in uncoupled situations (e.g.
        // decomposePar) srcMask, tgtMask can be zero size.
        if (srcMask_.size())
        {
            vectorField::subField Sf = faceAreas();
            scalarField::subField magSf = magFaceAreas();
            vectorField::subField noSf = nonOverlapPatch.faceAreas();
            scalarField::subField noMagSf = nonOverlapPatch.magFaceAreas();

            forAll(Sf, facei)
            {
                Sf[facei] *= srcMask_[facei];
                magSf[facei] = mag(Sf[facei]);
                noSf[facei] *= 1.0 - srcMask_[facei];
                noMagSf[facei] = mag(noSf[facei]);
            }
        }
        // Adapt slave side areas
        if (tgtMask_.size())
        {
            const cyclicACMIPolyPatch& cp =
                refCast<const cyclicACMIPolyPatch>(this->nbrPatch());
            const polyPatch& pp = cp.nonOverlapPatch();

            vectorField::subField Sf = cp.faceAreas();
            scalarField::subField magSf = cp.magFaceAreas();
            vectorField::subField noSf = pp.faceAreas();
            scalarField::subField noMagSf = pp.magFaceAreas();

            forAll(Sf, facei)
            {
                Sf[facei] *= tgtMask_[facei];
                magSf[facei] = mag(Sf[facei]);
                noSf[facei] *= 1.0 - tgtMask_[facei];
                noMagSf[facei] = mag(noSf[facei]);
            }
        }

        // Re-normalise the weights since the effect of overlap is already
        // accounted for in the area.
        {
            scalarListList& srcWeights = AMI.srcWeights();
            scalarField& srcWeightsSum = AMI.srcWeightsSum();
            forAll(srcWeights, i)
            {
                scalarList& wghts = srcWeights[i];
                if (wghts.size())
                {
                    scalar& sum = srcWeightsSum[i];

                    forAll(wghts, j)
                    {
                        wghts[j] /= sum;
                    }
                    sum = 1.0;
                }
            }
        }
        {
            scalarListList& tgtWeights = AMI.tgtWeights();
            scalarField& tgtWeightsSum = AMI.tgtWeightsSum();
            forAll(tgtWeights, i)
            {
                scalarList& wghts = tgtWeights[i];
                if (wghts.size())
                {
                    scalar& sum = tgtWeightsSum[i];
                    forAll(wghts, j)
                    {
                        wghts[j] /= sum;
                    }
                    sum = 1.0;
                }
            }
        }

        // Set the updated flag
        updated_ = true;
    }
}


void Foam::cyclicACMIPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    cyclicAMIPolyPatch::initCalcGeometry(pBufs);

    // Initialise the AMI
    resetAMI();
}


void Foam::cyclicACMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    static_cast<cyclicTransform&>(*this) =
        cyclicTransform
        (
            name(),
            faceAreas(),
            *this,
            nbrPatchName(),
            nbrPatch(),
            matchTolerance()
        );
}


void Foam::cyclicACMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    cyclicAMIPolyPatch::initMovePoints(pBufs, p);

    // Initialise the AMI
    resetAMI();
}


const Foam::scalarField& Foam::cyclicACMIPolyPatch::srcMask() const
{
    return srcMask_;
}


const Foam::scalarField& Foam::cyclicACMIPolyPatch::tgtMask() const
{
    return tgtMask_;
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicAMIPolyPatch
    (
        name,
        size,
        start,
        index,
        bm,
        patchType,
        false,
        AMIInterpolation::imPartialFaceAreaWeight
    ),
    nonOverlapPatchName_(word::null),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicAMIPolyPatch
    (
        name,
        dict,
        index,
        bm,
        patchType,
        false,
        AMIInterpolation::imPartialFaceAreaWeight
    ),
    nonOverlapPatchName_(dict.lookup("nonOverlapPatch")),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    if (nonOverlapPatchName_ == name)
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Non-overlapping patch name " << nonOverlapPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicAMIPolyPatch(pp, bm),
    nonOverlapPatchName_(pp.nonOverlapPatchName_),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName,
    const word& nonOverlapPatchName
)
:
    cyclicAMIPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    nonOverlapPatchName_(nonOverlapPatchName),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    if (nonOverlapPatchName_ == name())
    {
        FatalErrorInFunction
            << "Non-overlapping patch name " << nonOverlapPatchName_
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicAMIPolyPatch(pp, bm, index, mapAddressing, newStart),
    nonOverlapPatchName_(pp.nonOverlapPatchName_),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicACMIPolyPatch::~cyclicACMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::cyclicACMIPolyPatch& Foam::cyclicACMIPolyPatch::nbrPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[nbrPatchID()];
    return refCast<const cyclicACMIPolyPatch>(pp);
}


Foam::label Foam::cyclicACMIPolyPatch::nonOverlapPatchID() const
{
    if (nonOverlapPatchID_ == -1)
    {
        nonOverlapPatchID_ =
            this->boundaryMesh().findPatchID(nonOverlapPatchName_);

        if (nonOverlapPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal non-overlapping patch name " << nonOverlapPatchName_
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        if (nonOverlapPatchID_ < index())
        {
            FatalErrorInFunction
                << "Boundary ordering error: " << type()
                << " patch must be defined prior to its non-overlapping patch"
                << nl
                << type() << " patch: " << name() << ", ID:" << index() << nl
                << "Non-overlap patch: " << nonOverlapPatchName_
                << ", ID:" << nonOverlapPatchID_ << nl
                << exit(FatalError);
        }

        const polyPatch& noPp = this->boundaryMesh()[nonOverlapPatchID_];

        bool ok = true;

        if (size() == noPp.size())
        {
            const scalarField magSf(magFaceAreas());
            const scalarField noMagSf(noPp.magFaceAreas());

            forAll(magSf, facei)
            {
                scalar ratio = mag(magSf[facei]/(noMagSf[facei] + rootVSmall));

                if (ratio - 1 > tolerance_)
                {
                    ok = false;
                    break;
                }
            }
        }
        else
        {
            ok = false;
        }

        if (!ok)
        {
            FatalErrorInFunction
                << "Inconsistent ACMI patches " << name() << " and "
                << noPp.name() << ".  Patches should have identical topology"
                << exit(FatalError);
        }
    }

    return nonOverlapPatchID_;
}


void Foam::cyclicACMIPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{
    cyclicAMIPolyPatch::initOrder(pBufs, pp);
}


bool Foam::cyclicACMIPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    return cyclicAMIPolyPatch::order(pBufs, pp, faceMap, rotation);
}


void Foam::cyclicACMIPolyPatch::write(Ostream& os) const
{
    cyclicAMIPolyPatch::write(os);
    writeEntry(os, "nonOverlapPatch", nonOverlapPatchName_);
}


// ************************************************************************* //
