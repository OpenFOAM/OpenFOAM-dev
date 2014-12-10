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

const Foam::scalar Foam::cyclicACMIPolyPatch::tolerance_ = 1e-6;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::cyclicACMIPolyPatch::initPatchFaceAreas() const
{
    if (!empty() && faceAreas0_.empty())
    {
        faceAreas0_ = faceAreas();
    }

    const cyclicACMIPolyPatch& nbrACMI =
        refCast<const cyclicACMIPolyPatch>(this->neighbPatch());

    if (!nbrACMI.empty() && nbrACMI.faceAreas0().empty())
    {
        nbrACMI.initPatchFaceAreas();
    }
}


void Foam::cyclicACMIPolyPatch::resetAMI
(
    const AMIPatchToPatchInterpolation::interpolationMethod&
) const
{
    if (owner())
    {
        const polyPatch& nonOverlapPatch = this->nonOverlapPatch();

        initPatchFaceAreas();

        // reset patch face areas based on original patch for AMI calculation
        vectorField::subField Sf = faceAreas();
        vectorField::subField noSf = nonOverlapPatch.faceAreas();

        forAll(Sf, faceI)
        {
            Sf[faceI] = faceAreas0_[faceI];
            noSf[faceI] = faceAreas0_[faceI];
        }

        // calculate the AMI using partial face-area-weighted
        cyclicAMIPolyPatch::resetAMI
        (
            AMIPatchToPatchInterpolation::imPartialFaceAreaWeight
        );

        srcMask_ =
            min(scalar(1) - tolerance_, max(tolerance_, AMI().srcWeightsSum()));

        tgtMask_ =
            min(scalar(1) - tolerance_, max(tolerance_, AMI().tgtWeightsSum()));

        forAll(Sf, faceI)
        {
            Sf[faceI] *= srcMask_[faceI];
            noSf[faceI] *= 1.0 - srcMask_[faceI];
        }

        setNeighbourFaceAreas();

        // set the updated flag
        updated_ = true;
    }
}


void Foam::cyclicACMIPolyPatch::setNeighbourFaceAreas() const
{
    const cyclicACMIPolyPatch& cp =
        refCast<const cyclicACMIPolyPatch>(this->neighbPatch());
    const polyPatch& pp = cp.nonOverlapPatch();

    const vectorField& faceAreas0 = cp.faceAreas0();

    vectorField::subField Sf = cp.faceAreas();
    vectorField::subField noSf = pp.faceAreas();

    forAll(Sf, faceI)
    {
        Sf[faceI] = tgtMask_[faceI]*faceAreas0[faceI];
        noSf[faceI] = (1.0 - tgtMask_[faceI])*faceAreas0[faceI];
    }
}


void Foam::cyclicACMIPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    // initialise the AMI so that base geometry (e.g. cell volumes) are
    // correctly evaluated
    resetAMI();

    cyclicAMIPolyPatch::initGeometry(pBufs);
}


void Foam::cyclicACMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    cyclicAMIPolyPatch::calcGeometry(pBufs);
}


void Foam::cyclicACMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    cyclicAMIPolyPatch::initMovePoints(pBufs, p);
}


void Foam::cyclicACMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    cyclicAMIPolyPatch::movePoints(pBufs, p);
}


void Foam::cyclicACMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    cyclicAMIPolyPatch::initUpdateMesh(pBufs);
}


void Foam::cyclicACMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    cyclicAMIPolyPatch::updateMesh(pBufs);
}


void Foam::cyclicACMIPolyPatch::clearGeom()
{
    cyclicAMIPolyPatch::clearGeom();
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
    const word& patchType,
    const transformType transform
)
:
    cyclicAMIPolyPatch(name, size, start, index, bm, patchType, transform),
    faceAreas0_(),
    nonOverlapPatchName_(word::null),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    AMIRequireMatch_ = false;

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
    cyclicAMIPolyPatch(name, dict, index, bm, patchType),
    faceAreas0_(),
    nonOverlapPatchName_(dict.lookup("nonOverlapPatch")),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    AMIRequireMatch_ = false;

    if (nonOverlapPatchName_ == name)
    {
        FatalIOErrorIn
        (
            "cyclicACMIPolyPatch::cyclicACMIPolyPatch"
            "("
                "const word&, "
                "const dictionary&, "
                "const label, "
                "const polyBoundaryMesh&, "
                "const word&"
            ")",
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
    faceAreas0_(),
    nonOverlapPatchName_(pp.nonOverlapPatchName_),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    AMIRequireMatch_ = false;

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
    faceAreas0_(),
    nonOverlapPatchName_(nonOverlapPatchName),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    AMIRequireMatch_ = false;

    if (nonOverlapPatchName_ == name())
    {
        FatalErrorIn
        (
            "const cyclicACMIPolyPatch& "
            "const polyBoundaryMesh&, "
            "const label, "
            "const label, "
            "const label, "
            "const word&, "
            "const word&"
        )   << "Non-overlapping patch name " << nonOverlapPatchName_
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
    faceAreas0_(),
    nonOverlapPatchName_(pp.nonOverlapPatchName_),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    AMIRequireMatch_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicACMIPolyPatch::~cyclicACMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::cyclicACMIPolyPatch& Foam::cyclicACMIPolyPatch::neighbPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[neighbPatchID()];
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
            FatalErrorIn("cyclicPolyAMIPatch::neighbPatchID() const")
                << "Illegal non-overlapping patch name " << nonOverlapPatchName_
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        if (nonOverlapPatchID_ < index())
        {
            FatalErrorIn("cyclicPolyAMIPatch::neighbPatchID() const")
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
            const scalarField magSf(mag(faceAreas()));
            const scalarField noMagSf(mag(noPp.faceAreas()));

            forAll(magSf, faceI)
            {
                scalar ratio = mag(magSf[faceI]/(noMagSf[faceI] + ROOTVSMALL));

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
            FatalErrorIn
            (
                "Foam::label "
                "Foam::cyclicACMIPolyPatch::nonOverlapPatchID() const"
            )   << "Inconsistent ACMI patches " << name() << " and "
                << noPp.name() << ".  Patches should have identical topology"
                << exit(FatalError);
        }
    }

    return nonOverlapPatchID_;
}


void Foam::cyclicACMIPolyPatch::calcGeometry
(
    const primitivePatch& referPatch,
    const pointField& thisCtrs,
    const vectorField& thisAreas,
    const pointField& thisCc,
    const pointField& nbrCtrs,
    const vectorField& nbrAreas,
    const pointField& nbrCc
)
{
    cyclicAMIPolyPatch::calcGeometry
    (
        referPatch,
        thisCtrs,
        thisAreas,
        thisCc,
        nbrCtrs,
        nbrAreas,
        nbrCc
    );
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

    os.writeKeyword("nonOverlapPatch") << nonOverlapPatchName_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
