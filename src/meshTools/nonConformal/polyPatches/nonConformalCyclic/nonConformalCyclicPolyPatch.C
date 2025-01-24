/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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

#include "nonConformalCyclicPolyPatch.H"
#include "nonConformalErrorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "mappedPatchBaseBase.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "SubField.H"
#include "nonConformalBoundary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalCyclicPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, nonConformalCyclicPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        nonConformalCyclicPolyPatch,
        dictionary
    );
}


namespace Foam
{
    template<>
    const char*
        NamedEnum<nonConformalCyclicPolyPatch::moveUpdate, 3>::names[] =
        {"always", "detect", "never"};
}

const Foam::NamedEnum<Foam::nonConformalCyclicPolyPatch::moveUpdate, 3>
    Foam::nonConformalCyclicPolyPatch::moveUpdateNames_;


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::nonConformalCyclicPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    cyclicPolyPatch::initCalcGeometry(pBufs);
    intersectionIsValid_ = false;
    raysIsValid_ = false;
}


void Foam::nonConformalCyclicPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    static_cast<cyclicTransform&>(*this) =
        cyclicTransform
        (
            name(),
            origPatch().faceAreas(),
            *this,
            nbrPatchName(),
            nbrPatch(),
            matchTolerance()
        );
}


void Foam::nonConformalCyclicPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    cyclicPolyPatch::initMovePoints(pBufs, p);

    if (moveUpdate_ == moveUpdate::never)
    {
        // Do nothing
    }
    else if (moveUpdate_ == moveUpdate::detect)
    {
        intersectionIsValid_ = min(intersectionIsValid_, 1);
        raysIsValid_ = min(raysIsValid_, 1);
    }
    else
    {
        intersectionIsValid_ = 0;
        raysIsValid_ = 0;
    }
}


void Foam::nonConformalCyclicPolyPatch::initTopoChange(PstreamBuffers& pBufs)
{
    cyclicPolyPatch::initTopoChange(pBufs);
    intersectionIsValid_ = 0;
    raysIsValid_ = 0;
}


void Foam::nonConformalCyclicPolyPatch::clearGeom()
{
    cyclicPolyPatch::clearGeom();
    intersectionIsValid_ = 0;
    raysIsValid_ = 0;
}


void Foam::nonConformalCyclicPolyPatch::rename(const wordList& newNames)
{
    cyclicPolyPatch::rename(newNames);
    nonConformalCoupledPolyPatch::rename(newNames);
}


void Foam::nonConformalCyclicPolyPatch::reorder(const labelUList& newToOldIndex)
{
    cyclicPolyPatch::reorder(newToOldIndex);
    nonConformalCoupledPolyPatch::reorder(newToOldIndex);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalCyclicPolyPatch::nonConformalCyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicPolyPatch(name, size, start, index, bm, patchType),
    nonConformalCoupledPolyPatch(static_cast<const polyPatch&>(*this)),
    intersectionIsValid_(0),
    intersection_(false),
    raysIsValid_(0),
    rays_(false),
    moveUpdate_(moveUpdate::always)
{}


Foam::nonConformalCyclicPolyPatch::nonConformalCyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const word& nbrPatchName,
    const word& origPatchName,
    const cyclicTransform& transform
)
:
    cyclicPolyPatch
    (
        name,
        size,
        start,
        index,
        bm,
        patchType,
        nbrPatchName,
        transform
    ),
    nonConformalCoupledPolyPatch(*this, origPatchName),
    intersectionIsValid_(0),
    intersection_(false),
    raysIsValid_(0),
    rays_(false),
    moveUpdate_(moveUpdate::always)
{}


Foam::nonConformalCyclicPolyPatch::nonConformalCyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicPolyPatch(name, dict, index, bm, patchType, true),
    nonConformalCoupledPolyPatch(*this, dict),
    intersectionIsValid_(0),
    intersection_(false),
    raysIsValid_(0),
    rays_(false),
    moveUpdate_
    (
        dict.found("moveUpdate")
      ? moveUpdateNames_.read(dict.lookup("moveUpdate"))
      : dict.found("reMapAfterMove") // <-- backwards compatibility
      ? (
            dict.lookup<bool>("reMapAfterMove")
          ? moveUpdate::always
          : moveUpdate::never
        )
      : moveUpdate::always
    )
{}


Foam::nonConformalCyclicPolyPatch::nonConformalCyclicPolyPatch
(
    const nonConformalCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicPolyPatch(pp, bm),
    nonConformalCoupledPolyPatch(*this, pp),
    intersectionIsValid_(0),
    intersection_(false),
    raysIsValid_(0),
    rays_(false),
    moveUpdate_(pp.moveUpdate_)
{}


Foam::nonConformalCyclicPolyPatch::nonConformalCyclicPolyPatch
(
    const nonConformalCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName,
    const word& origPatchName
)
:
    cyclicPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    nonConformalCoupledPolyPatch(*this, origPatchName),
    intersectionIsValid_(0),
    intersection_(false),
    raysIsValid_(0),
    rays_(false),
    moveUpdate_(pp.moveUpdate_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalCyclicPolyPatch::~nonConformalCyclicPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::nonConformalCyclicPolyPatch&
Foam::nonConformalCyclicPolyPatch::nbrPatch() const
{
    return
        refCast<const nonConformalCyclicPolyPatch>
        (
            cyclicPolyPatch::nbrPatch()
        );
}


bool Foam::nonConformalCyclicPolyPatch::coupled() const
{
    return false;
}


const Foam::patchToPatches::intersection&
Foam::nonConformalCyclicPolyPatch::intersection() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "The non-conformal cyclic intersection is only available to "
            << "the owner patch" << abort(FatalError);
    }

    const bool intersectionIsValid =
        (intersectionIsValid_ == 2)
     || (
            intersectionIsValid_ == 1
         && !mappedPatchBaseBase::moving
            (
                origPatch(),
                nbrPatch().origPatch()
            )
        );

    if (!intersectionIsValid)
    {
        const polyMesh& mesh = boundaryMesh().mesh();

        const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh);

        const string inRegionName =
            mesh.name() == polyMesh::defaultRegion
          ? ""
          : " in region " + mesh.name();

        intersection_.update
        (
            origPatch(),
            ncb.patchPointNormals(origPatchIndex()),
            nbrPatch().origPatch(),
            transform(),
            {
                origPatchName() + inRegionName,
                nbrPatch().origPatchName() + inRegionName
            },
            cyclicTransform::str()
        );

        intersectionIsValid_ = 2;
    }

    return intersection_;
}


const Foam::patchToPatches::rays&
Foam::nonConformalCyclicPolyPatch::rays() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "The non-conformal cyclic rays is only available to "
            << "the owner patch" << abort(FatalError);
    }

    const bool raysIsValid =
        (raysIsValid_ == 2)
     || (
            raysIsValid_ == 1
         && !mappedPatchBaseBase::moving
            (
                origPatch(),
                nbrPatch().origPatch()
            )
        );

    if (!raysIsValid)
    {
        const polyMesh& mesh = boundaryMesh().mesh();

        const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh);

        const string inRegionName =
            mesh.name() == polyMesh::defaultRegion
          ? ""
          : " in region " + mesh.name();

        rays_.update
        (
            primitiveOldTimePatch
            (
                origPatch(),
                mesh.points(),
                mesh.oldPoints()
            ),
            ncb.patchPointNormals(origPatchIndex()),
            ncb.patchPointNormals0(origPatchIndex()),
            primitiveOldTimePatch
            (
                nbrPatch().origPatch(),
                mesh.points(),
                mesh.oldPoints()
            ),
            transform(),
            {
                origPatchName() + inRegionName,
                nbrPatch().origPatchName() + inRegionName
            },
            cyclicTransform::str()
        );

        raysIsValid_ = 2;
    }

    return rays_;
}


Foam::remote Foam::nonConformalCyclicPolyPatch::ray
(
    const scalar fraction,
    const label origFacei,
    const vector& p,
    const vector& n,
    point& nbrP
) const
{
    const polyMesh& mesh = boundaryMesh().mesh();

    const nonConformalCyclicPolyPatch& ownerPatch =
        owner() ? *this : nbrPatch();

    auto ownerRaysMethod =
        owner()
      ? &patchToPatches::rays::srcToTgtRay
      : &patchToPatches::rays::tgtToSrcRay;

    return
        (ownerPatch.rays().*ownerRaysMethod)
        (
            primitiveOldTimePatch
            (
                nbrPatch().origPatch(),
                mesh.points(),
                mesh.oldPoints()
            ),
            fraction,
            origFacei,
            transform().invTransformPosition(p),
            transform().invTransform(n),
            nbrP
        );
}


void Foam::nonConformalCyclicPolyPatch::write(Ostream& os) const
{
    cyclicPolyPatch::write(os);
    nonConformalCoupledPolyPatch::write(os);
    writeEntryIfDifferent<word>
    (
        os,
        "moveUpdate",
        moveUpdateNames_[moveUpdate::always],
        moveUpdateNames_[moveUpdate_]
    );
}


// ************************************************************************* //
