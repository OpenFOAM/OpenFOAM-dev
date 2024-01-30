/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "nonConformalMappedWallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalMappedWallPolyPatch, 0);

    addToRunTimeSelectionTable
    (
        polyPatch,
        nonConformalMappedWallPolyPatch,
        word
    );
    addToRunTimeSelectionTable
    (
        polyPatch,
        nonConformalMappedWallPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::nonConformalMappedWallPolyPatch::initCalcGeometry
(
    PstreamBuffers& pBufs
)
{
    wallPolyPatch::initCalcGeometry(pBufs);
    nonConformalMappedPatchBase::clearOut();
}


void Foam::nonConformalMappedWallPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)

{
    wallPolyPatch::initMovePoints(pBufs, p);
    nonConformalMappedPatchBase::clearOut();
}


void Foam::nonConformalMappedWallPolyPatch::initTopoChange
(
    PstreamBuffers& pBufs
)
{
    wallPolyPatch::initTopoChange(pBufs);
    nonConformalMappedPatchBase::clearOut();
}


void Foam::nonConformalMappedWallPolyPatch::clearGeom()
{
    wallPolyPatch::clearGeom();
    nonConformalMappedPatchBase::clearOut();
}


void Foam::nonConformalMappedWallPolyPatch::rename(const wordList& newNames)
{
    wallPolyPatch::rename(newNames);
    nonConformalPolyPatch::rename(newNames);
}


void Foam::nonConformalMappedWallPolyPatch::reorder
(
    const labelUList& newToOldIndex
)
{
    wallPolyPatch::reorder(newToOldIndex);
    nonConformalPolyPatch::reorder(newToOldIndex);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalMappedWallPolyPatch::nonConformalMappedWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    wallPolyPatch(name, size, start, index, bm, patchType),
    nonConformalPolyPatch(static_cast<const polyPatch&>(*this)),
    nonConformalMappedPatchBase
    (
        static_cast<const nonConformalPolyPatch&>(*this)
    )
{}


Foam::nonConformalMappedWallPolyPatch::nonConformalMappedWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& origPatchName,
    const word& nbrRegionName,
    const word& nbrPatchName,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(name, size, start, index, bm, typeName),
    nonConformalPolyPatch(*this, origPatchName),
    nonConformalMappedPatchBase
    (
        *this,
        nbrRegionName,
        nbrPatchName,
        cyclicTransform(true)
    )
{}


Foam::nonConformalMappedWallPolyPatch::nonConformalMappedWallPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    wallPolyPatch(name, dict, index, bm, patchType),
    nonConformalPolyPatch(*this, dict),
    nonConformalMappedPatchBase(*this, dict, true)
{}


Foam::nonConformalMappedWallPolyPatch::nonConformalMappedWallPolyPatch
(
    const nonConformalMappedWallPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(pp, bm),
    nonConformalPolyPatch(*this, pp),
    nonConformalMappedPatchBase(*this, pp)
{}


Foam::nonConformalMappedWallPolyPatch::nonConformalMappedWallPolyPatch
(
    const nonConformalMappedWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    wallPolyPatch(pp, bm, index, newSize, newStart),
    nonConformalPolyPatch(*this, pp),
    nonConformalMappedPatchBase(*this, pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalMappedWallPolyPatch::~nonConformalMappedWallPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nonConformalMappedWallPolyPatch::write(Ostream& os) const
{
    wallPolyPatch::write(os);
    nonConformalPolyPatch::write(os);
    nonConformalMappedPatchBase::write(os);
}


// ************************************************************************* //
