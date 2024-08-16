/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "mappedFilmSurfacePolyPatch.H"
#include "mappedPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedFilmSurfacePolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, mappedFilmSurfacePolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        mappedFilmSurfacePolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::mappedFilmSurfacePolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    filmSurfacePolyPatch::calcGeometry(pBufs);
    mappedPatchBase::clearOut(false);
}


void Foam::mappedFilmSurfacePolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    filmSurfacePolyPatch::movePoints(pBufs, p);
    mappedPatchBase::clearOut(true);
}


void Foam::mappedFilmSurfacePolyPatch::topoChange(PstreamBuffers& pBufs)
{
    filmSurfacePolyPatch::topoChange(pBufs);
    mappedPatchBase::clearOut(false);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedFilmSurfacePolyPatch::mappedFilmSurfacePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    filmSurfacePolyPatch(name, size, start, index, bm, patchType),
    mappedExtrudedPatchBase(static_cast<const polyPatch&>(*this))
{
    //  mapped is not constraint type so add mapped group explicitly
    if (findIndex(inGroups(), mappedPolyPatch::typeName) == -1)
    {
        inGroups().append(mappedPolyPatch::typeName);
    }
}


Foam::mappedFilmSurfacePolyPatch::mappedFilmSurfacePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& neighbourRegion,
    const word& neighbourPatch,
    const bool isExtrudedRegion,
    const polyBoundaryMesh& bm
)
:
    filmSurfacePolyPatch(name, size, start, index, bm, typeName),
    mappedExtrudedPatchBase
    (
        *this,
        neighbourRegion,
        neighbourPatch,
        isExtrudedRegion,
        cyclicTransform(true)
    )
{}


Foam::mappedFilmSurfacePolyPatch::mappedFilmSurfacePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    filmSurfacePolyPatch(name, dict, index, bm, patchType),
    mappedExtrudedPatchBase(*this, dict)
{
    //  mapped is not constraint type so add mapped group explicitly
    if (findIndex(inGroups(), mappedPolyPatch::typeName) == -1)
    {
        inGroups().append(mappedPolyPatch::typeName);
    }
}


Foam::mappedFilmSurfacePolyPatch::mappedFilmSurfacePolyPatch
(
    const mappedFilmSurfacePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    filmSurfacePolyPatch(pp, bm),
    mappedExtrudedPatchBase(*this, pp)
{}


Foam::mappedFilmSurfacePolyPatch::mappedFilmSurfacePolyPatch
(
    const mappedFilmSurfacePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    filmSurfacePolyPatch(pp, bm, index, newSize, newStart),
    mappedExtrudedPatchBase(*this, pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedFilmSurfacePolyPatch::~mappedFilmSurfacePolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedFilmSurfacePolyPatch::write(Ostream& os) const
{
    filmSurfacePolyPatch::write(os);
    mappedExtrudedPatchBase::write(os);
}


// ************************************************************************* //
