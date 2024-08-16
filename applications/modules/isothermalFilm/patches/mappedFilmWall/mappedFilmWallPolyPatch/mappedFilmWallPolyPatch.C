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

#include "mappedFilmWallPolyPatch.H"
#include "mappedPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedFilmWallPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, mappedFilmWallPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        mappedFilmWallPolyPatch,
        dictionary
    );
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::mappedFilmWallPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    filmWallPolyPatch::calcGeometry(pBufs);
    mappedPatchBase::clearOut(false);
}


void Foam::mappedFilmWallPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    filmWallPolyPatch::movePoints(pBufs, p);
    mappedPatchBase::clearOut(true);
}


void Foam::mappedFilmWallPolyPatch::topoChange(PstreamBuffers& pBufs)
{
    filmWallPolyPatch::topoChange(pBufs);
    mappedPatchBase::clearOut(false);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedFilmWallPolyPatch::mappedFilmWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    filmWallPolyPatch(name, size, start, index, bm, patchType),
    mappedPatchBase(static_cast<const polyPatch&>(*this))
{
    //  mapped is not constraint type so add mapped group explicitly
    if (findIndex(inGroups(), mappedPolyPatch::typeName) == -1)
    {
        inGroups().append(mappedPolyPatch::typeName);
    }
}


Foam::mappedFilmWallPolyPatch::mappedFilmWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& neighbourRegion,
    const word& neighbourPatch,
    const polyBoundaryMesh& bm
)
:
    filmWallPolyPatch(name, size, start, index, bm, typeName),
    mappedPatchBase
    (
        *this,
        neighbourRegion,
        neighbourPatch,
        cyclicTransform(true)
    )
{}


Foam::mappedFilmWallPolyPatch::mappedFilmWallPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    filmWallPolyPatch(name, dict, index, bm, patchType),
    mappedPatchBase(*this, dict, transformType::none)
{
    //  mapped is not constraint type so add mapped group explicitly
    if (findIndex(inGroups(), mappedPolyPatch::typeName) == -1)
    {
        inGroups().append(mappedPolyPatch::typeName);
    }
}


Foam::mappedFilmWallPolyPatch::mappedFilmWallPolyPatch
(
    const mappedFilmWallPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    filmWallPolyPatch(pp, bm),
    mappedPatchBase(*this, pp)
{}


Foam::mappedFilmWallPolyPatch::mappedFilmWallPolyPatch
(
    const mappedFilmWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    filmWallPolyPatch(pp, bm, index, newSize, newStart),
    mappedPatchBase(*this, pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedFilmWallPolyPatch::~mappedFilmWallPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedFilmWallPolyPatch::write(Ostream& os) const
{
    filmWallPolyPatch::write(os);
    mappedPatchBase::write(os);
}


// ************************************************************************* //
