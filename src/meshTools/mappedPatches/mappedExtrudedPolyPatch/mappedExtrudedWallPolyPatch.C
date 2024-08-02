/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "mappedExtrudedWallPolyPatch.H"
#include "mappedPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedExtrudedWallPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, mappedExtrudedWallPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        mappedExtrudedWallPolyPatch,
        dictionary
    );
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::mappedExtrudedWallPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    wallPolyPatch::calcGeometry(pBufs);
    mappedPatchBase::clearOut(false);
}


void Foam::mappedExtrudedWallPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    wallPolyPatch::movePoints(pBufs, p);
    mappedPatchBase::clearOut(true);
}


void Foam::mappedExtrudedWallPolyPatch::topoChange(PstreamBuffers& pBufs)
{
    wallPolyPatch::topoChange(pBufs);
    mappedPatchBase::clearOut(false);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedExtrudedWallPolyPatch::mappedExtrudedWallPolyPatch
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
    mappedExtrudedPatchBase(static_cast<const polyPatch&>(*this))
{
    //  mapped is not constraint type so add mapped group explicitly
    if (findIndex(inGroups(), mappedPolyPatch::typeName) == -1)
    {
        inGroups().append(mappedPolyPatch::typeName);
    }
}


Foam::mappedExtrudedWallPolyPatch::mappedExtrudedWallPolyPatch
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
    wallPolyPatch(name, size, start, index, bm, typeName),
    mappedExtrudedPatchBase
    (
        *this,
        neighbourRegion,
        neighbourPatch,
        isExtrudedRegion,
        cyclicTransform(true)
    )
{}


Foam::mappedExtrudedWallPolyPatch::mappedExtrudedWallPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    wallPolyPatch(name, dict, index, bm, patchType),
    mappedExtrudedPatchBase(*this, dict, true)
{
    //  mapped is not constraint type so add mapped group explicitly
    if (findIndex(inGroups(), mappedPolyPatch::typeName) == -1)
    {
        inGroups().append(mappedPolyPatch::typeName);
    }
}


Foam::mappedExtrudedWallPolyPatch::mappedExtrudedWallPolyPatch
(
    const mappedExtrudedWallPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(pp, bm),
    mappedExtrudedPatchBase(*this, pp)
{}


Foam::mappedExtrudedWallPolyPatch::mappedExtrudedWallPolyPatch
(
    const mappedExtrudedWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    wallPolyPatch(pp, bm, index, newSize, newStart),
    mappedExtrudedPatchBase(*this, pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedExtrudedWallPolyPatch::~mappedExtrudedWallPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedExtrudedWallPolyPatch::write(Ostream& os) const
{
    wallPolyPatch::write(os);
    mappedExtrudedPatchBase::write(os);
}


// ************************************************************************* //
