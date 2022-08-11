/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "nonConformalErrorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "SubField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalErrorPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, nonConformalErrorPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        nonConformalErrorPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::nonConformalErrorPolyPatch::rename(const wordList& newNames)
{
    nonConformalPolyPatch::rename(newNames);
}


void Foam::nonConformalErrorPolyPatch::reorder(const labelUList& newToOldIndex)
{
    nonConformalPolyPatch::reorder(newToOldIndex);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalErrorPolyPatch::nonConformalErrorPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType),
    nonConformalPolyPatch(static_cast<const polyPatch&>(*this))
{}


Foam::nonConformalErrorPolyPatch::nonConformalErrorPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const word& origPatchName
)
:
    polyPatch(name, size, start, index, bm, patchType),
    nonConformalPolyPatch(*this, origPatchName)
{}


Foam::nonConformalErrorPolyPatch::nonConformalErrorPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    nonConformalPolyPatch(*this, dict)
{}


Foam::nonConformalErrorPolyPatch::nonConformalErrorPolyPatch
(
    const nonConformalErrorPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    nonConformalPolyPatch(*this, pp)
{}


Foam::nonConformalErrorPolyPatch::nonConformalErrorPolyPatch
(
    const nonConformalErrorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& origPatchName
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    nonConformalPolyPatch(*this, origPatchName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalErrorPolyPatch::~nonConformalErrorPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nonConformalErrorPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    nonConformalPolyPatch::write(os);
}


// ************************************************************************* //
