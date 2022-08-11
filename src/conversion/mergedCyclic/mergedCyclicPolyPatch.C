/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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

#include "mergedCyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mergedCyclicPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, mergedCyclicPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, mergedCyclicPolyPatch, dictionary);
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mergedCyclicPolyPatch::mergedCyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType)
{}


Foam::mergedCyclicPolyPatch::mergedCyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType)
{}


Foam::mergedCyclicPolyPatch::mergedCyclicPolyPatch
(
    const mergedCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm)
{}


Foam::mergedCyclicPolyPatch::mergedCyclicPolyPatch
(
    const mergedCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart)
{}


// ************************************************************************* //
