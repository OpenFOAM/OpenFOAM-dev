/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "genericPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(genericPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, genericPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, genericPolyPatch, dictionary);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::genericPolyPatch::genericPolyPatch
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


Foam::genericPolyPatch::genericPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    actualTypeName_(dict.lookup("type")),
    dict_(dict)
{}


Foam::genericPolyPatch::genericPolyPatch
(
    const genericPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    actualTypeName_(pp.actualTypeName_),
    dict_(pp.dict_)
{}


Foam::genericPolyPatch::genericPolyPatch
(
    const genericPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    actualTypeName_(pp.actualTypeName_),
    dict_(pp.dict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::genericPolyPatch::~genericPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::genericPolyPatch::write(Ostream& os) const
{
    writeEntry(os, "type", actualTypeName_);
    patchIdentifier::write(os);
    writeEntry(os, "nFaces", size());
    writeEntry(os, "startFace", start());

    forAllConstIter(dictionary, dict_, iter)
    {
        if
        (
            iter().keyword() != "type"
         && iter().keyword() != "nFaces"
         && iter().keyword() != "startFace"
         && iter().keyword() != "inGroups"
         && iter().keyword() != "faces"
        )
        {
            iter().write(os);
        }
    }
}


// ************************************************************************* //
