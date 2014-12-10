/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "mappedVariableThicknessWallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedVariableThicknessWallPolyPatch, 0);

    addToRunTimeSelectionTable
    (
        polyPatch,
        mappedVariableThicknessWallPolyPatch,
        word
    );

    addToRunTimeSelectionTable
    (
        polyPatch,
        mappedVariableThicknessWallPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedVariableThicknessWallPolyPatch::mappedVariableThicknessWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    mappedWallPolyPatch(name, size, start, index, bm, patchType),
    thickness_(size)
{}


Foam::mappedVariableThicknessWallPolyPatch::mappedVariableThicknessWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const mappedPatchBase::sampleMode mode,
    const word& samplePatch,
    const vectorField& offset,
    const polyBoundaryMesh& bm
)
:
    mappedWallPolyPatch(name, size, start, index, bm, typeName),
    thickness_(size)
{}


Foam::mappedVariableThicknessWallPolyPatch::mappedVariableThicknessWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const mappedPatchBase::sampleMode mode,
    const word& samplePatch,
    const vector& offset,
    const polyBoundaryMesh& bm
)
:
    mappedWallPolyPatch(name, size, start, index, bm, typeName),
    thickness_(size)
{}


Foam::mappedVariableThicknessWallPolyPatch::mappedVariableThicknessWallPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    mappedWallPolyPatch(name, dict, index, bm, patchType),
    thickness_(scalarField("thickness", dict, this->size()))
{}


Foam::mappedVariableThicknessWallPolyPatch::
mappedVariableThicknessWallPolyPatch
(
    const mappedVariableThicknessWallPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    mappedWallPolyPatch(pp, bm),
    thickness_(pp.thickness_)
{}


Foam::mappedVariableThicknessWallPolyPatch::mappedVariableThicknessWallPolyPatch
(
    const mappedVariableThicknessWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    mappedWallPolyPatch(pp, bm, index, newSize, newStart),
    thickness_(newSize)
{}


Foam::mappedVariableThicknessWallPolyPatch::mappedVariableThicknessWallPolyPatch
(
    const mappedVariableThicknessWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    mappedWallPolyPatch(pp, bm, index, mapAddressing, newStart),
    thickness_(pp.size())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedVariableThicknessWallPolyPatch::
~mappedVariableThicknessWallPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedVariableThicknessWallPolyPatch::
write(Foam::Ostream& os) const
{
    os.writeKeyword("thickness") << thickness_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
