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

#include "nonConformalMappedWallFvPatch.H"
#include "nonConformalErrorFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "nonConformalMappedPolyFacesFvsPatchLabelField.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalMappedWallFvPatch, 0);
    addToRunTimeSelectionTable
    (
        fvPatch,
        nonConformalMappedWallFvPatch,
        polyPatch
    );
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalMappedWallFvPatch::nonConformalMappedWallFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    wallFvPatch(patch, bm),
    nonConformalFvPatch(static_cast<const fvPatch&>(*this)),
    nonConformalMappedFvPatchBase(static_cast<const fvPatch&>(*this)),
    nonConformalMappedWallPolyPatch_
    (
        refCast<const nonConformalMappedWallPolyPatch>(patch)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalMappedWallFvPatch::~nonConformalMappedWallFvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::nonConformalMappedWallPolyPatch&
Foam::nonConformalMappedWallFvPatch::nonConformalMappedWallPatch() const
{
    return nonConformalMappedWallPolyPatch_;
}


const Foam::nonConformalMappedWallFvPatch&
Foam::nonConformalMappedWallFvPatch::nbrPatch() const
{
    return refCast<const nonConformalMappedWallFvPatch>(nbrFvPatch());
}


bool Foam::nonConformalMappedWallFvPatch::owner() const
{
    return nonConformalMappedWallPolyPatch_.owner();
}


Foam::label Foam::nonConformalMappedWallFvPatch::start() const
{
    return nonConformalFvPatch::start();
}


Foam::label Foam::nonConformalMappedWallFvPatch::size() const
{
    return nonConformalFvPatch::size();
}


const Foam::labelUList& Foam::nonConformalMappedWallFvPatch::faceCells() const
{
    return nonConformalFvPatch::faceCells();
}


Foam::word Foam::nonConformalMappedWallFvPatch::polyFacesType() const
{
    return nonConformalMappedPolyFacesFvsPatchLabelField::typeName;
}


// ************************************************************************* //
