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

#include "nonConformalFvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalFvPatch, 0);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalFvPatch::nonConformalFvPatch
(
    const fvPatch& patch
)
:
    patch_(patch),
    nonConformalPolyPatch_(refCast<const nonConformalPolyPatch>(patch.patch())),
    faceCells_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalFvPatch::~nonConformalFvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::nonConformalPolyPatch&
Foam::nonConformalFvPatch::nonConformalPatch() const
{
    return nonConformalPolyPatch_;
}


const Foam::word& Foam::nonConformalFvPatch::origPatchName() const
{
    return nonConformalPolyPatch_.origPatchName();
}


Foam::label Foam::nonConformalFvPatch::origPatchID() const
{
    return nonConformalPolyPatch_.origPatchID();
}


const Foam::fvPatch& Foam::nonConformalFvPatch::origPatch() const
{
    return patch_.boundaryMesh()[origPatchID()];
}


const Foam::labelList& Foam::nonConformalFvPatch::polyFaces() const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();

    return
        mesh.conformal()
      ? labelList::null()
      : mesh.polyFacesBf()[patch_.index()];
}


Foam::label Foam::nonConformalFvPatch::start() const
{
    if (size())
    {
        FatalErrorInFunction
            << "The start face is not defined for a " << typeName
            << " patch with a non-zero number of faces"
            << exit(FatalError);
    }

    return patch_.patch().start();
}


Foam::label Foam::nonConformalFvPatch::size() const
{
    return polyFaces().size();
}


const Foam::labelUList& Foam::nonConformalFvPatch::faceCells() const
{
    // !!! This needs an update mechanism, rather than re-calculating the
    // face-cells every time

    const fvMesh& mesh = patch_.boundaryMesh().mesh();

    faceCells_ = UIndirectList<label>(mesh.faceOwner(), polyFaces());

    return faceCells_;
}


// ************************************************************************* //
