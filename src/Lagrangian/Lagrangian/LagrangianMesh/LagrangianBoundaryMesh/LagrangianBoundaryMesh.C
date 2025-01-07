/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianBoundaryMesh.H"
#include "LagrangianPatch.H"
#include "LagrangianMesh.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LagrangianBoundaryMesh::LagrangianBoundaryMesh
(
    const LagrangianMesh& mesh,
    const polyBoundaryMesh& pbm
)
:
    LagrangianBoundaryMesh(mesh, pbm, pbm.types())
{}


Foam::LagrangianBoundaryMesh::LagrangianBoundaryMesh
(
    const LagrangianMesh& mesh,
    const polyBoundaryMesh& pbm,
    const wordList& wantedPatchTypes
)
:
    PtrList<LagrangianPatch>(pbm.size()),
    mesh_(mesh)
{
    forAll(*this, patchi)
    {
        this->set
        (
            patchi,
            LagrangianPatch::New
            (
                wantedPatchTypes[patchi],
                pbm[patchi],
                *this
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LagrangianBoundaryMesh::~LagrangianBoundaryMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::LagrangianBoundaryMesh::findIndex
(
    const word& patchName
) const
{
    return mesh().mesh().boundaryMesh().findIndex(patchName);
}


Foam::labelList Foam::LagrangianBoundaryMesh::findIndices
(
    const wordRe& key,
    const bool useGroups
) const
{
    return mesh().mesh().boundaryMesh().findIndices(key, useGroups);
}


// ************************************************************************* //
