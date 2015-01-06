/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "reflectionVectors.H"
#include "wallFvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::reflectionVectors::reflectionVectors(const Foam::fvMesh& mesh)
:
    n_
    (
        IOobject
        (
            "reflectionVectors",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("n", dimless, vector::zero)
    )
{
    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::reflectionVectors::correct()
{
    const fvMesh& mesh = n_.mesh();
    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        // find the nearest face for every cell
        if (isA<wallFvPatch>(patches[patchi]))
        {
            n_.boundaryField()[patchi] =
                mesh.Sf().boundaryField()[patchi]
               /mesh.magSf().boundaryField()[patchi];
        }
    }
}


// ************************************************************************* //
