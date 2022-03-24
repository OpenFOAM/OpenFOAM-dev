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

#include "nearWallDist.H"
#include "fvMesh.H"
#include "patchDistFuncs.H"
#include "wallFvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearWallDist::nearWallDist(const Foam::fvMesh& mesh)
:
    volScalarField::Boundary
    (
        mesh.boundary(),
        volScalarField::Internal::null(),
        calculatedFvPatchScalarField::typeName
    ),
    mesh_(mesh)
{
    patchDistFuncs::correctBoundaryFaceFaceCells
    (
        mesh_,
        mesh_.boundaryMesh().findPatchIDs<wallPolyPatch>(),
        *this
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nearWallDist::~nearWallDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nearWallDist::correct()
{
    if (mesh_.topoChanging())
    {
        this->setSize(mesh_.boundary().size());

        forAll(*this, patchi)
        {
            this->set
            (
                patchi,
                fvPatchField<scalar>::New
                (
                    calculatedFvPatchScalarField::typeName,
                    mesh_.boundary()[patchi],
                    volScalarField::Internal::null()
                )
            );
        }
    }

    patchDistFuncs::correctBoundaryFaceFaceCells
    (
        mesh_,
        mesh_.boundaryMesh().findPatchIDs<wallPolyPatch>(),
        *this
    );
}


// ************************************************************************* //
