/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "nearWallDistNoSearch.H"
#include "fvMesh.H"
#include "wallFvPatch.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nearWallDistNoSearch::doAll()
{
    const volVectorField& cellCentres = mesh_.C();
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        fvPatchScalarField& ypatch = operator[](patchi);

        if (isA<wallFvPatch>(patches[patchi]))
        {
            const labelUList& faceCells = patches[patchi].faceCells();

            const fvPatchVectorField& patchCentres
                = cellCentres.boundaryField()[patchi];

            const fvsPatchVectorField& Apatch
                = mesh_.Sf().boundaryField()[patchi];

            const fvsPatchScalarField& magApatch
                = mesh_.magSf().boundaryField()[patchi];

            forAll(patchCentres, facei)
            {
                ypatch[facei] =
                (
                    Apatch[facei] &
                    (
                        patchCentres[facei]
                      - cellCentres[faceCells[facei]]
                    )
                )/magApatch[facei];
            }
        }
        else
        {
            ypatch = 0.0;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearWallDistNoSearch::nearWallDistNoSearch(const Foam::fvMesh& mesh)
:
    volScalarField::Boundary
    (
        mesh.boundary(),
        mesh.V(),           // Dummy internal field
        calculatedFvPatchScalarField::typeName
    ),
    mesh_(mesh)
{
    doAll();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nearWallDistNoSearch::~nearWallDistNoSearch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nearWallDistNoSearch::correct()
{
    if (mesh_.changing())
    {
        // Update size of Boundary
        forAll(mesh_.boundary(), patchi)
        {
            operator[](patchi).setSize(mesh_.boundary()[patchi].size());
        }
    }

    doAll();
}


// ************************************************************************* //
