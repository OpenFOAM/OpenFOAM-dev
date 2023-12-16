/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
#include "fvPatchDistWave.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nearWallDist, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearWallDist::nearWallDist(const Foam::fvMesh& mesh)
:
    DemandDrivenMeshObject
    <
        fvMesh,
        DeletableMeshObject,
        nearWallDist
    >(mesh),
    y_
    (
        mesh.boundary(),
        volScalarField::Internal::null(),
        calculatedFvPatchScalarField::typeName
    )
{
    volScalarField yVf(volScalarField::New("y", mesh, dimLength));

    fvPatchDistWave::correct
    (
        mesh,
        mesh.boundaryMesh().findIndices<wallPolyPatch>(),
        -vGreat,
        2,
        yVf
    );

    forAll(y_, patchi)
    {
        const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
        forAll(y_[patchi], patchFacei)
        {
            y_[patchi][patchFacei] = yVf[faceCells[patchFacei]];
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nearWallDist::~nearWallDist()
{}


// ************************************************************************* //
