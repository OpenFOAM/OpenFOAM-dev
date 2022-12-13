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
#include "fvPatchDistWave.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nearWallDist, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nearWallDist::resize()
{
    y_.setSize(mesh().boundary().size());

    forAll(y_, patchi)
    {
        y_.set
        (
            patchi,
            fvPatchField<scalar>::New
            (
                calculatedFvPatchScalarField::typeName,
                mesh().boundary()[patchi],
                volScalarField::Internal::null()
            )
        );
    }
}


void Foam::nearWallDist::correct()
{
    volScalarField yVf(volScalarField::New("y", mesh(), dimLength));

    fvPatchDistWave::correct
    (
        mesh(),
        mesh().boundaryMesh().findPatchIDs<wallPolyPatch>(),
        -vGreat,
        2,
        yVf
    );

    forAll(y_, patchi)
    {
        const labelUList& faceCells = mesh().boundary()[patchi].faceCells();
        forAll(y_[patchi], patchFacei)
        {
            y_[patchi][patchFacei] = yVf[faceCells[patchFacei]];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearWallDist::nearWallDist(const Foam::fvMesh& mesh)
:
    DemandDrivenMeshObject
    <
        fvMesh,
        UpdateableMeshObject,
        nearWallDist
    >(mesh),
    y_
    (
        mesh.boundary(),
        volScalarField::Internal::null(),
        calculatedFvPatchScalarField::typeName
    )
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nearWallDist::~nearWallDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::nearWallDist::movePoints()
{
    resize();

    correct();
    return true;
}


void Foam::nearWallDist::topoChange(const polyTopoChangeMap& map)
{
    resize();
    correct();
}


void Foam::nearWallDist::mapMesh(const polyMeshMap& map)
{
    resize();
    correct();
}


void Foam::nearWallDist::distribute(const polyDistributionMap& map)
{
    resize();
    correct();
}


// ************************************************************************* //
