/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "reconstructionTensors.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fviSurfaceIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvi
{
    defineTypeNameAndDebug(reconstructionTensors, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::fvi::reconstructionTensors::reconstructionTensors(const fvMesh& mesh)
:
    DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        reconstructionTensors
    >(mesh),
    tensors_
    (
        IOobject
        (
            "reconstructionTensors",
            mesh.pointsInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedSymmTensor(inv(dimensions::area), Zero)
    )
{
    calcReconstructTensors();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::fvi::reconstructionTensors::~reconstructionTensors()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvi::reconstructionTensors::calcReconstructTensors()
{
    if (debug)
    {
        InfoInFunction << "Calculating reconstruction tensors" << endl;
    }

    const fvMesh& mesh = this->mesh();

    tensors_ = inv(surfaceSum(sqr(mesh.Sf())/mesh.magSf()), mesh.solutionD());
}


bool Foam::fvi::reconstructionTensors::movePoints()
{
    calcReconstructTensors();
    return true;
}


// ************************************************************************* //
