/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "fvMeshToFvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshToFvMesh, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshToFvMesh::fvMeshToFvMesh
(
    const fvMesh& srcMesh,
    const fvMesh& tgtMesh,
    const word& engineType,
    const HashTable<word>& patchMap
)
:
    meshToMesh(srcMesh, tgtMesh, engineType, patchMap),
    srcMesh_(srcMesh),
    tgtMesh_(tgtMesh)
{
    if (debug)
    {
        Info<< typeName << ": Writing target coverage" << endl;

        volScalarField::Internal
        (
            "tgtCoverage",
            srcToTgt<scalar>
            (
                volScalarField::Internal::New
                (
                    "1",
                    srcMesh,
                    dimensionedScalar(dimless, scalar(1))
                )(),
                volScalarField::Internal::New
                (
                    "0",
                    srcMesh,
                    dimensionedScalar(dimless, scalar(0))
                )()
            )
        ).write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshToFvMesh::~fvMeshToFvMesh()
{}


// ************************************************************************* //
