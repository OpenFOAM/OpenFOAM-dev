/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshToFvMesh::writeTargetCoverage() const
{
    Info<< typeName << ": Writing target coverage" << endl;

    volInternalScalarField
    (
        "tgtCoverage",
        srcToTgt<scalar>
        (
            volInternalScalarField::New
            (
                "1",
                srcMesh_,
                dimensionedScalar(dimless, scalar(1))
            )(),
            volInternalScalarField::New
            (
                "0",
                tgtMesh_,
                dimensionedScalar(dimless, scalar(0))
            )()
        )
    ).write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshToFvMesh::fvMeshToFvMesh
(
    const fvMesh& srcMesh,
    const fvMesh& tgtMesh,
    const word& cellsInterpolationType,
    const word& cellsExtrapolationType,
    const word& patchInterpolationType,
    const word& patchExtrapolationType,
    const HashTable<word>& patchMap
)
:
    meshToMesh
    (
        srcMesh,
        tgtMesh,
        cellsInterpolationType,
        cellsExtrapolationType,
        patchInterpolationType,
        patchExtrapolationType,
        patchMap
    ),
    srcMesh_(srcMesh),
    tgtMesh_(tgtMesh)
{
    if (debug) writeTargetCoverage();
}


Foam::fvMeshToFvMesh::fvMeshToFvMesh
(
    const fvMesh& srcMesh,
    const fvMesh& tgtMesh,
    const word& cellsInterpolationType,
    const word& patchInterpolationType,
    const HashTable<word>& patchMap
)
:
    meshToMesh
    (
        srcMesh,
        tgtMesh,
        cellsInterpolationType,
        patchInterpolationType,
        patchMap
    ),
    srcMesh_(srcMesh),
    tgtMesh_(tgtMesh)
{
    if (debug) writeTargetCoverage();
}


Foam::fvMeshToFvMesh::fvMeshToFvMesh
(
    const fvMesh& srcMesh,
    const fvMesh& tgtMesh,
    const word& interpolationType,
    const HashTable<word>& patchMap
)
:
    meshToMesh(srcMesh, tgtMesh, interpolationType, patchMap),
    srcMesh_(srcMesh),
    tgtMesh_(tgtMesh)
{
    if (debug) writeTargetCoverage();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshToFvMesh::~fvMeshToFvMesh()
{}


// ************************************************************************* //
