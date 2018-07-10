/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "patchMeanVelocityForce.H"
#include "processorCyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(patchMeanVelocityForce, 0);

    addToRunTimeSelectionTable
    (
        option,
        patchMeanVelocityForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::patchMeanVelocityForce::patchMeanVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    meanVelocityForce(sourceName, modelType, dict, mesh),
    patch_(coeffs_.lookup("patch")),
    patchi_(mesh.boundaryMesh().findPatchID(patch_))
{
    if (patchi_ < 0)
    {
        FatalErrorInFunction
            << "Cannot find patch " << patch_
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fv::patchMeanVelocityForce::magUbarAve
(
    const volVectorField& U
) const
{
    vector2D sumAmagUsumA
    (
        sum
        (
            (flowDir_ & U.boundaryField()[patchi_])
           *mesh_.boundary()[patchi_].magSf()
        ),
        sum(mesh_.boundary()[patchi_].magSf())
    );


    // If the mean velocity force is applied to a cyclic patch
    // for parallel runs include contributions from processorCyclic patches
    // generated from the decomposition of the cyclic patch
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (Pstream::parRun() && isA<cyclicPolyPatch>(patches[patchi_]))
    {
        labelList processorCyclicPatches
        (
            processorCyclicPolyPatch::patchIDs(patch_, patches)
        );

        forAll(processorCyclicPatches, pcpi)
        {
            const label patchi = processorCyclicPatches[pcpi];

            sumAmagUsumA.x() +=
                sum
                (
                    (flowDir_ & U.boundaryField()[patchi])
                   *mesh_.boundary()[patchi].magSf()
                );

            sumAmagUsumA.y() += sum(mesh_.boundary()[patchi].magSf());
        }
    }

    mesh_.reduce(sumAmagUsumA, sumOp<vector2D>());

    return sumAmagUsumA.x()/sumAmagUsumA.y();
}


// ************************************************************************* //
