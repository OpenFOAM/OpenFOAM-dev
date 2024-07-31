/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2024 OpenFOAM Foundation
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
#include "volFields.H"
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
        fvConstraint,
        patchMeanVelocityForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::patchMeanVelocityForce::readCoeffs(const dictionary& dict)
{
    patch_ = dict.lookup<word>("patch");

    if (mesh().boundaryMesh().findIndex(patch_) < 0)
    {
        FatalErrorInFunction
            << "Cannot find patch " << patch_
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::patchMeanVelocityForce::patchMeanVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    meanVelocityForce(sourceName, modelType, mesh, dict),
    patch_(word::null)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fv::patchMeanVelocityForce::magUbarAve
(
    const volVectorField& U
) const
{
    const label patchi = mesh().boundaryMesh().findIndex(patch_);

    scalar sumA = sum(mesh().boundary()[patchi].magSf());
    scalar sumAmagU =
        sum
        (
            mesh().boundary()[patchi].magSf()
           *(normalised(Ubar()) & U.boundaryField()[patchi])
        );

    // If the mean velocity force is applied to a cyclic patch
    // for parallel runs include contributions from processorCyclic patches
    // generated from the decomposition of the cyclic patch
    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    if (Pstream::parRun() && isA<cyclicPolyPatch>(patches[patchi]))
    {
        const labelList processorCyclicPatches
        (
            processorCyclicPolyPatch::patchIDs(patch_, patches)
        );

        forAll(processorCyclicPatches, pcpi)
        {
            const label patchi = processorCyclicPatches[pcpi];

            sumA += sum(mesh().boundary()[patchi].magSf());
            sumAmagU +=
                sum
                (
                    mesh().boundary()[patchi].magSf()
                   *(normalised(Ubar()) & U.boundaryField()[patchi])
                );
        }
    }

    mesh().reduce(sumA, sumOp<scalar>());
    mesh().reduce(sumAmagU, sumOp<scalar>());

    return sumAmagU/sumA;
}


bool Foam::fv::patchMeanVelocityForce::read(const dictionary& dict)
{
    if (meanVelocityForce::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
