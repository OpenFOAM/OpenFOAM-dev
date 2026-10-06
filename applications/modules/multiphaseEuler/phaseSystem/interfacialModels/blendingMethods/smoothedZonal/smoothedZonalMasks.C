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

#include "smoothedZonalMasks.H"
#include "zeroGradientFvPatchField.H"
#include "fvmSup.H"
#include "gaussLaplacianScheme.H"
#include "smoothSolver.H"
#include "symGaussSeidelSmoother.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(smoothedZonal::masks, 0);
}
}


// * * * * * * * * * * * * * * Protected Constructors  * * * * * * * * * * * //

Foam::blendingMethods::smoothedZonal::masks::masks
(
    const word& name,
    const fvMesh& mesh,
    const labelList& zoneIndices,
    const scalar nCellLengthScales,
    const label nIter
)
:
    DemandDrivenMeshObject<fvMesh, MoveableMeshObject, masks>(name, mesh),
    defaultMask_
    (
        IOobject
        (
            name + ":maskDefault",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless, scalar(1)),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    zoneMasks_(zoneIndices.size())
{
    // Set the mask fields to one in their zones and zero elsewhere. Set the
    // default mask to be the opposite; i.e., zero in any specified zone and
    // one elsewhere. Don't worry about zones overlapping. Normalisation will
    // fix this.
    forAll(zoneIndices, zonei)
    {
        const cellZone& zoneCells =
            mesh.cellZones()[zoneIndices[zonei]];

        UIndirectList<scalar>(defaultMask_, zoneCells) = 0;

        zoneMasks_.set
        (
            zonei,
            new volScalarField
            (
                IOobject
                (
                    name + ":mask" + Foam::name(zonei),
                    mesh.time().name(),
                    mesh
                ),
                mesh,
                dimensionedScalar(dimless, scalar(0)),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );

        UIndirectList<scalar>(zoneMasks_[zonei], zoneCells) = 1;
    }

    // Smooth. Do a fixed, small number of iterations of the smooth solver.
    tmp<volInternalScalarField> size = mesh.V().clone();
    const Vector<label>& directions = mesh.geometricD();
    for (direction dir=0; dir<directions.nComponents; dir++)
    {
        if (directions[dir] == -1)
        {
            size.ref() /=
                dimensionedScalar
                (
                    dimensions::length,
                    mesh.bounds().span()[dir]
                );
        }
    }
    const volInternalScalarField invSqrLengthScale
    (
        1/sqr(nCellLengthScales*integerRoot(size, mesh.nGeometricD()))
    );
    const dictionary solveDict = dictionary::entries
    (
        "solver", smoothSolver::typeName,
        "smoother", symGaussSeidelSmoother::typeName,
        "maxIter", nIter,
        "tolerance", -1,
        "relTol", -1
    );
    const int scalarSolverPerformanceDebug0 = SolverPerformance<scalar>::debug;
    SolverPerformance<scalar>::debug = 0;
    auto solve = [&](volScalarField& mask)
    {
        (
            correction(fvm::Sp(invSqrLengthScale, mask))
          - fv::gaussLaplacianScheme<scalar, scalar>::fvmLaplacianUncorrected
            (
                mesh.magSf(),
                mesh.nonOrthDeltaCoeffs(),
                mask
            )
        )->solve(solveDict);
    };
    solve(defaultMask_);
    forAll(zoneIndices, zonei)
    {
        solve(zoneMasks_[zonei]);
    }
    SolverPerformance<scalar>::debug = scalarSolverPerformanceDebug0;

    // Clip. Assume a hyperbolic profile. Invert and then re-evaluate using an
    // approximation of tanh which reaches its -1/+1 limits at a finite value.
    auto clip = [](volScalarField& mask)
    {
        const volScalarField x
        (
            min(max(log(mask/max(1 - mask, vSmall))/2, -3), +3)
        );
        mask = (x*(27 + sqr(x))/(27 + 9*sqr(x)) + 1)/2;
    };
    clip(defaultMask_);
    forAll(zoneIndices, zonei)
    {
        clip(zoneMasks_[zonei]);
    }

    // Normalise
    volScalarField sumMask(defaultMask_);
    forAll(zoneIndices, zonei)
    {
        sumMask += zoneMasks_[zonei];
    }
    forAll(zoneIndices, zonei)
    {
        zoneMasks_[zonei] /= sumMask;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

const Foam::blendingMethods::smoothedZonal::masks&
Foam::blendingMethods::smoothedZonal::masks::New
(
    const fvMesh& mesh,
    const labelList& zoneIndices,
    const scalar nCellLengthScales,
    const label nIter
)
{
    return DemandDrivenMeshObject<fvMesh, MoveableMeshObject, masks>::New
    (
        Foam::name
        (
            Hash<labelList>()(zoneIndices)
          + Hash<scalar>()(nCellLengthScales)
          + Hash<label>()(nIter)
        ),
        mesh,
        zoneIndices,
        nCellLengthScales,
        nIter
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::smoothedZonal::masks::~masks()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::smoothedZonal::masks::movePoints()
{
    return true;
}


// ************************************************************************* //
