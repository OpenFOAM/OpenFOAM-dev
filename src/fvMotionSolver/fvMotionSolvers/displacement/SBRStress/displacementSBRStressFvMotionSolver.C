/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "displacementSBRStressFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "fvcLaplacian.H"
#include "mapPolyMesh.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementSBRStressFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementSBRStressFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementSBRStressFvMotionSolver::displacementSBRStressFvMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector
        (
            "cellDisplacement",
            pointDisplacement().dimensions(),
            Zero
        ),
        cellMotionBoundaryTypes<vector>(pointDisplacement().boundaryField())
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementSBRStressFvMotionSolver::
~displacementSBRStressFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::displacementSBRStressFvMotionSolver::curPoints() const
{
    volPointInterpolation::New(fvMesh_).interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    tmp<pointField> tcurPoints
    (
        points0() + pointDisplacement().primitiveField()
    );

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::displacementSBRStressFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the mtionSolver accordingly
    movePoints(fvMesh_.points());

    diffusivityPtr_->correct();
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    surfaceScalarField Df(diffusivityPtr_->operator()());

    volTensorField gradCd("gradCd", fvc::grad(cellDisplacement_));

    Foam::solve
    (
        fvm::laplacian
        (
            2*Df,
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )

      + fvc::div
        (
            Df
           *(
               fvc::dotInterpolate
               (
                   cellDisplacement_.mesh().Sf(),
                   gradCd.T() - gradCd
               )

               // Solid-body rotation "lambda" term
             - cellDisplacement_.mesh().Sf()*fvc::interpolate(tr(gradCd))
            )
        )

        /*
      - fvc::laplacian
        (
            2*Df,
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )

      + fvc::div
        (
            Df
           *(
               fvc::dotInterpolate
               (
                   cellDisplacement_.mesh().Sf(),
                   gradCd + gradCd.T()
               )

               // Solid-body rotation "lambda" term
             - cellDisplacement_.mesh().Sf()*fvc::interpolate(tr(gradCd))
           )
        )
        */
    );
}


void Foam::displacementSBRStressFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityPtr_ = motionDiffusivity::New
    (
        fvMesh_,
        coeffDict().lookup("diffusivity")
    );
}


// ************************************************************************* //
