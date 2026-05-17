/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "displacementSBRStress_fvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "fvcLaplacian.H"
#include "polyTopoChangeMap.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMotionSolvers
{
    defineTypeNameAndDebug(displacementSBRStress, 0);

    addToRunTimeSelectionTable
    (
        fvMeshMover,
        displacementSBRStress,
        fvMesh
    );

    addToRunTimeSelectionTable
    (
        pointMeshMover,
        displacementSBRStress,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMotionSolvers::displacementSBRStress::displacementSBRStress
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    fvMotionSolver(mesh),
    pointMeshMovers::displacement(mesh, dict, typeName),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMotionSolver::mesh(),
        dimensionedVector
        (
            "cellDisplacement",
            pointDisplacement().dimensions(),
            Zero
        ),
        cellMotionBoundaryTypes<vector>(pointDisplacement().boundaryField())
    ),
    diffusivityType_(dict.lookup("diffusivity")),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMotionSolver::mesh(), diffusivityType_)
    )
{}


Foam::fvMotionSolvers::displacementSBRStress::displacementSBRStress
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    displacementSBRStress(mesh.poly(), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMotionSolvers::displacementSBRStress::~displacementSBRStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::fvMotionSolvers::displacementSBRStress::newPoints()
{
    // The points have moved so before interpolation update
    // the mtionSolver accordingly
    movePoints(mesh().points());

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

    volPointInterpolation::New(mesh()).interpolate
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


void Foam::fvMotionSolvers::displacementSBRStress::topoChange
(
    const polyTopoChangeMap& map
)
{
    pointMeshMovers::displacement::topoChange(map);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityType_.rewind();
    diffusivityPtr_ = motionDiffusivity::New
    (
        mesh(),
        diffusivityType_
    );
}


void Foam::fvMotionSolvers::displacementSBRStress::mapMesh
(
    const polyMeshMap& map
)
{
    pointMeshMovers::displacement::mapMesh(map);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityType_.rewind();
    diffusivityPtr_ = motionDiffusivity::New
    (
        mesh(),
        diffusivityType_
    );
}


// ************************************************************************* //
