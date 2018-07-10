/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "displacementMeshMoverMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "localPointRegion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementMeshMoverMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementMeshMoverMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementMeshMoverMotionSolver::displacementMeshMoverMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName) // read pointDisplacement
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementMeshMoverMotionSolver::
~displacementMeshMoverMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::externalDisplacementMeshMover&
Foam::displacementMeshMoverMotionSolver::meshMover() const
{
    if (!meshMoverPtr_.valid())
    {
        const word moverType(coeffDict().lookup("meshMover"));

        meshMoverPtr_ = externalDisplacementMeshMover::New
        (
            moverType,
            coeffDict().optionalSubDict(moverType + "Coeffs"),
            localPointRegion::findDuplicateFacePairs(mesh()),
            pointDisplacement_
        );
    }
    return meshMoverPtr_();
}


Foam::tmp<Foam::pointField>
Foam::displacementMeshMoverMotionSolver::curPoints() const
{
    // Return actual points. Cannot do a reference since complains about
    // assignment to self in polyMesh::movePoints
    return tmp<pointField>(new pointField(mesh().points()));
}


void Foam::displacementMeshMoverMotionSolver::solve()
{
    // The points have moved so before calculation update
    // the mesh and motionSolver accordingly
    movePoints(mesh().points());

    // Update any point motion bcs (e.g. timevarying)
    pointDisplacement().boundaryFieldRef().updateCoeffs();

    label nAllowableErrors = 0;
    labelList checkFaces(identity(mesh().nFaces()));
    meshMover().move
    (
        coeffDict().optionalSubDict(meshMover().type() + "Coeffs"),
        nAllowableErrors,
        checkFaces
    );

    // This will have updated the mesh and implicitly the pointDisplacement
    pointDisplacement().correctBoundaryConditions();
}


void Foam::displacementMeshMoverMotionSolver::movePoints(const pointField& p)
{
    displacementMotionSolver::movePoints(p);

    // Update meshMover for new geometry
    if (meshMoverPtr_.valid())
    {
        meshMover().movePoints(p);
    }
}


void Foam::displacementMeshMoverMotionSolver::updateMesh
(
    const mapPolyMesh& map
)
{
    displacementMotionSolver::updateMesh(map);

    // Update meshMover for new topology
    meshMoverPtr_.clear();
}


// ************************************************************************* //
