/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "points0MotionSolver.H"
#include "polyDistributionMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(points0MotionSolver, 0);
}


Foam::pointVectorField Foam::points0MotionSolver::readPoints0
(
    const polyMesh& mesh
)
{
    const word instance
    (
        mesh.time().findInstance
        (
            ".",
            "points0",
            IOobject::READ_IF_PRESENT
        )
    );

    if (instance != mesh.time().constant())
    {
        // Points0 written to a time folder

        return pointVectorField
        (
            IOobject
            (
                "points0",
                instance,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(mesh)
        );
    }
    else
    {
        // Return copy of original mesh points

        pointIOField points
        (
            IOobject
            (
                "points",
                mesh.time().constant(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        pointVectorField points0
        (
            IOobject
            (
                "points0",
                instance,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(mesh),
            dimensionedVector(dimLength, Zero)
        );

        points0.primitiveFieldRef() = points;

        return points0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::points0MotionSolver::points0MotionSolver
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict,
    const word& type
)
:
    motionSolver(name, mesh, type),
    points0_(readPoints0(mesh))
{
    if (points0_.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of points in mesh " << mesh.nPoints()
            << " differs from number of points " << points0_.size()
            << " read from file "
            << typeIOobject<pointIOField>
               (
                   "points",
                   mesh.time().constant(),
                   polyMesh::meshSubDir,
                   mesh,
                   IOobject::MUST_READ,
                   IOobject::NO_WRITE,
                   false
               ).filePath()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::points0MotionSolver::~points0MotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::points0MotionSolver::movePoints(const pointField&)
{}


void Foam::points0MotionSolver::topoChange(const polyTopoChangeMap& map)
{
    NotImplemented;
}


void Foam::points0MotionSolver::mapMesh(const polyMeshMap& map)
{
    points0_.primitiveFieldRef() = mesh().points();

    // The processor boundaries may have changed, so we need to update the
    // boundary field. There is no data in this field, so we don't need to map
    // anything. We can just reset it to a freshly created calculated field.
    points0_.boundaryFieldRef().reset
    (
        pointVectorField::Boundary
        (
            points0_.mesh().boundary(),
            points0_.internalField(),
            calculatedPointPatchVectorField::typeName
        )
    );
}


void Foam::points0MotionSolver::distribute
(
    const polyDistributionMap& map
)
{
    map.distributePointData(points0_.primitiveFieldRef());

    // See above
    points0_.boundaryFieldRef().reset
    (
        pointVectorField::Boundary
        (
            points0_.mesh().boundary(),
            points0_.internalField(),
            calculatedPointPatchVectorField::typeName
        )
    );
}


bool Foam::points0MotionSolver::write() const
{
    if (points0_.writeOpt() == IOobject::AUTO_WRITE)
    {
        points0_.write();
    }

    return motionSolver::write();
}


// ************************************************************************* //
