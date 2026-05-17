/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2026 OpenFOAM Foundation
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

#include "linearDisplacement_pointMeshMover.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{
    defineTypeNameAndDebug(linearDisplacement, 0);

    addToRunTimeSelectionTable
    (
        pointMeshMover,
        linearDisplacement,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMeshMovers::linearDisplacement::linearDisplacement
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacementPoints0(mesh, dict, typeName),
    axis_(normalised(dict.lookup<vector>("axis", dimless))),
    xFixed_(dict.lookup<scalar>("xFixed", dimLength)),
    xMoving_(dict.lookup<scalar>("xMoving", dimLength)),
    displacement_
    (
        Function1<scalar>::New
        (
            "displacement",
            mesh.time().userUnits(),
            dimLength,
            dict
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMeshMovers::linearDisplacement::~linearDisplacement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::pointMeshMovers::linearDisplacement::newPoints()
{
    tmp<pointField> tcurPoints(new pointField(points0()));
    pointField& curPoints = tcurPoints.ref();

    const scalar displacement = displacement_->value(poly().time().value());

    forAll(curPoints, i)
    {
        const scalar lambda =
            (xFixed_ - (axis_ & curPoints[i]))/(xFixed_ - xMoving_);

        if (lambda > 1)
        {
            curPoints[i] += axis_*displacement;
        }
        else if (lambda > 0)
        {
            curPoints[i] += axis_*lambda*displacement;
        }
    }

    return tcurPoints;
}


void Foam::pointMeshMovers::linearDisplacement::mapMesh(const polyMeshMap& map)
{
    FatalErrorInFunction
        << "Mesh-to-mesh mapping in not implemented for displacement solvers"
        << nl
        << "    velocity based motion solvers are preferable for cases in which"
           " the mesh is reset periodically avoiding accumulation of error."
        << exit(FatalError);
}


// ************************************************************************* //
