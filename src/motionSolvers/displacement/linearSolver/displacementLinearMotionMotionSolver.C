/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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

#include "displacementLinearMotionMotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementLinearMotionMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementLinearMotionMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementLinearMotionMotionSolver::
displacementLinearMotionMotionSolver
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    points0MotionSolver(name, mesh, dict, typeName),
    axis_(normalised(vector(coeffDict().lookup("axis")))),
    xFixed_(coeffDict().lookup<scalar>("xFixed")),
    xMoving_(coeffDict().lookup<scalar>("xMoving")),
    displacement_(Function1<scalar>::New("displacement", coeffDict()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementLinearMotionMotionSolver::
~displacementLinearMotionMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::displacementLinearMotionMotionSolver::curPoints() const
{
    tmp<pointField> tcurPoints(new pointField(points0()));
    pointField& curPoints = tcurPoints.ref();

    const scalar t = mesh().time().userTimeValue();

    const scalar displacement = displacement_->value(t);

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


// ************************************************************************* //
