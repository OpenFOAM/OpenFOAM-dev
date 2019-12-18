/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "dynamicInterpolatedFvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicInterpolatedFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicInterpolatedFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicInterpolatedFvMesh::dynamicInterpolatedFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_(dynamicMeshDict().optionalSubDict(typeName + "Coeffs")),
    pointInterpolator_(*this, dynamicMeshCoeffs_),
    displacement_(dynamicMeshCoeffs_.lookup("displacement")),
    points0_
    (
        displacement_
      ? new pointIOField(points0IO(*this))
      : nullptr
    ),
    velocityMotionCorrection_(*this, dynamicMeshDict())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicInterpolatedFvMesh::~dynamicInterpolatedFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicInterpolatedFvMesh::update()
{
    if (displacement_)
    {
        fvMesh::movePoints(points0_() + pointInterpolator_.curPointField()());
    }
    else
    {
        fvMesh::movePoints(pointInterpolator_.curPointField());
    }

    velocityMotionCorrection_.update();

    return true;
}


// ************************************************************************* //
