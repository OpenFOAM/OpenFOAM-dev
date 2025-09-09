/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "path.H"
#include "edgeMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{
    defineTypeNameAndDebug(path, 0);
    addToRunTimeSelectionTable(extrudeModel, path, dictionary);
}
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::vector Foam::extrudeModels::path::orthogonal(const vector& dirn)
{
    scalar maxOrthogonality = -1;
    direction maxOrthogonalityDir = 0;

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        vector dir(Zero);
        dir[cmpt] = 1;

        const scalar orthogonality = magSqr(dir ^ dirn);

        if (orthogonality > maxOrthogonality)
        {
            maxOrthogonality = orthogonality;
            maxOrthogonalityDir = cmpt;
        }
    }

    vector dir(Zero);
    dir[maxOrthogonalityDir] = 1;

    return dir;
}


Foam::tensor Foam::extrudeModels::path::orthonormalBasis
(
    const vector& T,
    const vector& N
)
{
    const vector B(normalised(T ^ N));
    const vector TNT(normalised(B ^ T));

    return tensor
    (
        T.x(), TNT.x(), B.x(),
        T.y(), TNT.y(), B.y(),
        T.z(), TNT.z(), B.z()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extrudeModels::path::path(const dictionary& dict)
:
    extrudeModel(dict),
    eMeshPtr_(edgeMesh::New(dict.lookup("path"))),
    distances_(eMeshPtr_->points().size(), scalar(0)),
    directions_(distances_.size()),
    normals_(directions_.size()),
    R0T_(Zero)
{
    const pointField& points = eMeshPtr_->points();

    directions_[0] = normalised(points[1] - points[0]);
    normals_[0] = orthogonal(directions_[0]);

    for (label i = 0; i < distances_.size() - 1; i++)
    {
        const vector d(points[i + 1] - points[i]);

        distances_[i + 1] = distances_[i] + mag(d);
        directions_[i + 1] = normalised(d);
        normals_[i + 1] =
            normals_[i] - 2*(directions_[i] & normals_[i])*directions_[i];
    }

    R0T_ = orthonormalBasis(directions_.first(), normals_.first()).T();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extrudeModels::path::~path()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::extrudeModels::path::findFrameIndex
(
    const scalar distance
) const
{
    forAll(distances_, i)
    {
        if (distances_[i] <= distance && distance < distances_[i + 1])
        {
            return i;
        }
    }

    return distances_.size() - 1;
}


// * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * //

Foam::point Foam::extrudeModels::path::operator()
(
    const point& surfacePoint,
    const vector& surfaceNormal,
    const label layer
) const
{
    const pointField& points = eMeshPtr_->points();

    const scalar d = sumThickness(layer)*distances_.last();
    const label frame = findFrameIndex(d);

    vector localCoords(R0T_ & (surfacePoint - points.first()));
    localCoords.x() = 0;

    const tensor R(orthonormalBasis( directions_[frame], normals_[frame]));

    return
    (
        points[frame]
      + directions_[frame]*(d - distances_[frame])
      + (R & localCoords)
    );
}


// ************************************************************************* //
