/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "searchableSurfacesQueries.H"
#include "projectEdge.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"
#include "pointConstraint.H"
#include "plane.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(projectEdge, 0);
    addToRunTimeSelectionTable(blockEdge, projectEdge, Istream);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::projectEdge::findNearest
(
    const point& pt,
    point& near,
    pointConstraint& constraint
) const
{
    if (surfaces_.size())
    {
        const scalar distSqr = magSqr(points_[end_]-points_[start_]);

        pointField boundaryNear(1);
        List<pointConstraint> boundaryConstraint(1);

        searchableSurfacesQueries::findNearest
        (
            geometry_,
            surfaces_,
            pointField(1, pt),
            scalarField(1, distSqr),
            boundaryNear,
            boundaryConstraint
        );
        near = boundaryNear[0];
        constraint = boundaryConstraint[0];
    }
    else
    {
        near = pt;
        constraint = pointConstraint();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::projectEdge::projectEdge
(
    const searchableSurfaces& geometry,
    const pointField& points,
    Istream& is
)
:
    blockEdge(points, is),
    geometry_(geometry)
{
    wordList names(is);
    surfaces_.setSize(names.size());
    forAll(names, i)
    {
        surfaces_[i] = geometry_.findSurfaceID(names[i]);

        if (surfaces_[i] == -1)
        {
            FatalIOErrorInFunction(is)
                << "Cannot find surface " << names[i] << " in geometry"
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::projectEdge::position(const scalar lambda) const
{
    // Initial guess
    const point start(points_[start_]+lambda*(points_[end_]-points_[start_]));

    point near(start);

    if (lambda >= SMALL && lambda < 1.0-SMALL)
    {
        pointConstraint constraint;
        findNearest(start, near, constraint);
    }

    return near;
}


Foam::tmp<Foam::pointField>
Foam::projectEdge::position(const scalarList& lambdas) const
{
    tmp<pointField> tpoints(new pointField(lambdas.size()));
    pointField& points = tpoints.ref();

    const point& startPt = points_[start_];
    const point& endPt = points_[end_];
    const vector d = endPt-startPt;

    forAll(lambdas, i)
    {
        points[i] = startPt+lambdas[i]*d;
    }

    forAll(lambdas, i)
    {
        if (lambdas[i] >= SMALL && lambdas[i] < 1.0-SMALL)
        {
            pointConstraint constraint;
            findNearest(points[i], points[i], constraint);
        }
    }

    return tpoints;
}


// ************************************************************************* //
