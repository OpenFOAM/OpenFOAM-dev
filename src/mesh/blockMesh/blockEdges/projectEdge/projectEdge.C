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
#include "OBJstream.H"

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
    static label iter = 0;

    tmp<pointField> tpoints(new pointField(lambdas.size()));
    pointField& points = tpoints.ref();

    const point& startPt = points_[start_];
    const point& endPt = points_[end_];
    const vector d = endPt-startPt;

    // Initial guess
    forAll(lambdas, i)
    {
        points[i] = startPt+lambdas[i]*d;
    }


    for (label i = 0; i < 3; i++)
    {
        // Do projection
        {
            List<pointConstraint> constraints(lambdas.size());
            searchableSurfacesQueries::findNearest
            (
                geometry_,
                surfaces_,
                pointField(points),
                scalarField(points.size(), magSqr(d)),
                points,
                constraints
            );

            // Reset start and end point
            if (lambdas[0] < SMALL)
            {
                points[0] = startPt;
            }
            if (lambdas.last() > 1.0-SMALL)
            {
                points.last() = endPt;
            }
        }

        // Calculate distances
        scalarField nearLength(points.size());
        {
            nearLength[0] = 0.0;
            for(label i = 1; i < points.size(); i++)
            {
                nearLength[i] = nearLength[i-1] + mag(points[i]-points[i-1]);
            }
        }

        // Compare actual distances and move points (along straight line;
        // not along surface)
        for(label i = 1; i < points.size() - 1; i++)
        {
            scalar nearDelta = mag(points[i]-points[i-1])/nearLength.last();
            scalar wantedDelta = lambdas[i]-lambdas[i-1];

            vector v(points[i]-points[i-1]);
            points[i] = points[i-1]+wantedDelta/nearDelta*v;
        }
    }


    if (debug)
    {
        OBJstream str("projectEdge_" + Foam::name(iter++) + ".obj");
        Info<< "Writing lines from straight-line start points"
            << " to projected points to " << str.name() << endl;

        pointField startPts(lambdas.size());
        forAll(lambdas, i)
        {
            startPts[i] = startPt+lambdas[i]*d;
        }

        pointField nearPts(lambdas.size());
        List<pointConstraint> nearConstraints(lambdas.size());
        {
            const scalar distSqr = magSqr(d);
            searchableSurfacesQueries::findNearest
            (
                geometry_,
                surfaces_,
                startPts,
                scalarField(startPts.size(), distSqr),
                nearPts,
                nearConstraints
            );
        }

        forAll(startPts, i)
        {
            str.write(linePointRef(startPts[i], nearPts[i]));
            str.write(linePointRef(nearPts[i], points[i]));
        }
    }

    return tpoints;
}


// ************************************************************************* //
