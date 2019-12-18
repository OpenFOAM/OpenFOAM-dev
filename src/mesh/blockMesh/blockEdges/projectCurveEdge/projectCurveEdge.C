/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
#include "projectCurveEdge.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"
#include "pointConstraint.H"
#include "OBJstream.H"
#include "linearInterpolationWeights.H"
#include "searchableExtrudedCircle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(projectCurveEdge, 0);
    addToRunTimeSelectionTable(blockEdge, projectCurveEdge, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::projectCurveEdge::projectCurveEdge
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces& geometry,
    const pointField& points,
    Istream& is
)
:
    blockEdge(dict, index, points, is),
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

        if (isA<searchableExtrudedCircle>(geometry_[surfaces_[i]]))
        {
            Info<< type() << " : Using curved surface "
                << geometry_[surfaces_[i]].name()
                << " to predict starting points." << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::projectCurveEdge::position(const scalarList& lambdas) const
{
    // For debugging to tag the output
    static label eIter = 0;

    autoPtr<OBJstream> debugStr;
    if (debug)
    {
        debugStr.reset
        (
            new OBJstream("projectCurveEdge_" + Foam::name(eIter++) + ".obj")
        );
        Info<< "Writing lines from straight-line start points"
            << " to projected points to " << debugStr().name() << endl;
    }


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

    // Use special interpolation to keep initial guess on same position on
    // surface
    forAll(surfaces_, i)
    {
        if (isA<searchableExtrudedCircle>(geometry_[surfaces_[i]]))
        {
            const searchableExtrudedCircle& s =
            refCast<const searchableExtrudedCircle>
            (
                geometry_[surfaces_[i]]
            );
            List<pointIndexHit> nearInfo;
            s.findParametricNearest
            (
                points[0],
                points.last(),
                scalarField(lambdas),
                scalarField(points.size(), magSqr(d)),
                nearInfo
            );
            forAll(nearInfo, i)
            {
                if (nearInfo[i].hit())
                {
                    points[i] = nearInfo[i].hitPoint();
                }
            }

            break;
        }
    }



    // Upper limit for number of iterations
    const label maxIter = 10;

    // Residual tolerance
    const scalar relTol = 0.1;
    const scalar absTol = 1e-4;

    scalar initialResidual = 0;

    for (label iter = 0; iter < maxIter; iter++)
    {
        // Do projection
        {
            List<pointConstraint> constraints(lambdas.size());
            pointField start(points);
            searchableSurfacesQueries::findNearest
            (
                geometry_,
                surfaces_,
                start,
                scalarField(start.size(), magSqr(d)),
                points,
                constraints
            );

            // Reset start and end point
            if (lambdas[0] < small)
            {
                points[0] = startPt;
            }
            if (lambdas.last() > 1 - small)
            {
                points.last() = endPt;
            }

            if (debugStr.valid())
            {
                forAll(points, i)
                {
                    debugStr().write(linePointRef(start[i], points[i]));
                }
            }
        }

        // Calculate lambdas (normalised coordinate along edge)
        scalarField projLambdas(points.size());
        {
            projLambdas[0] = 0;
            for (label i = 1; i < points.size(); i++)
            {
                projLambdas[i] =
                    projLambdas[i-1] + mag(points[i] - points[i-1]);
            }
            projLambdas /= projLambdas.last();
        }
        linearInterpolationWeights interpolator(projLambdas);

        // Compare actual distances and move points (along straight line;
        // not along surface)
        vectorField residual(points.size(), vector::zero);
        labelList indices;
        scalarField weights;
        for (label i = 1; i < points.size() - 1; i++)
        {
            interpolator.valueWeights(lambdas[i], indices, weights);

            point predicted = vector::zero;
            forAll(indices, indexi)
            {
                predicted += weights[indexi]*points[indices[indexi]];
            }
            residual[i] = predicted - points[i];
        }

        const scalar scalarResidual = sum(mag(residual));

        if (debug)
        {
            Pout<< "Iter:" << iter << " initialResidual:" << initialResidual
                << " residual:" << scalarResidual << endl;
        }

        if (scalarResidual < absTol*0.5*lambdas.size())
        {
            break;
        }
        else if (iter == 0)
        {
            initialResidual = scalarResidual;
        }
        else if (scalarResidual/initialResidual < relTol)
        {
            break;
        }


        if (debugStr.valid())
        {
            forAll(points, i)
            {
                const point predicted(points[i] + residual[i]);
                debugStr().write(linePointRef(points[i], predicted));
            }
        }

        points += residual;
    }

    return tpoints;
}


// ************************************************************************* //
