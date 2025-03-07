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

#include "Spline.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Derived>
Foam::Spline<Derived>::Spline
(
    const pointField& knots,
    const scalar tol,
    const label nIter
)
:
    knots_(knots),
    tol_(tol),
    nIter_(nIter)
{}


template<class Derived>
Foam::Spline<Derived>::Spline(const pointField& knots)
:
    knots_(knots),
    tol_(rootSmall),
    nIter_(16)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Derived>
Foam::point Foam::Spline<Derived>::position(const scalar lambda) const
{
    const Derived& derived = static_cast<const Derived&>(*this);

    // Quick return if the parameter is out of range
    if (lambda <= 0) return start();
    if (lambda >= 1) return end();

    // Evaluate the positions for the given local parameters
    auto evaluatePoint = [&derived](const scalar& lambda)
    {
        const scalar l = lambda*derived.nSegments();
        const label segmenti = floor(l);
        const scalar segmentLambda = l - segmenti;
        return derived.position(segmenti, segmentLambda);
    };

    // Algorithm for iteratively improving adherence to the desired grading
    scalar lambdaStar = lambda;
    scalar error = vGreat;
    for (label iter = 0; error > tol_ && iter < nIter_; ++ iter)
    {
        // Evaluate the position at the current local parameter
        const point p = evaluatePoint(lambdaStar);

        // Calculate the actual local parameter based on cumulative distance
        const scalar dist0 = mag(p - start());
        const scalar dist1 = mag(p - end());
        const scalar lambdaDist = dist0/(dist0 + dist1);

        // Update the error
        error = mag(lambda - lambdaDist);

        // Interpolate a new parameter to try and space the cumulative distance
        // as required by the original parameter
        const scalar lambdaStar0 = lambda < lambdaDist ? 0 : lambdaStar;
        const scalar lambdaStar1 = lambda < lambdaDist ? lambdaStar : 1;
        const scalar lambdaDist0 = lambda < lambdaDist ? 0 : lambdaDist;
        const scalar lambdaDist1 = lambda < lambdaDist ? lambdaDist : 1;
        const scalar x =
            (lambda - lambdaDist0)/max(lambdaDist1 - lambdaDist0, vSmall);
        lambdaStar = (1 - x)*lambdaStar0 + x*lambdaStar1;
    }

    // Evaluate the position at the final local parameter and return
    return evaluatePoint(lambdaStar);
}


template<class Derived>
Foam::tmp<Foam::pointField> Foam::Spline<Derived>::position
(
    const scalarList& lambdasArg
) const
{
    const Derived& derived = static_cast<const Derived&>(*this);

    const label n = lambdasArg.size();

    // Clip the parameters to the bounds
    const scalarField lambdas
    (
        min(max(scalarField(lambdasArg), scalar(0)), scalar(1))
    );

    // Workspace
    tmp<scalarField> tlambdasStar(new scalarField(lambdas));
    tmp<scalarField> tlambdasStarNew(new scalarField(n));
    scalarField lambdasDist(n);

    // Allocate the result
    tmp<pointField> tpoints(new pointField(n));
    pointField& points = tpoints.ref();

    // Evaluate the positions for the given local parameters
    auto evaluatePoints = [&derived,&points](const scalarField& lambdas)
    {
        forAll(lambdas, i)
        {
            const scalar l = lambdas[i]*derived.nSegments();
            const label segmenti = floor(l);
            const scalar segmentLambda = l - segmenti;
            points[i] = derived.position(segmenti, segmentLambda);
        }
    };

    // Algorithm for iteratively improving adherence to the desired grading
    scalar error = vGreat;
    for (label iter = 0; error > tol_ && iter < nIter_; ++ iter)
    {
        scalarField& lambdasStar = tlambdasStar.ref();
        scalarField& lambdasStarNew = tlambdasStarNew.ref();

        // Evaluate the positions at the current local parameters
        evaluatePoints(lambdasStar);

        // Calculate the actual local parameters based on cumulative distance
        scalarField lambdasDist(n);
        lambdasDist.first() = mag(points.first() - start());
        for (label i = 1; i < n; ++ i)
        {
            lambdasDist[i] =
                lambdasDist[i - 1] + mag(points[i] - points[i - 1]);
        }
        lambdasDist /= lambdasDist.last() + mag(end() - points.last());

        // Update the error
        error = -vGreat;
        forAll(lambdas, i)
        {
            error = max(error, mag(lambdas[i] - lambdasDist[i]));
        }

        // Interpolate new parameters to try and space the cumulative distance
        // as required by the original parameters
        label iNew = 0;
        for (label i = 0; i <= n; ++ i)
        {
            const scalar lambdaStar0 = i == 0 ? 0 : lambdasStar[i - 1];
            const scalar lambdaStar1 = i < n ? lambdasStar[i] : 1;

            const scalar lambdaDist0 = i == 0 ? 0 : lambdasDist[i - 1];
            const scalar lambdaDist1 = i < n ? lambdasDist[i] : 1;

            while (iNew < n && lambdas[iNew] <= lambdaDist1)
            {
                const scalar x =
                    (lambdas[iNew] - lambdaDist0)
                   /max(lambdaDist1 - lambdaDist0, vSmall);

                lambdasStarNew[iNew] = (1 - x)*lambdaStar0 + x*lambdaStar1;

                iNew ++;
            }
        }

        // Set the old parameters to the new ones and iterate
        Swap(tlambdasStar, tlambdasStarNew);
    }

    // Evaluate the positions at the final local parameters
    evaluatePoints(tlambdasStar());

    return tpoints;
}


// ************************************************************************* //
