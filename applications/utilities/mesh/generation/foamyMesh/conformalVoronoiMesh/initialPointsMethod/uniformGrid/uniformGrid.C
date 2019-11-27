/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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

#include "uniformGrid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(uniformGrid, 0);
addToRunTimeSelectionTable(initialPointsMethod, uniformGrid, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uniformGrid::uniformGrid
(
    const dictionary& initialPointsDict,
    const Time& runTime,
    Random& rndGen,
    const conformationSurfaces& geometryToConformTo,
    const cellShapeControl& cellShapeControls,
    const autoPtr<backgroundMeshDecomposition>& decomposition
)
:
    initialPointsMethod
    (
        typeName,
        initialPointsDict,
        runTime,
        rndGen,
        geometryToConformTo,
        cellShapeControls,
        decomposition
    ),
    initialCellSize_(detailsDict().lookup<scalar>("initialCellSize")),
    randomiseInitialGrid_(detailsDict().lookup("randomiseInitialGrid")),
    randomPerturbationCoeff_
    (
        detailsDict().lookup<scalar>("randomPerturbationCoeff")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

List<Vb::Point> uniformGrid::initialPoints() const
{
    boundBox bb;

    // Pick up the bounds of this processor, or the whole geometry, depending
    // on whether this is a parallel run.
    if (Pstream::parRun())
    {
        bb = decomposition().procBounds();
    }
    else
    {
        bb = geometryToConformTo().globalBounds();
    }

    scalar x0 = bb.min().x();
    scalar xR = bb.max().x() - x0;
    label ni = label(xR/initialCellSize_);

    scalar y0 = bb.min().y();
    scalar yR = bb.max().y() - y0;
    label nj = label(yR/initialCellSize_);

    scalar z0 = bb.min().z();
    scalar zR = bb.max().z() - z0;
    label nk = label(zR/initialCellSize_);

    vector delta(xR/ni, yR/nj, zR/nk);

    delta *= pow((1.0),-(1.0/3.0));

    scalar pert = randomPerturbationCoeff_*cmptMin(delta);

    // Initialise points list
    DynamicList<Vb::Point> initialPoints(scalar(ni)*nj*nk/10);

    for (label i = 0; i < ni; i++)
    {
        for (label j = 0; j < nj; j++)
        {
            // Generating, testing and adding points one line at a time to
            // reduce the memory requirement for cases with bounding boxes that
            // are very large in comparison to the volume to be filled

            label pI = 0;

            pointField points(nk);

            for (label k = 0; k < nk; k++)
            {
                point p
                (
                    x0 + (i + 0.5)*delta.x(),
                    y0 + (j + 0.5)*delta.y(),
                    z0 + (k + 0.5)*delta.z()
                );

                if (randomiseInitialGrid_)
                {
                    p.x() += pert*(rndGen().scalar01() - 0.5);
                    p.y() += pert*(rndGen().scalar01() - 0.5);
                    p.z() += pert*(rndGen().scalar01() - 0.5);
                }

                if
                (
                    Pstream::parRun()
                 && !decomposition().positionOnThisProcessor(p)
                )
                {
                    // Skip this point if, in parallel, this position is not on
                    // this processor.
                    continue;
                }

                points[pI++] = p;
            }

            points.setSize(pI);

            Field<bool> insidePoints =
                geometryToConformTo().wellInside
                (
                    points,
                    minimumSurfaceDistanceCoeffSqr_
                   *sqr
                    (
                        cellShapeControls().cellSize(points)
                    )
                );

            forAll(insidePoints, i)
            {
                if (insidePoints[i])
                {
                    const point& p(points[i]);

                    initialPoints.append(Vb::Point(p.x(), p.y(), p.z()));
                }
            }
        }
    }

    return initialPoints.shrink();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
