/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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

#include "pointFile.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointFile, 0);
addToRunTimeSelectionTable(initialPointsMethod, pointFile, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pointFile::pointFile
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
    pointFileName_(detailsDict().lookup("pointFile")),
    insideOutsideCheck_(detailsDict().lookup("insideOutsideCheck")),
    randomiseInitialGrid_(detailsDict().lookup("randomiseInitialGrid")),
    randomPerturbationCoeff_
    (
        detailsDict().lookup<scalar>("randomPerturbationCoeff")
    )
{
    Info<< "    Inside/Outside check is " << insideOutsideCheck_.asText()
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

List<Vb::Point> pointFile::initialPoints() const
{
    pointIOField points
    (
        IOobject
        (
            pointFileName_.name(),
            time().name(),
            time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "    Inserting points from file " << pointFileName_ << endl;

    if (points.empty())
    {
        FatalErrorInFunction
            << "Point file contain no points"
            << exit(FatalError) << endl;
    }

    if (Pstream::parRun())
    {
        // Testing filePath to see if the file originated in a processor
        // directory, if so, assume that the points in each processor file
        // are unique.  They are unlikely to belong on the current
        // processor as the background mesh is unlikely to be the same.

        const bool isParentFile = (points.objectPath() != points.filePath());

        if (!isParentFile)
        {
            decomposition().distributePoints(points);
        }
        else
        {
            // Otherwise, this is assumed to be points covering the whole
            // domain, so filter the points to be only those on this processor
            boolList procPt(decomposition().positionOnThisProcessor(points));

            List<boolList> allProcPt(Pstream::nProcs());

            allProcPt[Pstream::myProcNo()] = procPt;

            Pstream::gatherList(allProcPt);

            Pstream::scatterList(allProcPt);

            forAll(procPt, ptI)
            {
                bool foundAlready = false;

                forAll(allProcPt, proci)
                {
                    // If a processor with a lower index has found this point
                    // to insert already, defer to it and don't insert.
                    if (foundAlready)
                    {
                        allProcPt[proci][ptI] = false;
                    }
                    else if (allProcPt[proci][ptI])
                    {
                        foundAlready = true;
                    }
                }
            }

            procPt = allProcPt[Pstream::myProcNo()];

            inplaceSubset(procPt, points);
        }
    }

    Field<bool> insidePoints(points.size(), true);

    if (insideOutsideCheck_)
    {
        insidePoints = geometryToConformTo().wellInside
        (
            points,
            minimumSurfaceDistanceCoeffSqr_
           *sqr(cellShapeControls().cellSize(points))
        );
    }

    DynamicList<Vb::Point> initialPoints(insidePoints.size()/10);

    forAll(insidePoints, i)
    {
        if (insidePoints[i])
        {
            point& p = points[i];

            if (randomiseInitialGrid_)
            {
                p.x() += randomPerturbationCoeff_*(rndGen().scalar01() - 0.5);
                p.y() += randomPerturbationCoeff_*(rndGen().scalar01() - 0.5);
                p.z() += randomPerturbationCoeff_*(rndGen().scalar01() - 0.5);
            }

            initialPoints.append(Vb::Point(p.x(), p.y(), p.z()));
        }
    }

    initialPoints.shrink();

    label nPointsRejected = points.size() - initialPoints.size();

    if (Pstream::parRun())
    {
        reduce(nPointsRejected, sumOp<label>());
    }

    if (nPointsRejected)
    {
        Info<< "    " << nPointsRejected << " points rejected from "
            << pointFileName_.name() << endl;
    }

    return Foam::move(initialPoints);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
