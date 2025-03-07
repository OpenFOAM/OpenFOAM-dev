/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "argList.H"
#include "IFstream.H"
#include "BSpline.H"
#include "CatmullRomSpline.H"
#include "vtkWritePolyData.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.insert("dictionary");
    argList args(argc, argv);

    const word dictName(args.argRead<word>(1));
    Info<< "Reading " << dictName << nl << endl;
    dictionary dict = IFstream(dictName)();

    // Controls
    const pointField points(dict.lookup<pointField>("points"));
    const scalar tol = dict.lookupOrDefault<scalar>("tol", rootSmall);
    const label nIter = dict.lookupOrDefault<label>("nIter", 16);
    const scalar grading = dict.lookup<scalar>("grading");
    const label nSplinePoints = dict.lookup<label>("n");

    // Spline classes
    const BSpline bSpline(points, tol, nIter);
    const CatmullRomSpline crSpline(points, tol, nIter);

    // Parameters
    const scalarField is(scalarField(scalarList(identityMap(nSplinePoints))));
    const scalarField lambdas
    (
        mag(grading - 1) < small
      ? is/(is.size() - 1)
      : (1 - pow(grading, is/(is.size() - 2)))
       /(1 - Foam::pow(grading, (is.size() - 1)/scalar(is.size() - 2)))
    );

    // Evaluate
    const pointField bSplinePoints(bSpline.position(lambdas));
    const pointField crSplinePoints(crSpline.position(lambdas));

    // Write
    pointField allPoints;
    allPoints.append(points);
    allPoints.append(bSplinePoints);
    allPoints.append(crSplinePoints);
    labelListList allLines;
    allLines.append(identityMap(points.size()));
    allLines.append(points.size() + identityMap(nSplinePoints));
    allLines.append(points.size() + nSplinePoints + identityMap(nSplinePoints));
    vtkWritePolyData::write
    (
        dictName + ".vtk",
        dictName,
        false,
        allPoints,
        labelList(),
        allLines,
        labelListList()
    );

    return 0;
}


// ************************************************************************* //
