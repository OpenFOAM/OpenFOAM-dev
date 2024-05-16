/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "refinementParameters.H"
#include "unitConversion.H"
#include "polyMesh.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementParameters::refinementParameters(const dictionary& dict)
:
    maxGlobalCells_(dict.lookup<label>("maxGlobalCells")),
    maxLocalCells_(dict.lookup<label>("maxLocalCells")),
    minRefineCells_(dict.lookup<label>("minRefinementCells")),
    planarAngle_
    (
        dict.lookupOrDefault
        (
            "planarAngle",
            unitDegrees,
            dict.lookup<scalar>("resolveFeatureAngle", unitDegrees)
        )
    ),
    nBufferLayers_(dict.lookup<label>("nCellsBetweenLevels")),
    selectionPoints_(dict),
    allowFreeStandingZoneFaces_(dict.lookup("allowFreeStandingZoneFaces")),
    useTopologicalSnapDetection_
    (
        dict.lookupOrDefault<bool>("useTopologicalSnapDetection", true)
    ),
    maxLoadUnbalance_(dict.lookupOrDefault<scalar>("maxLoadUnbalance", 0)),
    handleSnapProblems_
    (
        dict.lookupOrDefault<Switch>("handleSnapProblems", true)
    )
{
    scalar featAngle(dict.lookup<scalar>("resolveFeatureAngle", unitDegrees));

    if (featAngle < 0 || featAngle > degToRad(180))
    {
        curvature_ = -great;
    }
    else
    {
        curvature_ = Foam::cos(featAngle);
    }
}


Foam::refinementParameters::cellSelectionPoints::cellSelectionPoints
(
    const dictionary& dict
)
:
    inside_
    (
        dict.found("insidePoints")
      ? dict.lookup<List<point>>("insidePoints", dimLength)
      : dict.found("insidePoint")
      ? List<point>(1, dict.lookup<point>("insidePoint", dimLength))
      : dict.found("locationInMesh")
      ? List<point>(1, dict.lookup<point>("locationInMesh", dimLength))
      : List<point>::null()
    ),
    outside_
    (
        dict.found("outsidePoints")
      ? dict.lookup<List<point>>("outsidePoints", dimLength)
      : dict.found("outsidePoint")
      ? List<point>(1, dict.lookup<point>("outsidePoint", dimLength))
      : List<point>::null()
    )
{
    if (inside_.size())
    {
        Info << "Cell selection insidePoints: " << inside_ << endl;
    }

    if (outside_.size())
    {
        Info << "Cell selection outsidePoints: " << outside_ << endl;
    }

    if (!inside_.size() && !outside_.size())
    {
        FatalErrorInFunction
            << "Neither insidePoint/insidePoints nor "
               "outsidePoint/outsidePoints specified: "
            << "cannot select any cells."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::refinementParameters::findCells(const polyMesh& mesh)
const
{
    // Force calculation of tet-diag decomposition (for use in findCell)
    (void)mesh.tetBasePtIs();

    // Global calculation engine
    globalIndex globalCells(mesh.nCells());

    // Cell label per point
    labelList cellLabels(selectionPoints_.inside().size());

    forAll(selectionPoints_.inside(), i)
    {
        const point& insidePoint = selectionPoints_.inside()[i];

        label localCelli = mesh.findCell(insidePoint);

        label globalCelli = -1;

        if (localCelli != -1)
        {
            globalCelli = globalCells.toGlobal(localCelli);
        }

        reduce(globalCelli, maxOp<label>());

        if (globalCelli == -1)
        {
            FatalErrorInFunction
                << "Point " << insidePoint
                << " is not inside the mesh or on a face or edge." << nl
                << "Bounding box of the mesh:" << mesh.bounds()
                << exit(FatalError);
        }


        label proci = globalCells.whichProcID(globalCelli);
        label procCelli = globalCells.toLocal(proci, globalCelli);

        Info<< "Found point " << insidePoint << " in cell " << procCelli
            << " on processor " << proci << endl;


        if (globalCells.isLocal(globalCelli))
        {
            cellLabels[i] = localCelli;
        }
        else
        {
            cellLabels[i] = -1;
        }
    }
    return cellLabels;
}


// ************************************************************************* //
