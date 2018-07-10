/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    maxGlobalCells_(readLabel(dict.lookup("maxGlobalCells"))),
    maxLocalCells_(readLabel(dict.lookup("maxLocalCells"))),
    minRefineCells_(readLabel(dict.lookup("minRefinementCells"))),
    planarAngle_
    (
        dict.lookupOrDefault
        (
            "planarAngle",
            readScalar(dict.lookup("resolveFeatureAngle"))
        )
    ),
    nBufferLayers_(readLabel(dict.lookup("nCellsBetweenLevels"))),
    keepPoints_(pointField(1, dict.lookup("locationInMesh"))),
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
    scalar featAngle(readScalar(dict.lookup("resolveFeatureAngle")));

    if (featAngle < 0 || featAngle > 180)
    {
        curvature_ = -great;
    }
    else
    {
        curvature_ = Foam::cos(degToRad(featAngle));
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
    labelList cellLabels(keepPoints_.size());

    forAll(keepPoints_, i)
    {
        const point& keepPoint = keepPoints_[i];

        label localCelli = mesh.findCell(keepPoint);

        label globalCelli = -1;

        if (localCelli != -1)
        {
            globalCelli = globalCells.toGlobal(localCelli);
        }

        reduce(globalCelli, maxOp<label>());

        if (globalCelli == -1)
        {
            FatalErrorInFunction
                << "Point " << keepPoint
                << " is not inside the mesh or on a face or edge." << nl
                << "Bounding box of the mesh:" << mesh.bounds()
                << exit(FatalError);
        }


        label proci = globalCells.whichProcID(globalCelli);
        label procCelli = globalCells.toLocal(proci, globalCelli);

        Info<< "Found point " << keepPoint << " in cell " << procCelli
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
