/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "regionSplit.H"
#include "FaceCellWave.H"
#include "globalIndex.H"
#include "minData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionSplit, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionSplit::calcNonCompactRegionSplit
(
    const globalIndex& globalFaces,
    const boolList& blockedFace,
    const List<labelPair>& explicitConnections,
    labelList& cellRegion,
    const bool global
) const
{
    // Initialise wave data. Blocked faces are given a negative number, so it
    // is always less than that of the global face indices, so the walk will
    // not cross them
    List<minData> faceData(mesh().nFaces());
    label nUnblocked = 0;
    if (blockedFace.empty())
    {
        nUnblocked += mesh().nFaces();
    }
    else
    {
        forAll(faceData, facei)
        {
            if (blockedFace[facei])
            {
                faceData[facei] = minData(-2);
            }
            else
            {
                nUnblocked++;
            }
        }
    }

    List<minData> cellData(mesh().nCells());

    // Seed all unblocked faces with a globally unique number
    labelList seedFaces(nUnblocked);
    List<minData> seedData(nUnblocked);
    auto& td = FaceCellWave<minData>::defaultTrackingData_;
    nUnblocked = 0;
    forAll(faceData, facei)
    {
        if (!faceData[facei].valid(td))
        {
            seedFaces[nUnblocked] = facei;
            seedData[nUnblocked] = minData(globalFaces.toGlobal(facei));
            nUnblocked++;
        }
    }

    // Propagate information inwards
    FaceCellWave<minData> deltaCalc
    (
        mesh(),
        explicitConnections,
        seedFaces,
        seedData,
        faceData,
        cellData,
        mesh().globalData().nTotalCells() + 1,
        FaceCellWave<minData>::defaultTrackingData_,
        global
    );

    // Extract into the cells
    cellRegion.setSize(mesh().nCells());
    forAll(cellRegion, celli)
    {
        if (cellData[celli].valid(deltaCalc.data()))
        {
            cellRegion[celli] = cellData[celli].data();
        }
        else
        {
            // Unvisited cell. This is only possible if it is surrounded by
            // blocked faces. If so make up region from any of the faces.

            const label facei = mesh().cells()[celli][0];

            if (blockedFace.size() && !blockedFace[facei])
            {
                FatalErrorInFunction
                    << "Unblocked face " << facei
                    << " at " << mesh().faceCentres()[facei]
                    << " on unassigned cell " << celli
                    << mesh().cellCentres()[celli]
                    << exit(FatalError);
            }

            cellRegion[celli] = globalFaces.toGlobal(facei);
        }
    }
}


Foam::label Foam::regionSplit::calcLocalRegionSplit
(
    const boolList& blockedFace,
    const List<labelPair>& explicitConnections,
    labelList& cellRegion
) const
{
    // Minimise across locally connected cells. This leaves cellRegion
    // containing the minimum face index in the contiguous region of mesh.
    calcNonCompactRegionSplit
    (
        localGlobalIndex(mesh().nFaces())(),
        blockedFace,
        explicitConnections,
        cellRegion,
        false
    );

    // Compact and return
    return compactLocalRegionSplit(cellRegion);
}


Foam::label Foam::regionSplit::calcGlobalRegionSplit
(
    const boolList& blockedFace,
    const List<labelPair>& explicitConnections,
    labelList& cellRegion
) const
{
    // Create global face index
    const globalIndex globalFaces(mesh().nFaces());

    // Minimise regions across connected cells. This leaves cellRegion
    // containing the minimum face index in the contiguous region of mesh.
    calcNonCompactRegionSplit
    (
        globalFaces,
        blockedFace,
        explicitConnections,
        cellRegion,
        true
    );

    // Compact and return
    return compactGlobalRegionSplit(globalFaces, cellRegion);
}


Foam::label Foam::regionSplit::calcRegionSplit
(
    const boolList& blockedFace,
    const List<labelPair>& explicitConnections,
    labelList& cellRegion,
    const bool global
) const
{
    return
        global
      ? calcGlobalRegionSplit(blockedFace, explicitConnections, cellRegion)
      : calcLocalRegionSplit(blockedFace, explicitConnections, cellRegion);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSplit::regionSplit(const polyMesh& mesh, const bool global)
:
    regionSplitBase(mesh.nCells()),
    mesh_(mesh)
{
    nRegions_ =
        calcRegionSplit
        (
            boolList(0), // blockedFaces
            List<labelPair>(0), // explicitConnections,
            *this,
            global
        );
}


Foam::regionSplit::regionSplit
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const bool global
)
:
    regionSplitBase(mesh.nCells()),
    mesh_(mesh)
{
    nRegions_ =
        calcRegionSplit
        (
            blockedFace,
            List<labelPair>(0), // explicitConnections,
            *this,
            global
        );
}


Foam::regionSplit::regionSplit
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const List<labelPair>& explicitConnections,
    const bool global
)
:
    regionSplitBase(mesh.nCells()),
    mesh_(mesh)
{
    nRegions_ =
        calcRegionSplit
        (
            blockedFace,
            explicitConnections,
            *this,
            global
        );
}


// ************************************************************************* //
