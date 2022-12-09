/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

Description
    Update the polyMesh corresponding to the given map.

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"
#include "polyTopoChangeMap.H"
#include "Time.H"
#include "globalMeshData.H"
#include "pointMesh.H"
#include "indexedOctree.H"
#include "treeDataCell.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMesh::topoChange(const polyTopoChangeMap& map)
{
    if (debug)
    {
        InfoInFunction
            << "Updating addressing and (optional) pointMesh/pointFields"
            << endl;
    }

    // Update boundaryMesh (note that patches themselves already ok)
    boundary_.topoChange();

    // Update zones
    pointZones_.clearAddressing();
    faceZones_.clearAddressing();
    cellZones_.clearAddressing();

    // Remove the stored tet base points
    tetBasePtIsPtr_.clear();

    // Remove the cell tree
    cellTreePtr_.clear();

    // Update parallel data
    if (globalMeshDataPtr_.valid())
    {
        globalMeshDataPtr_->topoChange();
    }

    setInstance(time().name());

    // Map the old motion points if present
    if (oldPointsPtr_.valid())
    {
        // Make a copy of the original points
        pointField oldMotionPoints = oldPointsPtr_();

        pointField& newMotionPoints = oldPointsPtr_();

        // Resize the list to new size
        newMotionPoints.setSize(points_.size());

        // Map the list
        newMotionPoints.map(oldMotionPoints, map.pointMap());

        // Any points created out-of-nothing get set to the current coordinate
        // for lack of anything better.
        forAll(map.pointMap(), newPointi)
        {
            if (map.pointMap()[newPointi] == -1)
            {
                newMotionPoints[newPointi] = points_[newPointi];
            }
        }
    }

    if (oldCellCentresPtr_.valid())
    {
        // Make a copy of the original cell-centres
        pointField oldMotionCellCentres = oldCellCentresPtr_();

        pointField& newMotionCellCentres = oldCellCentresPtr_();

        // Resize the list to new size
        newMotionCellCentres.setSize(cellCentres().size());

        // Map the list
        newMotionCellCentres.map(oldMotionCellCentres, map.cellMap());

        // Any points created out-of-nothing get set to the current coordinate
        // for lack of anything better.
        forAll(map.cellMap(), newCelli)
        {
            if (map.cellMap()[newCelli] == -1)
            {
                newMotionCellCentres[newCelli] = cellCentres()[newCelli];
            }
        }
    }

    meshObjects::topoChange<polyMesh>(*this, map);
    meshObjects::topoChange<pointMesh>(*this, map);

    // Reset valid directions (could change by faces put into empty patches)
    geometricD_ = Zero;
    solutionD_ = Zero;
}


void Foam::polyMesh::mapMesh(const polyMeshMap& map)
{
    meshObjects::mapMesh<polyMesh>(*this, map);
    meshObjects::mapMesh<pointMesh>(*this, map);
}


void Foam::polyMesh::distribute(const polyDistributionMap& map)
{
    meshObjects::distribute<polyMesh>(*this, map);
    meshObjects::distribute<pointMesh>(*this, map);
}


// ************************************************************************* //
