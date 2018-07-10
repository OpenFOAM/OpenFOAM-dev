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

#include "duplicatePoints.H"
#include "localPointRegion.H"
#include "polyTopoChange.H"
#include "polyAddPoint.H"
#include "polyModifyFace.H"
#include "polyMesh.H"
#include "OFstream.H"
#include "meshTools.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(duplicatePoints, 0);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::duplicatePoints::duplicatePoints(const polyMesh& mesh)
:
    mesh_(mesh),
    duplicates_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::duplicatePoints::setRefinement
(
    const localPointRegion& regionSide,
    polyTopoChange& meshMod
)
{
    const Map<label>& meshPointMap = regionSide.meshPointMap();
    const labelListList& pointRegions = regionSide.pointRegions();
    const Map<label>& meshFaceMap = regionSide.meshFaceMap();
    const faceList& faceRegions = regionSide.faceRegions();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Create duplicates for points. One for each region.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Per point-to-be-duplicated, in order of the regions the point added.
    duplicates_.setSize(meshPointMap.size());

    forAllConstIter(Map<label>, meshPointMap, iter)
    {
        label pointi = iter.key();
        label localI = iter();
        const labelList& regions = pointRegions[localI];

        duplicates_[localI].setSize(regions.size());
        duplicates_[localI][0] = pointi;
        for (label i = 1; i < regions.size(); i++)
        {
            duplicates_[localI][i] = meshMod.addPoint
            (
                mesh_.points()[pointi],  // point
                pointi,                 // master point
                -1,                     // zone for point
                true                    // supports a cell
            );
        }

        // Pout<< "For point:" << pointi << " coord:" << mesh_.points()[pointi]
        //    << endl;
        // forAll(duplicates_[localI], i)
        //{
        //    Pout<< "    region:" << regions[i]
        //        << "  addedpoint:" << duplicates_[localI][i]
        //        << endl;
        //}
    }



    // Modfify faces according to face region
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    face newFace;

    forAllConstIter(Map<label>, meshFaceMap, iter)
    {
        label facei = iter.key();
        label localI = iter();

        // Replace points with duplicated ones.
        const face& fRegion = faceRegions[localI];
        const face& f = mesh_.faces()[facei];

        newFace.setSize(f.size());
        forAll(f, fp)
        {
            label pointi = f[fp];

            Map<label>::const_iterator iter = meshPointMap.find(pointi);

            if (iter != meshPointMap.end())
            {
                // Point has been duplicated. Find correct one for my
                // region.

                // Get the regions and added points for this point
                const labelList& regions = pointRegions[iter()];
                const labelList& dupPoints = duplicates_[iter()];

                // Look up index of my region in the regions for this point
                label index = findIndex(regions, fRegion[fp]);
                // Get the corresponding added point
                newFace[fp] = dupPoints[index];
            }
            else
            {
                newFace[fp] = pointi;
            }
        }

        // Get current zone info
        label zoneID = mesh_.faceZones().whichZone(facei);
        bool zoneFlip = false;
        if (zoneID >= 0)
        {
            const faceZone& fZone = mesh_.faceZones()[zoneID];
            zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
        }


        if (mesh_.isInternalFace(facei))
        {
            meshMod.modifyFace
            (
                newFace,                    // modified face
                facei,                      // label of face being modified
                mesh_.faceOwner()[facei],   // owner
                mesh_.faceNeighbour()[facei],   // neighbour
                false,                      // face flip
                -1,                         // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            );
        }
        else
        {
            meshMod.modifyFace
            (
                newFace,                    // modified face
                facei,                      // label of face being modified
                mesh_.faceOwner()[facei],   // owner
                -1,                         // neighbour
                false,                      // face flip
                patches.whichPatch(facei),  // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            );
        }
    }


    if (debug)
    {
        // Output duplicated points
        {
            OFstream str(mesh_.time().path()/"duplicatedPoints.obj");
            forAllConstIter(Map<label>, meshPointMap, iter)
            {
                label localI = iter();
                const labelList& dups = duplicates_[localI];

                forAll(dups, i)
                {
                    meshTools::writeOBJ(str, meshMod.points()[dups[i]]);
                }
            }
        }
    }
}


void Foam::duplicatePoints::updateMesh(const mapPolyMesh& map)
{
    forAll(duplicates_, masterI)
    {
        inplaceRenumber(map.reversePointMap(), duplicates_[masterI]);
    }
}


// ************************************************************************* //
