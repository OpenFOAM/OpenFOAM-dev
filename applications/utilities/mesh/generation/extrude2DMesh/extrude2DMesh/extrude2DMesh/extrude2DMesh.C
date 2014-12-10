/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "extrude2DMesh.H"
#include "polyMesh.H"
#include "polyTopoChange.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extrude2DMesh, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::extrude2DMesh::check2D() const
{
    const faceList& faces = mesh_.faces();
    forAll(faces, faceI)
    {
        if (faces[faceI].size() != 2)
        {
            FatalErrorIn("void Foam::extrude2DMesh::check2D() const")
                << "Face " << faceI << " size " << faces[faceI].size()
                << " is not of size 2: mesh is not a valid two-dimensional "
                << "mesh" << exit(FatalError);
        }
    }
}


//void Foam::extrude2DMesh::findExtrudeDirection()
//{
//    scalar minRange = GREAT;

//    for (direction dir = 0; dir < 3; dir++)
//    {
//        scalarField cmpts(mesh_.points().component(dir));

//        scalar range = max(cmpts)-min(cmpts);

//        Info<< "Direction:" << dir << " range:" << range << endl;

//        if (range < minRange)
//        {
//            minRange = range;
//            extrudeDir_ = dir;
//        }
//    }

//    Info<< "Extruding in direction " << extrudeDir_
//        << " with thickness " << thickness_ << nl
//        << endl;
//}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extrude2DMesh::extrude2DMesh
(
    polyMesh& mesh,
    const dictionary& dict,
    const extrudeModel& model
)
:
    mesh_(mesh),
    dict_(dict),
    //patchDict_(dict.subDict("patchInfo")),
    model_(model),
    modelType_(dict.lookup("extrudeModel")),
    patchType_(dict.lookup("patchType")),
    frontPatchI_(-1),
    backPatchI_(-1)
{
    check2D();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extrude2DMesh::~extrude2DMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::extrude2DMesh::addFrontBackPatches()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    frontPatchI_ = patches.findPatchID("front");
    backPatchI_ = patches.findPatchID("back");

    // Add patch.
    List<polyPatch*> newPatches(patches.size() + 2);

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        newPatches[patchI] =
            pp.clone
            (
                patches,
                newPatches.size(),
                pp.size(),
                pp.start()
            ).ptr();
    }

    if (frontPatchI_ == -1)
    {
        frontPatchI_ = patches.size();

        newPatches[frontPatchI_] =
            polyPatch::New
            (
                patchType_,
                "front",
                0,
                mesh_.nFaces(),
                frontPatchI_,
                patches
            ).ptr();

//        newPatches[frontPatchI_] = polyPatch::New
//        (
//            "front",
//            patchDict_,
//            frontPatchI_,
//            patches
//        ).ptr();

        Info<< "Adding patch " << newPatches[frontPatchI_]->name()
            << " at index " << frontPatchI_
            << " for front faces." << nl << endl;
    }

    if (backPatchI_ == -1)
    {
        backPatchI_ = patches.size() + 1;

        newPatches[backPatchI_] =
            polyPatch::New
            (
                patchType_,
                "back",
                0,
                mesh_.nFaces(),
                backPatchI_,
                patches
            ).ptr();

//        newPatches[frontPatchI_] = polyPatch::New
//        (
//            "back",
//            patchDict_,
//            backPatchI_,
//            patches
//        ).ptr();

        Info<< "Adding patch " << newPatches[backPatchI_]->name()
            << " at index " << backPatchI_
            << " for back faces." << nl << endl;
    }

    mesh_.removeBoundary();
    mesh_.addPatches(newPatches);
}


void Foam::extrude2DMesh::setRefinement
(
    polyTopoChange& meshMod
)
{
    const label nLayers = model_.nLayers();
    const pointField& points = mesh_.points();
    label nFaces = 0;

    for (label layer = 0; layer < nLayers; ++layer)
    {
        label offset = layer * mesh_.nCells();

        forAll(mesh_.cells(), cellI)
        {
            meshMod.addCell
            (
                -1,     //masterPointID,
                -1,     //masterEdgeID,
                -1,     //masterFaceID,
                cellI + offset,  //masterCellID,
                mesh_.cellZones().whichZone(cellI)  //zoneID
            );
        }
    }


    // Generate points
    // ~~~~~~~~~~~~~~~

    for (label layer = 0; layer <= nLayers; ++layer)
    {
        label offset = layer * points.size();

        forAll(points, pointI)
        {
            // Don't need the surface normal for either linearDirection or
            // wedge. Will need to add to be able to use others.
            point newPoint = model_
            (
                points[pointI],
                vector(),
                layer
            );

            meshMod.addPoint
            (
                newPoint,
                pointI + offset,
                -1,             // zoneID
                true            // inCell
            );
        }

        Pout<< "Added " << points.size() << " points to layer "
            << layer << endl;
    }


    // Generate faces
    // ~~~~~~~~~~~~~~

    const faceList& faces = mesh_.faces();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    for (label layer = 0; layer < nLayers; ++layer)
    {
        label currentLayerOffset = layer * mesh_.nPoints();
        label nextLayerOffset = currentLayerOffset + mesh_.nPoints();

        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label zoneID = mesh_.faceZones().whichZone(faceI);
            bool zoneFlip = false;
            if (zoneID != -1)
            {
                const faceZone& fZone = mesh_.faceZones()[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
            }

            face newFace(4);
            const face& f = faces[faceI];
            newFace[0] = f[0] + currentLayerOffset;
            newFace[1] = f[1] + currentLayerOffset;
            newFace[2] = f[1] + nextLayerOffset;
            newFace[3] = f[0] + nextLayerOffset;

//{
//    vector n = newFace.normal(pointField(meshMod.points()));
//    label own = mesh_.faceOwner()[faceI];
//    const labelList& ownPoints = mesh_.cellPoints()[own];
//    point ownCc = sum(pointField(mesh_.points(), ownPoints))/ownPoints.size();
//    label nei = mesh_.faceNeighbour()[faceI];
//    const labelList& neiPoints = mesh_.cellPoints()[nei];
//    point neiCc = sum(pointField(mesh_.points(), neiPoints))/neiPoints.size();
//    vector d = neiCc - ownCc;

//    Pout<< "face:" << faceI << " at:" << f.centre(mesh_.points()) << endl
//        << "    own:" << own << " at:" << ownCc << endl
//        << "    nei:" << nei << " at:" << neiCc << endl
//        << "    sign:" << (n & d) << endl
//        << endl;
//}

            label offset = layer * mesh_.nCells();

            meshMod.addFace
            (
                newFace,
                mesh_.faceOwner()[faceI] + offset,       // own
                mesh_.faceNeighbour()[faceI] + offset,   // nei
                -1,                             // masterPointID
                -1,                             // masterEdgeID
                nFaces++,    // masterFaceID
                false,                          // flipFaceFlux
                -1,                             // patchID
                zoneID,                         // zoneID
                zoneFlip                        // zoneFlip
            );

            if (debug)
            {
                Info<< newFace << " "
                    << mesh_.faceOwner()[faceI] + offset << " "
                    << mesh_.faceNeighbour()[faceI] + offset << " "
                    << nFaces - 1
                    << endl;
            }
        }
    }

    forAll(patches, patchI)
    {
        for (label layer=0; layer < nLayers; layer++)
        {
            label currentLayerOffset = layer*mesh_.nPoints();
            label nextLayerOffset = currentLayerOffset + mesh_.nPoints();

            label startFaceI = patches[patchI].start();
            label endFaceI = startFaceI + patches[patchI].size();

            for (label faceI = startFaceI; faceI < endFaceI; faceI++)
            {
                label zoneID = mesh_.faceZones().whichZone(faceI);
                bool zoneFlip = false;
                if (zoneID != -1)
                {
                    const faceZone& fZone = mesh_.faceZones()[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                }

                face newFace(4);
                const face& f = faces[faceI];
                newFace[0] = f[0] + currentLayerOffset;
                newFace[1] = f[1] + currentLayerOffset;
                newFace[2] = f[1] + nextLayerOffset;
                newFace[3] = f[0] + nextLayerOffset;

                label offset = layer * mesh_.nCells();

                meshMod.addFace
                (
                    newFace,
                    mesh_.faceOwner()[faceI] + offset,       // own
                    -1,                                      // nei
                    -1,                                      // masterPointID
                    -1,                                      // masterEdgeID
                    nFaces++,                                // masterFaceID
                    false,                                   // flipFaceFlux
                    patchI,                                  // patchID
                    zoneID,                                  // zoneID
                    zoneFlip                                 // zoneFlip
                );

                if (debug)
                {
                    Info<< newFace << " "
                        << mesh_.faceOwner()[faceI] + offset << " "
                        << nFaces - 1
                        << endl;
                }
            }
        }
    }

    // Add extra internal faces that need special treatment for owners and
    // neighbours.
    forAll(mesh_.cells(), cellI)
    {
        const cell& cFaces = mesh_.cells()[cellI];

        face frontFace(cFaces.size());

        // Make a loop out of faces.
        label nextFaceI = cFaces[0];

        const face& f = faces[nextFaceI];

        label nextPointI;
        if (mesh_.faceOwner()[nextFaceI] == cellI)
        {
            frontFace[0] = f[0];
            nextPointI = f[1];
        }
        else
        {
            frontFace[0] = f[1];
            nextPointI = f[0];
        }


        for (label i = 1; i < frontFace.size(); i++)
        {
            frontFace[i] = nextPointI;

            // Find face containing pointI
            forAll(cFaces, cFaceI)
            {
                label faceI = cFaces[cFaceI];
                if (faceI != nextFaceI)
                {
                    const face& f = faces[faceI];

                    if (f[0] == nextPointI)
                    {
                        nextPointI = f[1];
                        nextFaceI = faceI;
                        break;
                    }
                    else if (f[1] == nextPointI)
                    {
                        nextPointI = f[0];
                        nextFaceI = faceI;
                        break;
                    }
                }
            }
        }

        for (label layer = 0; layer < nLayers - 1; ++layer)
        {
            // Offset to create front face.
            forAll(frontFace, fp)
            {
                frontFace[fp] += mesh_.nPoints();
            }

            label offset = layer * mesh_.nCells();

            label nei = -1;
            if (layer != nLayers - 1)
            {
                nei = cellI + offset + mesh_.nCells();
            }

            meshMod.addFace
            (
                frontFace,
                cellI + offset,                 // own
                nei,                            // nei
                -1,                             // masterPointID
                -1,                             // masterEdgeID
                nFaces++,                       // masterFaceID
                false,                          // flipFaceFlux
                -1,                             // patchID
                -1,                             // zoneID
                false                           // zoneFlip
            );

            if (debug)
            {
                Info<< frontFace << " "
                    << cellI + offset << " "
                    << nei << " "
                    << nFaces - 1
                    << endl;
            }
        }
    }

    // Generate front and back faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(mesh_.cells(), cellI)
    {
        const cell& cFaces = mesh_.cells()[cellI];

        face frontFace(cFaces.size());

        // Make a loop out of faces.
        label nextFaceI = cFaces[0];

        const face& f = faces[nextFaceI];

        label nextPointI;
        if (mesh_.faceOwner()[nextFaceI] == cellI)
        {
            frontFace[0] = f[0];
            nextPointI = f[1];
        }
        else
        {
            frontFace[0] = f[1];
            nextPointI = f[0];
        }


        for (label i = 1; i < frontFace.size(); i++)
        {
            frontFace[i] = nextPointI;

            // Find face containing pointI
            forAll(cFaces, cFaceI)
            {
                label faceI = cFaces[cFaceI];
                if (faceI != nextFaceI)
                {
                    const face& f = faces[faceI];

                    if (f[0] == nextPointI)
                    {
                        nextPointI = f[1];
                        nextFaceI = faceI;
                        break;
                    }
                    else if (f[1] == nextPointI)
                    {
                        nextPointI = f[0];
                        nextFaceI = faceI;
                        break;
                    }
                }
            }
        }

        // Add back face.
        meshMod.addFace
        (
            frontFace.reverseFace(),
            cellI,                          // own
            -1,                             // nei
            -1,                             // masterPointID
            -1,                             // masterEdgeID
            nFaces++,                       // masterFaceID
            false,                          // flipFaceFlux
            backPatchI_,                    // patchID
            -1,                             // zoneID
            false                           // zoneFlip
        );

        if (debug)
        {
            Info<< nl<<frontFace.reverseFace() << " "
                << cellI << " "
                << nFaces - 1
                << endl;
        }

        // Offset to create front face.
        forAll(frontFace, fp)
        {
            frontFace[fp] += mesh_.nPoints()* (nLayers);
        }

        label offset = (nLayers - 1) * mesh_.nCells();

        meshMod.addFace
        (
            frontFace,
            cellI + offset,                 // own
            -1,                             // nei
            -1,                             // masterPointID
            -1,                             // masterEdgeID
            nFaces++,                       // masterFaceID
            false,                          // flipFaceFlux
            frontPatchI_,                   // patchID
            -1,                             // zoneID
            false                           // zoneFlip
        );

        if (debug)
        {
            Info<< frontFace << " "
                << cellI + offset << " "
                << nFaces - 1
                << endl;
        }
    }
}


// ************************************************************************* //
