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

Application
    mergeOrSplitBaffles

Description
    Detects faces that share points (baffles). Either merge them or
    duplicate the points.

    Notes:
    - can only handle pairwise boundary faces. So three faces using
      the same points is not handled (is illegal mesh anyway)

    - there is no option to only split/merge some baffles.

    - surfaces consisting of duplicate faces can be topologically split
    if the points on the interior of the surface cannot walk to all the
    cells that use them in one go.

    - Parallel operation (where duplicate face is perpendicular to a coupled
    boundary) is supported but not really tested.
    (Note that coupled faces themselves are not seen as duplicate faces)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "syncTools.H"
#include "faceSet.H"
#include "pointSet.H"
#include "meshTools.H"
#include "polyTopoChange.H"
#include "polyRemoveFace.H"
#include "polyModifyFace.H"
#include "indirectPrimitivePatch.H"
#include "processorPolyPatch.H"
#include "localPointRegion.H"
#include "duplicatePoints.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void insertDuplicateMerge
(
    const polyMesh& mesh,
    const labelList& duplicates,
    polyTopoChange& meshMod
)
{
    const faceList& faces = mesh.faces();
    const labelList& faceOwner = mesh.faceOwner();
    const faceZoneMesh& faceZones = mesh.faceZones();

    forAll(duplicates, bFacei)
    {
        label otherFacei = duplicates[bFacei];

        if (otherFacei != -1 && otherFacei > bFacei)
        {
            // Two duplicate faces. Merge.

            label face0 = mesh.nInternalFaces() + bFacei;
            label face1 = mesh.nInternalFaces() + otherFacei;

            label own0 = faceOwner[face0];
            label own1 = faceOwner[face1];

            if (own0 < own1)
            {
                // Use face0 as the new internal face.
                label zoneID = faceZones.whichZone(face0);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
                }

                meshMod.setAction(polyRemoveFace(face1));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face0],           // modified face
                        face0,                  // label of face being modified
                        own0,                   // owner
                        own1,                   // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
            else
            {
                // Use face1 as the new internal face.
                label zoneID = faceZones.whichZone(face1);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
                }

                meshMod.setAction(polyRemoveFace(face0));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face1],           // modified face
                        face1,                  // label of face being modified
                        own1,                   // owner
                        own0,                   // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
        }
    }
}


labelList findBaffles(const polyMesh& mesh, const labelList& boundaryFaces)
{
    // Get all duplicate face labels (in boundaryFaces indices!).
    labelList duplicates = localPointRegion::findDuplicateFaces
    (
        mesh,
        boundaryFaces
    );


    // Check that none are on processor patches
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(duplicates, bFacei)
    {
        if (duplicates[bFacei] != -1)
        {
            label facei = mesh.nInternalFaces() + bFacei;
            label patchi = patches.whichPatch(facei);

            if (isA<processorPolyPatch>(patches[patchi]))
            {
                FatalErrorInFunction
                    << "Duplicate face " << facei
                    << " is on a processorPolyPatch."
                    << "This is not allowed." << nl
                    << "Face:" << facei
                    << " is on patch:" << patches[patchi].name()
                    << abort(FatalError);
            }
        }
    }


    // Write to faceSet for ease of postprocessing.
    {
        faceSet duplicateSet
        (
            mesh,
            "duplicateFaces",
            (mesh.nFaces() - mesh.nInternalFaces())/256
        );

        forAll(duplicates, bFacei)
        {
            label otherFacei = duplicates[bFacei];

            if (otherFacei != -1 && otherFacei > bFacei)
            {
                duplicateSet.insert(mesh.nInternalFaces() + bFacei);
                duplicateSet.insert(mesh.nInternalFaces() + otherFacei);
            }
        }

        Pout<< "Writing " << duplicateSet.size()
            << " duplicate faces to faceSet " << duplicateSet.objectPath()
            << nl << endl;
        duplicateSet.write();
    }

    return duplicates;
}




int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Detect faces that share points (baffles).\n"
        "Merge them or duplicate the points."
    );

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "detectOnly",
        "find baffles only, but do not merge or split them"
    );
    argList::addBoolOption
    (
        "split",
        "topologically split duplicate surfaces"
    );
    argList::addBoolOption
    (
        "fields",
        "update fields"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createNamedMesh.H"

    const word oldInstance = mesh.pointsInstance();

    const bool split      = args.optionFound("split");
    const bool overwrite  = args.optionFound("overwrite");
    const bool detectOnly = args.optionFound("detectOnly");
    const bool fields     = args.optionFound("fields");

    // Collect all boundary faces
    labelList boundaryFaces(mesh.nFaces() - mesh.nInternalFaces());

    forAll(boundaryFaces, i)
    {
        boundaryFaces[i] = i+mesh.nInternalFaces();
    }


    if (detectOnly)
    {
        findBaffles(mesh, boundaryFaces);
        return 0;
    }


    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    if (fields) Info<< "Reading geometric fields" << nl << endl;

    PtrList<volScalarField> vsFlds;
    if (fields) ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    if (fields) ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    if (fields) ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    if (fields) ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    if (fields) ReadFields(mesh, objects, vtFlds);

    PtrList<surfaceScalarField> ssFlds;
    if (fields) ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    if (fields) ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    if (fields) ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    if (fields) ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    if (fields) ReadFields(mesh, objects, stFlds);


    // Mesh change engine
    polyTopoChange meshMod(mesh);


    if (split)
    {
        Pout<< "Topologically splitting duplicate surfaces"
            << ", i.e. duplicating points internal to duplicate surfaces."
            << nl << endl;

        // Analyse which points need to be duplicated
        localPointRegion regionSide(mesh);

        // Point duplication engine
        duplicatePoints pointDuplicator(mesh);

        // Insert topo changes
        pointDuplicator.setRefinement(regionSide, meshMod);
    }
    else
    {
        Pout<< "Merging duplicate faces."
            << nl << endl;

        // Get all duplicate face labels (in boundaryFaces indices!).
        labelList duplicates(findBaffles(mesh, boundaryFaces));

        // Merge into internal faces.
        insertDuplicateMerge(mesh, duplicates, meshMod);
    }

    if (!overwrite)
    {
        runTime++;
    }

    // Change the mesh. No inflation.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    Pout<< "Writing mesh to time " << runTime.timeName() << endl;
    mesh.write();

    // Dump duplicated points (if any)
    if (split)
    {
        const labelList& pointMap = map().pointMap();

        labelList nDupPerPoint(map().nOldPoints(), 0);

        pointSet dupPoints(mesh, "duplicatedPoints", 100);

        forAll(pointMap, pointi)
        {
            label oldPointi = pointMap[pointi];

            nDupPerPoint[oldPointi]++;

            if (nDupPerPoint[oldPointi] > 1)
            {
                dupPoints.insert(map().reversePointMap()[oldPointi]);
                dupPoints.insert(pointi);
            }
        }

        Pout<< "Writing " << dupPoints.size()
            << " duplicated points to pointSet "
            << dupPoints.objectPath() << nl << endl;

        dupPoints.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
