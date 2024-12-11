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

#include "mergePolyMesh.H"
#include "Time.H"
#include "polyTopoChangeMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mergePolyMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::mergePolyMesh::patchIndex(const polyPatch& p)
{
    // Find the patch name on the list.  If the patch is already there
    // and patch types match, return index
    const word& pType = p.type();
    const word& pName = p.name();

    bool nameFound = false;

    forAll(patchNames_, patchi)
    {
        if (patchNames_[patchi] == pName)
        {
            if (word(patchDicts_[patchi]["type"]) == pType)
            {
                // Found name and types match
                return patchi;
            }
            else
            {
                // Found the name, but type is different
                nameFound = true;
            }
        }
    }

    // Patch not found.  Append to the list
    {
        OStringStream os;
        p.write(os);
        patchDicts_.append(dictionary(IStringStream(os.str())()));
    }

    if (nameFound)
    {
        // Duplicate name is not allowed.  Create a composite name from the
        // patch name and case name
        const word& caseName = p.boundaryMesh().mesh().time().caseName();

        patchNames_.append(pName + "_" + caseName);

        Info<< "label patchIndex(const polyPatch& p) : "
            << "Patch " << p.index() << " named "
            << pName << " in mesh " << caseName
            << " already exists, but patch types "
            << " do not match.\nCreating a composite name as "
            << patchNames_.last() << endl;
    }
    else
    {
        patchNames_.append(pName);
    }

    return patchNames_.size() - 1;
}


Foam::label Foam::mergePolyMesh::zoneIndex
(
    DynamicList<word>& names,
    const word& curName
)
{
    forAll(names, zonei)
    {
        if (names[zonei] == curName)
        {
            return zonei;
        }
    }

    // Not found.  Add new name to the list
    names.append(curName);

    return names.size() - 1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mergePolyMesh::mergePolyMesh(polyMesh& mesh)
:
    mesh_(mesh),
    meshMod_(mesh_),
    patchNames_(2*mesh_.boundaryMesh().size()),
    patchDicts_(2*mesh_.boundaryMesh().size()),
    pointZoneNames_(),
    faceZoneNames_(),
    cellZoneNames_()
{
    // Insert the original patches into the list
    wordList curPatchNames = mesh_.boundaryMesh().names();

    forAll(mesh_.boundaryMesh(), patchi)
    {
        patchNames_.append(mesh_.boundaryMesh()[patchi].name());

        OStringStream os;
        mesh_.boundaryMesh()[patchi].write(os);
        patchDicts_.append(dictionary(IStringStream(os.str())()));
    }

    // Insert point, face and cell zones into the list

    // Point zones
    wordList curPointZoneNames = mesh_.pointZones().toc();
    if (curPointZoneNames.size())
    {
        pointZoneNames_.setCapacity(2*curPointZoneNames.size());
    }

    forAll(curPointZoneNames, zonei)
    {
        pointZoneNames_.append(curPointZoneNames[zonei]);
    }

    pointZonesAddedPoints_.setSize(pointZoneNames_.size());

    // Face zones
    wordList curFaceZoneNames = mesh_.faceZones().toc();

    if (curFaceZoneNames.size())
    {
        faceZoneNames_.setCapacity(2*curFaceZoneNames.size());
    }

    forAll(curFaceZoneNames, zonei)
    {
        faceZoneNames_.append(curFaceZoneNames[zonei]);
    }

    faceZonesAddedFaces_.setSize(faceZoneNames_.size());

    // Cell zones
    wordList curCellZoneNames = mesh_.cellZones().toc();

    if (curCellZoneNames.size())
    {
        cellZoneNames_.setCapacity(2*curCellZoneNames.size());
    }

    forAll(curCellZoneNames, zonei)
    {
        cellZoneNames_.append(curCellZoneNames[zonei]);
    }

    cellZonesAddedCells_.setSize(cellZoneNames_.size());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mergePolyMesh::addMesh(const polyMesh& m)
{
    // Add all the points, faces and cells of the new mesh

    // Add points

    const pointField& p = m.points();
    labelList renumberPoints(p.size());

    const pointZoneList& pz = m.pointZones();
    labelList pointZoneIndices(pz.size());

    forAll(pz, zonei)
    {
        pointZoneIndices[zonei] = zoneIndex(pointZoneNames_, pz[zonei].name());
        pointZonesAddedPoints_.setSize(pointZoneNames_.size());
    }

    forAll(p, pointi)
    {
        renumberPoints[pointi] = meshMod_.addPoint
        (
            p[pointi],            // Point to add
            -1,                   // Master point (straight addition)
            pointi < m.nPoints()  // Is in cell?
        );

        const labelList zones(pz.whichZones(pointi));
        forAll(zones, zonei)
        {
            pointZonesAddedPoints_[pointZoneIndices[zones[zonei]]]
           .insert(renumberPoints[pointi]);
        }
    }

    // Add cells

    const cellList& c = m.cells();
    labelList renumberCells(c.size());

    const cellZoneList& cz = m.cellZones();
    labelList cellZoneIndices(cz.size());

    forAll(cz, zonei)
    {
        cellZoneIndices[zonei] = zoneIndex(cellZoneNames_, cz[zonei].name());
        cellZonesAddedCells_.setSize(cellZoneNames_.size());
    }

    forAll(c, celli)
    {
        renumberCells[celli] = meshMod_.addCell(-1);

        const labelList zones(cz.whichZones(celli));
        forAll(zones, zonei)
        {
            cellZonesAddedCells_[cellZoneIndices[zones[zonei]]]
           .insert(renumberCells[celli]);
        }
    }

    // Add faces
    const polyBoundaryMesh& bm = m.boundaryMesh();

    // Gather the patch indices
    labelList patchIndices(bm.size());

    forAll(patchIndices, patchi)
    {
        patchIndices[patchi] = patchIndex(bm[patchi]);
    }

    // Temporary: update number of allowable patches. This should be
    // determined at the top - before adding anything.
    meshMod_.setNumPatches(patchNames_.size());



    const faceZoneList& fz = m.faceZones();
    labelList faceZoneIndices(fz.size());

    forAll(fz, zonei)
    {
        faceZoneIndices[zonei] = zoneIndex(faceZoneNames_, fz[zonei].name());
        faceZonesAddedFaces_.setSize(faceZoneNames_.size());
    }

    const faceList& f = m.faces();
    labelList renumberFaces(f.size());

    const labelList& own = m.faceOwner();
    const labelList& nei = m.faceNeighbour();

    label newOwn, newNei, newPatch;

    forAll(f, facei)
    {
        const face& curFace = f[facei];

        face newFace(curFace.size());

        forAll(curFace, pointi)
        {
            newFace[pointi] = renumberPoints[curFace[pointi]];
        }

        if (debug)
        {
            // Check that the face is valid
            if (min(newFace) < 0)
            {
                FatalErrorInFunction
                    << "Error in point mapping for face " << facei
                    << ".  Old face: " << curFace << " New face: " << newFace
                    << abort(FatalError);
            }
        }

        if (facei < m.nInternalFaces() || facei >= m.nFaces())
        {
            newPatch = -1;
        }
        else
        {
            newPatch = patchIndices[bm.whichPatch(facei)];
        }

        newOwn = own[facei];
        if (newOwn > -1) newOwn = renumberCells[newOwn];

        if (newPatch > -1)
        {
            newNei = -1;
        }
        else
        {
            newNei = nei[facei];
            newNei = renumberCells[newNei];
        }

        renumberFaces[facei] = meshMod_.addFace
        (
            newFace,
            newOwn,
            newNei,
            -1,
            false,
            newPatch
        );

        const labelList zones(fz.whichZones(facei));
        forAll(zones, zonei)
        {
            const faceZone& fzi = fz[zones[zonei]];
            const bool flip = fzi.flipMap()[fzi.localIndex(facei)];

            faceZonesAddedFaces_[faceZoneIndices[zones[zonei]]]
           .insert(renumberFaces[facei], flip);
        }
    }
}


void Foam::mergePolyMesh::merge()
{
    if (debug)
    {
        Info<< "patch names: " << patchNames_ << nl
            << "patch dicts: " << patchDicts_ << nl
            << "point zone names: " << pointZoneNames_ << nl
            << "face zone names: " << faceZoneNames_ << nl
            << "cell zone names: " << cellZoneNames_ << endl;
    }

    // Add the patches if necessary
    if (patchNames_.size() != mesh_.boundaryMesh().size())
    {
        if (debug)
        {
            Info<< "Copying old patches" << endl;
        }

        List<polyPatch*> newPatches(patchNames_.size());

        const polyBoundaryMesh& oldPatches = mesh_.boundaryMesh();

        // Note.  Re-using counter in two for loops
        label patchi = 0;

        for (patchi = 0; patchi < oldPatches.size(); patchi++)
        {
            newPatches[patchi] = oldPatches[patchi].clone(oldPatches).ptr();
        }

        if (debug)
        {
            Info<< "Adding new patches. " << endl;
        }

        label endOfLastPatch =
            patchi == 0
          ? 0
          : oldPatches[patchi - 1].start() + oldPatches[patchi - 1].size();

        for (; patchi < patchNames_.size(); patchi++)
        {
            // Add a patch
            dictionary dict(patchDicts_[patchi]);
            dict.set("nFaces", 0);
            dict.set("startFace", endOfLastPatch);

            newPatches[patchi] =
            (
                polyPatch::New
                (
                    patchNames_[patchi],
                    dict,
                    patchi,
                    oldPatches
                ).ptr()
            );
        }

        mesh_.removeBoundary();
        mesh_.addPatches(newPatches);
    }

    // Add the zones if necessary
    if (pointZoneNames_.size() > mesh_.pointZones().size())
    {
        if (debug)
        {
            Info<< "Adding new pointZones. " << endl;
        }

        label nZones = mesh_.pointZones().size();

        mesh_.pointZones().setSize(pointZoneNames_.size());

        for (label zonei = nZones; zonei < pointZoneNames_.size(); zonei++)
        {
            mesh_.pointZones().set
            (
                zonei,
                new pointZone
                (
                    pointZoneNames_[zonei],
                    labelList(),
                    mesh_.pointZones()
                )
            );
        }
    }

    if (cellZoneNames_.size() > mesh_.cellZones().size())
    {
        if (debug)
        {
            Info<< "Adding new cellZones. " << endl;
        }

        label nZones = mesh_.cellZones().size();

        mesh_.cellZones().setSize(cellZoneNames_.size());

        for (label zonei = nZones; zonei < cellZoneNames_.size(); zonei++)
        {
            mesh_.cellZones().set
            (
                zonei,
                new cellZone
                (
                    cellZoneNames_[zonei],
                    labelList(),
                    mesh_.cellZones()
                )
            );
        }
    }

    if (faceZoneNames_.size() > mesh_.faceZones().size())
    {
        if (debug)
        {
            Info<< "Adding new faceZones. " << endl;
        }

        label nZones = mesh_.faceZones().size();

        mesh_.faceZones().setSize(faceZoneNames_.size());

        for (label zonei = nZones; zonei < faceZoneNames_.size(); zonei++)
        {
            mesh_.faceZones().set
            (
                zonei,
                new faceZone
                (
                    faceZoneNames_[zonei],
                    labelList(),
                    boolList(),
                    mesh_.faceZones()
                )
            );
        }
    }

    // Change mesh
    autoPtr<polyTopoChangeMap> map(meshMod_.changeMesh(mesh_));

    // Add the new points to the pointZones in the merged mesh
    mesh_.pointZones().insert(pointZonesAddedPoints_);

    // Add the new faces to the faceZones in the merged mesh
    mesh_.faceZones().insert(faceZonesAddedFaces_);

    // Add the new cells to the cellZones in the merged mesh
    mesh_.cellZones().insert(cellZonesAddedCells_);

    mesh_.topoChange(map);

    // Clear topo change for the next operation
    meshMod_.clear();
}


// ************************************************************************* //
