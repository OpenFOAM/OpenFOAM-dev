/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    boundary faces
    - use pointCells when searching for connectivity
    - initialize the cell connectivity with '-1'
    - find both cell faces corresponding to the baffles and mark them
      to prevent a connection
    - standard connectivity checks

    - added baffle and monitoring support

\*---------------------------------------------------------------------------*/

#include "meshReader.H"
#include "Time.H"
#include "polyPatch.H"
#include "emptyPolyPatch.H"
#include "preservePatchTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::meshReader::addPolyBoundaryFace
(
    const label cellId,
    const label cellFaceId,
    const label nCreatedFaces
)
{
#ifdef DEBUG_BOUNDARY
    Info<< nCreatedFaces
        << " add bnd for cell " << cellId
        << " face " << cellFaceId
        << " (original cell " << origCellId_[cellId] << ")"
        << endl;
#endif

    // standard case: volume cells
    const face& thisFace = cellFaces_[cellId][cellFaceId];

    // Debugging
    if (cellPolys_[cellId][cellFaceId] > nInternalFaces_)
    {
        Info<< "meshReader::createPolyBoundary(): "
            << "Problem with face: " << thisFace << endl
            << "Probably multiple definitions "
            << "of a single boundary face." << endl
            << endl;
    }
    else if (cellPolys_[cellId][cellFaceId] >= 0)
    {
        Info<< "meshReader::createPolyBoundary(): "
            << "Problem with face: " << thisFace << endl
            << "Probably trying to define a boundary face "
            << "on a previously matched internal face." << endl
            << "Internal face: "
            << meshFaces_[cellPolys_[cellId][cellFaceId]]
            << endl;
    }

    meshFaces_[nCreatedFaces] = thisFace;
    cellPolys_[cellId][cellFaceId] = nCreatedFaces;
}


void Foam::meshReader::addPolyBoundaryFace
(
    const cellFaceIdentifier& identifier,
    const label nCreatedFaces
)
{
    addPolyBoundaryFace(identifier.cell, identifier.face, nCreatedFaces);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::meshReader::createPolyBoundary()
{
    label nBoundaryFaces = 0;
    label nMissingFaces = 0;
    label nInterfaces = 0;

    const faceListList& cFaces = cellFaces();

    // determine number of non-patched faces:
    forAll(cellPolys_, cellI)
    {
        cell& curCell = cellPolys_[cellI];

        forAll(curCell, fI)
        {
            if (curCell[fI] < 0)
            {
                nMissingFaces++;
            }
        }
    }

    forAll(boundaryIds_, patchI)
    {
        nBoundaryFaces += boundaryIds_[patchI].size();
    }

    Info<< nl
        << "There are " << nMissingFaces
        << " faces to be patched and " << nBoundaryFaces
        << " specified - collect missed boundaries to final patch" << endl;

    patchStarts_.setSize(boundaryIds_.size());
    patchSizes_.setSize(boundaryIds_.size());

    label nCreatedFaces = nInternalFaces_;
    label baffleOffset  = cFaces.size();
    interfaces_.setSize(baffleIds_.size());
    nBoundaryFaces = 0;

    forAll(boundaryIds_, patchI)
    {
        const List<cellFaceIdentifier>& idList = boundaryIds_[patchI];

        patchStarts_[patchI] = nCreatedFaces;

        // write each baffle side separately
        if (patchPhysicalTypes_[patchI] == "baffle")
        {
            label count = 0;

            for (label side = 0; side < 2; ++side)
            {
                label position = nInterfaces;

                forAll(idList, bndI)
                {
                    label baffleI = idList[bndI].cell - baffleOffset;

                    if
                    (
                        baffleI >= 0
                     && baffleI < baffleFaces_.size()
                     && baffleIds_[baffleI].size()
                    )
                    {
                        addPolyBoundaryFace
                        (
                            baffleIds_[baffleI][side],
                            nCreatedFaces
                        );

                        // remove applied boundaries (2nd pass)
                        if (side == 1)
                        {
                            baffleIds_[baffleI].clear();
                        }

                        interfaces_[position][side] = nCreatedFaces;

                        nBoundaryFaces++;
                        nCreatedFaces++;
                        position++;
                        count++;
                    }
                }
            }

            nInterfaces += (count - (count % 2)) / 2;
        }
        else if (patchPhysicalTypes_[patchI] == "monitoring")
        {
            // translate the "monitoring" pseudo-boundaries to face sets
            List<label> monitoring(idList.size());

            label monitorI = 0;
            forAll(idList, bndI)
            {
                label cellId = idList[bndI].cell;
                label faceId = idList[bndI].face;

                // standard case: volume cells
                if (cellId < baffleOffset)
                {
                    label faceNr = cellPolys_[cellId][faceId];
                    if (faceNr >= 0)
                    {
                        monitoring[monitorI++] = faceNr;
                    }
                }
            }

            monitoringSets_.insert(patchNames_[patchI], monitoring);
        }
        else
        {
            forAll(idList, bndI)
            {
                // standard case: volume cells
                if (idList[bndI].cell < baffleOffset)
                {
                    addPolyBoundaryFace
                    (
                        idList[bndI],
                        nCreatedFaces
                    );

                    nBoundaryFaces++;
                    nCreatedFaces++;
                }
            }
        }

        patchSizes_[patchI] = nCreatedFaces - patchStarts_[patchI];
    }

    // add in missing faces
    Info<< "Missing faces added to patch after face "
        << nCreatedFaces << ":" <<endl;
    nMissingFaces = 0;

    // look for baffles first - keep them together at the start of the patch
    for (label side = 0; side < 2; ++side)
    {
        label position = nInterfaces;

        forAll(baffleIds_, baffleI)
        {
            if (baffleIds_[baffleI].size())
            {
                // add each side for each baffle
                addPolyBoundaryFace
                (
                    baffleIds_[baffleI][side],
                    nCreatedFaces
                );

                interfaces_[position][side] = nCreatedFaces;

                // remove applied boundaries (2nd pass)
                if (side == 1)
                {
                    baffleIds_[baffleI].clear();
                }

                nMissingFaces++;
                nCreatedFaces++;
                position++;
            }
        }
    }

    nInterfaces += (nMissingFaces - (nMissingFaces % 2)) / 2;

    // scan for any other missing faces
    forAll(cellPolys_, cellI)
    {
        const labelList& curFaces = cellPolys_[cellI];

        forAll(curFaces, cellFaceI)
        {
            if (curFaces[cellFaceI] < 0)
            {
                // just report the first few
                if (nMissingFaces < 4)
                {
                    const face& thisFace = cFaces[cellI][cellFaceI];

                    Info<< "  cell " << cellI << " face " << cellFaceI
                        << " (original cell " << origCellId_[cellI] << ")"
                        << " face: " << thisFace
                        << endl;
                }
                else if (nMissingFaces == 5)
                {
                    Info<< "  ..." << nl << endl;
                }

                addPolyBoundaryFace(cellI, cellFaceI, nCreatedFaces);
                nMissingFaces++;
                nCreatedFaces++;
            }
        }
    }

    Info<< "Added " << nMissingFaces << " unmatched faces" << endl;

    // Add missing faces to last patch ('Default_Empty' etc.)
    if (nMissingFaces > 0)
    {
        patchSizes_.last() = nMissingFaces;
    }


    // reset the size of the face list
    meshFaces_.setSize(nCreatedFaces);

    // check the mesh for face mismatch
    // (faces addressed once or more than twice)
    labelList markupFaces(meshFaces_.size(), 0);

    forAll(cellPolys_, cellI)
    {
        const labelList& curFaces = cellPolys_[cellI];

        forAll(curFaces, faceI)
        {
            markupFaces[curFaces[faceI]]++;
        }
    }

    for (label i = nInternalFaces_; i < markupFaces.size(); i++)
    {
        markupFaces[i]++;
    }

    label nProblemFaces = 0;

    forAll(markupFaces, faceI)
    {
        if (markupFaces[faceI] != 2)
        {
            const face& problemFace = meshFaces_[faceI];

            Info<< "meshReader::createPolyBoundary() : "
                << "problem with face " << faceI << ": addressed "
                << markupFaces[faceI] << " times (should be 2!). Face: "
                << problemFace << endl;

            nProblemFaces++;
        }
    }

    if (nProblemFaces > 0)
    {
        Info<< "Number of incorrectly matched faces: "
            << nProblemFaces << endl;
    }

    // adjust for missing members
    if (nInterfaces < interfaces_.size())
    {
        interfaces_.setSize(nInterfaces);
    }

    Info<< "Number of boundary faces: " << nBoundaryFaces << nl
        << "Total number of faces: " << nCreatedFaces << nl
        << "Number of interfaces: " << nInterfaces << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::List<Foam::polyPatch*>
Foam::meshReader::polyBoundaryPatches(const polyMesh& mesh)
{
    label nUsed = 0, nEmpty = 0;
    label nPatches = patchStarts_.size();

    // avoid empty patches - move to the end of the lists and truncate
    labelList oldToNew = identity(nPatches);
    forAll(patchSizes_, patchI)
    {
        if (patchSizes_[patchI] > 0)
        {
            oldToNew[patchI] = nUsed++;
        }
        else
        {
            nEmpty++;
            oldToNew[patchI] = nPatches - nEmpty;
        }
    }

    nPatches = nUsed;

    if (nEmpty)
    {
        Info<< "Removing " << nEmpty << " empty patches" << endl;

        inplaceReorder(oldToNew, patchTypes_);
        inplaceReorder(oldToNew, patchNames_);
        inplaceReorder(oldToNew, patchStarts_);
        inplaceReorder(oldToNew, patchSizes_);
    }

    patchTypes_.setSize(nPatches);
    patchNames_.setSize(nPatches);
    patchStarts_.setSize(nPatches);
    patchSizes_.setSize(nPatches);


    List<polyPatch*> p(nPatches);

    // All patch dictionaries
    PtrList<dictionary> patchDicts(patchNames_.size());
    // Default boundary patch types
    word defaultFacesType(emptyPolyPatch::typeName);

    // we could consider dropping this entirely
    preservePatchTypes
    (
        mesh,
        mesh.instance(),
        mesh.meshDir(),
        patchNames_,
        patchDicts,
        "defaultFaces",
        defaultFacesType
    );
    forAll(patchDicts, patchI)
    {
        if (!patchDicts.set(patchI))
        {
            patchDicts.set(patchI, new dictionary());
        }
        dictionary& patchDict = patchDicts[patchI];

        // add but not overwrite type
        patchDict.add("type", patchTypes_[patchI], false);
        if (patchPhysicalTypes_.size() && patchPhysicalTypes_[patchI].size())
        {
            patchDict.add("startFace", patchPhysicalTypes_[patchI], false);
        }

        // overwrite sizes and start
        patchDict.add("nFaces", patchSizes_[patchI], true);
        patchDict.add("startFace", patchStarts_[patchI], true);
    }


    forAll(patchStarts_, patchI)
    {
        p[patchI] = polyPatch::New
        (
            patchNames_[patchI],
            patchDicts[patchI],
            patchI,
            mesh.boundaryMesh()
        ).ptr();
    }

    return p;
}


// ************************************************************************* //
