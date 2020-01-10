/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "STARCDMeshReader.H"
#include "mergedCyclicPolyPatch.H"
#include "emptyPolyPatch.H"
#include "wallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "cellModeller.H"
#include "ListOps.H"
#include "IFstream.H"
#include "IOMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::meshReaders::STARCD::defaultBoundaryName =
    "Default_Boundary_Region";

const char* const Foam::meshReaders::STARCD::defaultSolidBoundaryName =
    "Default_Boundary_Solid";

bool Foam::meshReaders::STARCD::keepSolids = false;

const int Foam::meshReaders::STARCD::starToFoamFaceAddr[4][6] =
{
    { 4, 5, 2, 3, 0, 1 },     // 11 = pro-STAR hex
    { 0, 1, 4, -1, 2, 3 },    // 12 = pro-STAR prism
    { 3, -1, 2, -1, 1, 0 },   // 13 = pro-STAR tetra
    { 0, -1, 4, 2, 1, 3 }     // 14 = pro-STAR pyramid
};


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::meshReaders::STARCD::readToNewline(IFstream& is)
{
    char ch = '\n';
    do
    {
        (is).get(ch);
    }
    while ((is) && ch != '\n');
}


bool Foam::meshReaders::STARCD::readHeader(IFstream& is, word fileSignature)
{
    if (!is.good())
    {
        FatalErrorInFunction
            << abort(FatalError);
    }

    word header;
    label majorVersion;

    is >> header;
    is >> majorVersion;

    // skip the rest of the line
    readToNewline(is);

    // add other checks ...
    if (header != fileSignature)
    {
        Info<< "header mismatch " << fileSignature << "  " << is.name()
            << endl;
    }

    return true;
}


void Foam::meshReaders::STARCD::readAux(const objectRegistry& registry)
{
    boundaryRegion_.readDict(registry);
    cellTable_.readDict(registry);
}


// read in the points from the .vrt file
//
/*---------------------------------------------------------------------------*\
Line 1:
  PROSTAR_VERTEX [newline]

Line 2:
  <version> 0 0 0 0 0 0 0 [newline]

Body:
  <vertexId>  <x>  <y>  <z> [newline]

\*---------------------------------------------------------------------------*/
void Foam::meshReaders::STARCD::readPoints
(
    const fileName& inputName,
    const scalar scaleFactor
)
{
    const word fileSignature = "PROSTAR_VERTEX";
    label nPoints = 0, maxId = 0;

    // Pass 1:
    // get # points and maximum vertex label
    {
        IFstream is(inputName);
        readHeader(is, fileSignature);

        label lineLabel;
        scalar x, y, z;

        while ((is >> lineLabel).good())
        {
            nPoints++;
            maxId = max(maxId, lineLabel);
            is >> x >> y >> z;
        }
    }

    Info<< "Number of points  = " << nPoints << endl;

    // set sizes and reset to invalid values

    points_.setSize(nPoints);
    mapToFoamPointId_.setSize(maxId+1);

    //- Original Point number for a given vertex
    // might need again in the future
    ////     labelList origPointId(nPoints);
    ////     origPointId = -1;

    mapToFoamPointId_ = -1;

    // Pass 2:
    // construct pointList and conversion table
    // from Star vertex numbers to Foam point labels
    if (nPoints > 0)
    {
        IFstream is(inputName);
        readHeader(is, fileSignature);

        label lineLabel;

        label pointi = 0;
        while ((is >> lineLabel).good())
        {
            is  >> points_[pointi].x()
                >> points_[pointi].y()
                >> points_[pointi].z();

            // might need again in the future
            ////  origPointId[pointi] = lineLabel;
            mapToFoamPointId_[lineLabel] = pointi;
            pointi++;
        }

        if (nPoints > pointi)
        {
            nPoints = pointi;
            points_.setSize(nPoints);
            // might need again in the future
            //// origPointId.setSize(nPoints);
        }

        if (scaleFactor > 1.0 + small || scaleFactor < 1.0 - small)
        {
            points_ *= scaleFactor;
        }
    }
    else
    {
        FatalErrorInFunction
            << "no points in file " << inputName
            << abort(FatalError);
    }

}


// read in the cells from the .cel file
//
/*---------------------------------------------------------------------------*\
Line 1:
  PROSTAR_CELL [newline]

Line 2:
  <version> 0 0 0 0 0 0 0 [newline]

Body:
  <cellId>  <shapeId>  <nLabels>  <cellTableId>  <typeId> [newline]
  <cellId>  <int1> .. <int8>
  <cellId>  <int9> .. <int16>

 with shapeId:
 *   1 = point
 *   2 = line
 *   3 = shell
 *  11 = hexa
 *  12 = prism
 *  13 = tetra
 *  14 = pyramid
 * 255 = polyhedron

 with typeId
 *   1 = fluid
 *   2 = solid
 *   3 = baffle
 *   4 = shell
 *   5 = line
 *   6 = point

For primitive cell shapes, the number of vertices will never exceed 8 (hexa)
and corresponds to <nLabels>.
For polyhedral, <nLabels> includes an index table comprising beg/end pairs
for each cell face.

Strictly speaking, we only need the cellModeller for adding boundaries.
\*---------------------------------------------------------------------------*/

void Foam::meshReaders::STARCD::readCells(const fileName& inputName)
{
    const word fileSignature = "PROSTAR_CELL";
    label nFluids = 0, nSolids = 0, nBaffles = 0, nShells = 0;
    label maxId = 0;

    bool unknownVertices = false;


    // Pass 1:
    // count nFluids, nSolids, nBaffle, nShell and maxId
    // also see if polyhedral cells were used
    {
        IFstream is(inputName);
        readHeader(is, fileSignature);

        label lineLabel, shapeId, nLabels, cellTableId, typeId;

        while ((is >> lineLabel).good())
        {
            label starCellId = lineLabel;
            is  >> shapeId
                >> nLabels
                >> cellTableId
                >> typeId;

            // skip the rest of the line
            readToNewline(is);

            // max 8 indices per line
            while (nLabels > 0)
            {
                readToNewline(is);
                nLabels -= 8;
            }

            if (typeId == starcdFluidType)
            {
                nFluids++;
                maxId = max(maxId, starCellId);

                if (!cellTable_.found(cellTableId))
                {
                    cellTable_.setName(cellTableId);
                    cellTable_.setMaterial(cellTableId, "fluid");
                }
            }
            else if (typeId == starcdSolidType)
            {
                nSolids++;
                if (keepSolids)
                {
                    maxId = max(maxId, starCellId);
                }

                if (!cellTable_.found(cellTableId))
                {
                    cellTable_.setName(cellTableId);
                    cellTable_.setMaterial(cellTableId, "solid");
                }

            }
            else if (typeId == starcdBaffleType)
            {
                // baffles have no cellTable entry
                nBaffles++;
                maxId = max(maxId, starCellId);
            }
            else if (typeId == starcdShellType)
            {
                nShells++;
                if (!cellTable_.found(cellTableId))
                {
                    cellTable_.setName(cellTableId);
                    cellTable_.setMaterial(cellTableId, "shell");
                }
            }

        }
    }

    Info<< "Number of fluids  = " << nFluids << nl
        << "Number of baffles = " << nBaffles << nl;
    if (keepSolids)
    {
        Info<< "Number of solids  = " << nSolids << nl;
    }
    else
    {
        Info<< "Ignored   solids  = " << nSolids << nl;
    }
    Info<< "Ignored   shells  = " << nShells << endl;


    label nCells;
    if (keepSolids)
    {
        nCells = nFluids + nSolids;
    }
    else
    {
        nCells = nFluids;
    }

    cellFaces_.setSize(nCells);
    cellShapes_.setSize(nCells);
    cellTableId_.setSize(nCells);

    // information for the interfaces
    baffleFaces_.setSize(nBaffles);

    // extra space for baffles
    origCellId_.setSize(nCells + nBaffles);
    mapToFoamCellId_.setSize(maxId+1);
    mapToFoamCellId_ = -1;


    // avoid undefined shapes for polyhedra
    cellShape genericShape(*unknownModel, labelList(0));

    // Pass 2:
    // construct cellFaces_ and possibly cellShapes_
    if (nCells <= 0)
    {
        FatalErrorInFunction
            << "no cells in file " << inputName
            << abort(FatalError);
    }
    else
    {
        IFstream is(inputName);
        readHeader(is, fileSignature);

        labelList starLabels(64);
        label lineLabel, shapeId, nLabels, cellTableId, typeId;

        label celli = 0;
        label baffleI = 0;

        while ((is >> lineLabel).good())
        {
            label starCellId = lineLabel;
            is  >> shapeId
                >> nLabels
                >> cellTableId
                >> typeId;

            if (nLabels > starLabels.size())
            {
                starLabels.setSize(nLabels);
            }
            starLabels = -1;

            // read indices - max 8 per line
            for (label i = 0; i < nLabels; ++i)
            {
                if ((i % 8) == 0)
                {
                    is >> lineLabel;
                }
                is >> starLabels[i];
            }

            // skip solid cells
            if (typeId == starcdSolidType && !keepSolids)
            {
                continue;
            }

            // determine the foam cell shape
            const cellModel* curModelPtr = nullptr;

            // fluid/solid cells
            switch (shapeId)
            {
                case starcdHex:
                    curModelPtr = hexModel;
                    break;
                case starcdPrism:
                    curModelPtr = prismModel;
                    break;
                case starcdTet:
                    curModelPtr = tetModel;
                    break;
                case starcdPyr:
                    curModelPtr = pyrModel;
                    break;
            }

            if (curModelPtr)
            {
                // primitive cell - use shapes

                // convert orig vertex Id to point label
                bool isBad = false;
                for (label i=0; i < nLabels; ++i)
                {
                    label pointId = mapToFoamPointId_[starLabels[i]];
                    if (pointId < 0)
                    {
                        Info<< "Cells inconsistent with vertex file. "
                            << "Star vertex " << starLabels[i]
                            << " does not exist" << endl;
                        isBad = true;
                        unknownVertices = true;
                    }
                    starLabels[i] = pointId;
                }

                if (isBad)
                {
                    continue;
                }

                // record original cell number and lookup
                origCellId_[celli] = starCellId;
                mapToFoamCellId_[starCellId] = celli;

                cellTableId_[celli] = cellTableId;
                cellShapes_[celli] = cellShape
                (
                    *curModelPtr,
                    SubList<label>(starLabels, nLabels)
                );

                cellFaces_[celli] = cellShapes_[celli].faces();
                celli++;
            }
            else if (shapeId == starcdPoly)
            {
                // polyhedral cell
                label nFaces = starLabels[0] - 1;

                // convert orig vertex id to point label
                // start with offset (skip the index table)
                bool isBad = false;
                for (label i=starLabels[0]; i < nLabels; ++i)
                {
                    label pointId = mapToFoamPointId_[starLabels[i]];
                    if (pointId < 0)
                    {
                        Info<< "Cells inconsistent with vertex file. "
                            << "Star vertex " << starLabels[i]
                            << " does not exist" << endl;
                        isBad = true;
                        unknownVertices = true;
                    }
                    starLabels[i] = pointId;
                }

                if (isBad)
                {
                    continue;
                }

                // traverse beg/end indices
                faceList faces(nFaces);
                label facei = 0;
                for (label i=0; i < nFaces; ++i)
                {
                    label beg = starLabels[i];
                    label n   = starLabels[i+1] - beg;

                    face f
                    (
                        SubList<label>(starLabels, n, beg)
                    );

                    f.collapse();

                    // valid faces only
                    if (f.size() >= 3)
                    {
                        faces[facei++] = f;
                    }
                }

                if (nFaces > facei)
                {
                    Info<< "star cell " << starCellId << " has "
                        << (nFaces - facei)
                        << " empty faces - could cause boundary "
                        << "addressing problems"
                        << endl;

                    nFaces = facei;
                    faces.setSize(nFaces);
                }

                if (nFaces < 4)
                {
                    FatalErrorInFunction
                        << "star cell " << starCellId << " has " << nFaces
                        << abort(FatalError);
                }

                // record original cell number and lookup
                origCellId_[celli] = starCellId;
                mapToFoamCellId_[starCellId] = celli;

                cellTableId_[celli] = cellTableId;
                cellShapes_[celli]  = genericShape;
                cellFaces_[celli]   = faces;
                celli++;
            }
            else if (typeId == starcdBaffleType)
            {
                // baffles

                // convert orig vertex id to point label
                bool isBad = false;
                for (label i=0; i < nLabels; ++i)
                {
                    label pointId = mapToFoamPointId_[starLabels[i]];
                    if (pointId < 0)
                    {
                        Info<< "Baffles inconsistent with vertex file. "
                            << "Star vertex " << starLabels[i]
                            << " does not exist" << endl;
                        isBad = true;
                        unknownVertices = true;
                    }
                    starLabels[i] = pointId;
                }

                if (isBad)
                {
                    continue;
                }


                face f
                (
                    SubList<label>(starLabels, nLabels)
                );

                f.collapse();

                // valid faces only
                if (f.size() >= 3)
                {
                    baffleFaces_[baffleI] = f;
                    // insert lookup addressing in normal list
                    mapToFoamCellId_[starCellId]  = nCells + baffleI;
                    origCellId_[nCells + baffleI] = starCellId;
                    baffleI++;
                }
            }
        }

        baffleFaces_.setSize(baffleI);
    }

    if (unknownVertices)
    {
        FatalErrorInFunction
            << "cells with unknown vertices"
            << abort(FatalError);
    }

    // truncate lists

#ifdef DEBUG_READING
    Info<< "CELLS READ" << endl;
#endif

    // cleanup
    mapToFoamPointId_.clear();
}


// read in the boundaries from the .bnd file
//
/*---------------------------------------------------------------------------*\
Line 1:
  PROSTAR_BOUNDARY [newline]

Line 2:
  <version> 0 0 0 0 0 0 0 [newline]

Body:
  <boundId>  <cellId>  <cellFace>  <regionId>  0  <boundaryType> [newline]

where boundaryType is truncated to 4 characters from one of the following:
INLET
PRESSURE
OUTLET
BAFFLE
etc,
\*---------------------------------------------------------------------------*/

void Foam::meshReaders::STARCD::readBoundary(const fileName& inputName)
{
    const word fileSignature = "PROSTAR_BOUNDARY";
    label nPatches = 0, nFaces = 0, nBafflePatches = 0, maxId = 0;
    label lineLabel, starCellId, cellFaceId, starRegion, configNumber;
    word patchType;

    labelList mapToFoamPatchId(1000, label(-1));
    labelList nPatchFaces(1000, label(0));
    labelList origRegion(1000, label(0));
    patchTypes_.setSize(1000);

    // this is what we seem to need
    // these MUST correspond to starToFoamFaceAddr
    //
    Map<label> faceLookupIndex;

    faceLookupIndex.insert(hexModel->index(), 0);
    faceLookupIndex.insert(prismModel->index(), 1);
    faceLookupIndex.insert(tetModel->index(), 2);
    faceLookupIndex.insert(pyrModel->index(), 3);

    // Pass 1:
    // collect
    // no. of faces (nFaces), no. of patches (nPatches)
    // and for each of these patches the number of faces
    // (nPatchFaces[patchLabel])
    //
    // and a conversion table from Star regions to (Foam) patchLabels
    //
    // additionally note the no. of baffle patches (nBafflePatches)
    // so that we sort these to the end of the patch list
    // - this makes it easier to transfer them to an adjacent patch if reqd
    {
        IFstream is(inputName);

        if (is.good())
        {
            readHeader(is, fileSignature);

            while ((is >> lineLabel).good())
            {
                nFaces++;
                is  >> starCellId
                    >> cellFaceId
                    >> starRegion
                    >> configNumber
                    >> patchType;

                // Build translation table to convert star patch to foam patch
                label patchLabel = mapToFoamPatchId[starRegion];
                if (patchLabel == -1)
                {
                    patchLabel = nPatches;
                    mapToFoamPatchId[starRegion] = patchLabel;
                    origRegion[patchLabel] = starRegion;
                    patchTypes_[patchLabel] = patchType;

                    maxId = max(maxId, starRegion);

                    // should actually be case-insensitive
                    if (patchType == "BAFF")
                    {
                        nBafflePatches++;
                    }
                    nPatches++;
                }

                nPatchFaces[patchLabel]++;
            }

            if (nPatches == 0)
            {
                Info<< "No boundary faces in file " << inputName << endl;
            }
        }
        else
        {
            Info<< "Could not read boundary file " << inputName << endl;
        }
    }

    // keep empty patch region in reserve
    nPatches++;
    Info<< "Number of patches = " << nPatches
        << " (including extra for missing)" << endl;

    // resize
    origRegion.setSize(nPatches);
    patchTypes_.setSize(nPatches);
    patchNames_.setSize(nPatches);
    nPatchFaces.setSize(nPatches);

    // add our empty patch
    origRegion[nPatches-1] = 0;
    nPatchFaces[nPatches-1] = 0;
    patchTypes_[nPatches-1] = "none";

    // create names
    // - use 'Label' entry from "constant/boundaryRegion" dictionary
    forAll(patchTypes_, patchi)
    {
        bool foundName = false, foundType = false;

        Map<dictionary>::const_iterator
            iter = boundaryRegion_.find(origRegion[patchi]);

        if
        (
            iter != boundaryRegion_.end()
        )
        {
            foundType = iter().readIfPresent
            (
                "BoundaryType",
                patchTypes_[patchi]
            );

            foundName = iter().readIfPresent
            (
                "Label",
                patchNames_[patchi]
            );
        }

        // consistent names, in long form and in lowercase
        if (!foundType)
        {
            // transform
            forAllIter(string, patchTypes_[patchi], i)
            {
                *i = tolower(*i);
            }

            if (patchTypes_[patchi] == "symp")
            {
                patchTypes_[patchi] = "symplane";
            }
            else if (patchTypes_[patchi] == "cycl")
            {
                patchTypes_[patchi] = "cyclic";
            }
            else if (patchTypes_[patchi] == "baff")
            {
                patchTypes_[patchi] = "baffle";
            }
            else if (patchTypes_[patchi] == "moni")
            {
                patchTypes_[patchi] = "monitoring";
            }
        }

        // create a name if needed
        if (!foundName)
        {
            patchNames_[patchi] =
                patchTypes_[patchi] + "_" + name(origRegion[patchi]);
        }
    }

    // enforce name "Default_Boundary_Region"
    patchNames_[nPatches-1] = defaultBoundaryName;

    // sort according to ascending region numbers, but leave
    // Default_Boundary_Region as the final patch
    {
        labelList sortedIndices;
        sortedOrder(SubList<label>(origRegion, nPatches-1), sortedIndices);

        labelList oldToNew = identity(nPatches);
        forAll(sortedIndices, i)
        {
            oldToNew[sortedIndices[i]] = i;
        }

        inplaceReorder(oldToNew, origRegion);
        inplaceReorder(oldToNew, patchTypes_);
        inplaceReorder(oldToNew, patchNames_);
        inplaceReorder(oldToNew, nPatchFaces);
    }

    // re-sort to have baffles near the end
    nBafflePatches = 1;
    if (nBafflePatches)
    {
        labelList oldToNew = identity(nPatches);
        label newIndex = 0;
        label baffleIndex = (nPatches-1 - nBafflePatches);

        for (label i=0; i < oldToNew.size()-1; ++i)
        {
            if (patchTypes_[i] == "baffle")
            {
                oldToNew[i] = baffleIndex++;
            }
            else
            {
                oldToNew[i] = newIndex++;
            }
        }

        inplaceReorder(oldToNew, origRegion);
        inplaceReorder(oldToNew, patchTypes_);
        inplaceReorder(oldToNew, patchNames_);
        inplaceReorder(oldToNew, nPatchFaces);
    }

    mapToFoamPatchId.setSize(maxId+1, -1);
    forAll(origRegion, patchi)
    {
        mapToFoamPatchId[origRegion[patchi]] = patchi;
    }

    boundaryIds_.setSize(nPatches);
    forAll(boundaryIds_, patchi)
    {
        boundaryIds_[patchi].setSize(nPatchFaces[patchi]);
        nPatchFaces[patchi] = 0;
    }


    // Pass 2:
    //
    if (nPatches > 1 && mapToFoamCellId_.size() > 1)
    {
        IFstream is(inputName);
        readHeader(is, fileSignature);

        while ((is >> lineLabel).good())
        {
            is
                >> starCellId
                >> cellFaceId
                >> starRegion
                >> configNumber
                >> patchType;

            label patchi = mapToFoamPatchId[starRegion];

            // zero-based indexing
            cellFaceId--;

            label cellId = -1;

            // convert to FOAM cell number
            if (starCellId < mapToFoamCellId_.size())
            {
                cellId = mapToFoamCellId_[starCellId];
            }

            if (cellId < 0)
            {
                Info
                    << "Boundaries inconsistent with cell file. "
                    << "Star cell " << starCellId << " does not exist"
                    << endl;
            }
            else
            {
                // restrict lookup to volume cells (no baffles)
                if (cellId < cellShapes_.size())
                {
                    label index = cellShapes_[cellId].model().index();
                    if (faceLookupIndex.found(index))
                    {
                        index = faceLookupIndex[index];
                        cellFaceId = starToFoamFaceAddr[index][cellFaceId];
                    }
                }
                else
                {
                    // we currently use cellId >= nCells to tag baffles,
                    // we can also use a negative face number
                    cellFaceId = -1;
                }

                boundaryIds_[patchi][nPatchFaces[patchi]] =
                    cellFaceIdentifier(cellId, cellFaceId);

#ifdef DEBUG_BOUNDARY
                Info<< "bnd " << cellId << " " << cellFaceId << endl;
#endif
                // increment counter of faces in current patch
                nPatchFaces[patchi]++;
            }
        }
    }

    // retain original information in patchPhysicalTypes_ - overwrite latter
    patchPhysicalTypes_.setSize(patchTypes_.size());


    forAll(boundaryIds_, patchi)
    {
        // resize - avoid invalid boundaries
        if (nPatchFaces[patchi] < boundaryIds_[patchi].size())
        {
            boundaryIds_[patchi].setSize(nPatchFaces[patchi]);
        }

        word origType = patchTypes_[patchi];
        patchPhysicalTypes_[patchi] = origType;

        if (origType == "symplane")
        {
            patchTypes_[patchi] = symmetryPolyPatch::typeName;
            patchPhysicalTypes_[patchi] = patchTypes_[patchi];
        }
        else if (origType == "wall")
        {
            patchTypes_[patchi] = wallPolyPatch::typeName;
            patchPhysicalTypes_[patchi] = patchTypes_[patchi];
        }
        else if (origType == "cyclic")
        {
            patchTypes_[patchi] = mergedCyclicPolyPatch::typeName;
            patchPhysicalTypes_[patchi] = patchTypes_[patchi];
        }
        else if (origType == "baffle")
        {
            patchTypes_[patchi] = mergedCyclicPolyPatch::typeName;
            patchPhysicalTypes_[patchi] = "baffle";
        }
        else
        {
            patchTypes_[patchi] = polyPatch::typeName;
        }

        Info<< "patch " << patchi
            << " (region " << origRegion[patchi]
            << ": " << origType << ") type: '" << patchTypes_[patchi]
            << "' name: " << patchNames_[patchi] << endl;
    }

    // cleanup
    mapToFoamCellId_.clear();
    cellShapes_.clear();
}


//
// remove unused points
//
void Foam::meshReaders::STARCD::cullPoints()
{
    label nPoints = points_.size();
    labelList oldToNew(nPoints, -1);

    // loop through cell faces and note which points are being used
    forAll(cellFaces_, celli)
    {
        const faceList& faces = cellFaces_[celli];
        forAll(faces, i)
        {
            const labelList& labels = faces[i];
            forAll(labels, j)
            {
                oldToNew[labels[j]]++;
            }
        }
    }

    // the new ordering and the count of unused points
    label pointi = 0;
    forAll(oldToNew, i)
    {
        if (oldToNew[i] >= 0)
        {
            oldToNew[i] = pointi++;
        }
    }

    // report unused points
    if (nPoints > pointi)
    {
        Info<< "Unused    points  = " << (nPoints - pointi) << endl;
        nPoints = pointi;

        // adjust points and truncate
        inplaceReorder(oldToNew, points_);
        points_.setSize(nPoints);

        // adjust cellFaces - with mesh shapes this might be faster
        forAll(cellFaces_, celli)
        {
            faceList& faces = cellFaces_[celli];
            forAll(faces, i)
            {
                inplaceRenumber(oldToNew, faces[i]);
            }
        }

        // adjust baffles
        forAll(baffleFaces_, facei)
        {
            inplaceRenumber(oldToNew, baffleFaces_[facei]);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::meshReaders::STARCD::readGeometry(const scalar scaleFactor)
{
    readPoints(geometryFile_ + ".vrt", scaleFactor);
    readCells(geometryFile_ + ".cel");
    cullPoints();
    readBoundary(geometryFile_ + ".bnd");

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshReaders::STARCD::STARCD
(
    const fileName& prefix,
    const objectRegistry& registry,
    const scalar scaleFactor
)
:
    meshReader(prefix, scaleFactor),
    cellShapes_(0),
    mapToFoamPointId_(0),
    mapToFoamCellId_(0)
{
    readAux(registry);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshReaders::STARCD::~STARCD()
{}


// ************************************************************************* //
