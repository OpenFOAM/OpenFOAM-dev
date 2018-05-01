/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

Description
    Reads CCM files as written by Prostar/ccm using ccm 2.6 (not 2.4)

    - does polyhedral mesh
    - does not handle 'interfaces' (couples)
    - does not handle cyclics. Use createPatch to recreate these
    - does not do data
    - does patch names only if they are in the problem description

\*---------------------------------------------------------------------------*/

#include "ListOps.H"
#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "emptyPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wallPolyPatch.H"
#include "SortableList.H"
#include "cellSet.H"

#include <ccmio.h>
#include <vector>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static char const kDefaultState[] = "default";
static int const kVertOffset = 2;


// Determine upper-triangular order for internal faces:
labelList getInternalFaceOrder
(
    const cellList& cells,
    const labelList& owner,
    const labelList& neighbour
)
{
    labelList oldToNew(owner.size(), -1);

    // First unassigned face
    label newFacei = 0;

    forAll(cells, celli)
    {
        const labelList& cFaces = cells[celli];

        SortableList<label> nbr(cFaces.size());

        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            label nbrCelli = neighbour[facei];

            if (nbrCelli != -1)
            {
                // Internal face. Get cell on other side.
                if (nbrCelli == celli)
                {
                    nbrCelli = owner[facei];
                }

                if (celli < nbrCelli)
                {
                    // Celli is master
                    nbr[i] = nbrCelli;
                }
                else
                {
                    // nbrCell is master. Let it handle this face.
                    nbr[i] = -1;
                }
            }
            else
            {
                // External face. Do later.
                nbr[i] = -1;
            }
        }

        nbr.sort();

        forAll(nbr, i)
        {
            if (nbr[i] != -1)
            {
                oldToNew[cFaces[nbr.indices()[i]]] = newFacei++;
            }
        }
    }

    // Keep boundary faces in same order.
    for (label facei = newFacei; facei < owner.size(); facei++)
    {
        oldToNew[facei] = facei;
    }

    return oldToNew;
}


void storeCellInZone
(
    const label celli,
    const label cellType,
    Map<label>& typeToZone,
    List<DynamicList<label>>& zoneCells
)
{
    if (cellType >= 0)
    {
        Map<label>::iterator zoneFnd = typeToZone.find(cellType);

        if (zoneFnd == typeToZone.end())
        {
            // New type. Allocate zone for it.
            zoneCells.setSize(zoneCells.size() + 1);

            label zoneI = zoneCells.size()-1;

            Info<< "Mapping type " << cellType << " to Foam cellZone "
                << zoneI << endl;
            typeToZone.insert(cellType, zoneI);

            zoneCells[zoneI].append(celli);
        }
        else
        {
            // Existing zone for type
            zoneCells[zoneFnd()].append(celli);
        }
    }
}


void CheckError(CCMIOError const &err, const Foam::string& str)
{
    if (err != kCCMIONoErr)
    {
        FatalErrorInFunction
            << str << exit(FatalError);
    }
}


void ReadVertices
(
    CCMIOError &err,
    CCMIOID &vertices,
    labelList& foamPointMap,
    pointField& foamPoints
)
{

    // Read the vertices.  This involves reading both the vertex data and
    // the map, which maps the index into the data array with the ID number.
    // As we process the vertices we need to be sure to scale them by the
    // appropriate scaling factor.  The offset is just to show you can read
    // any chunk.  Normally this would be in a for loop.

    CCMIOSize nVertices;
    CCMIOEntitySize(&err, vertices, &nVertices, nullptr);

    List<int> mapData(nVertices, 0);
    List<float> verts(3*nVertices, 0);

    int offset = 0;
    int offsetPlusSize = offset+nVertices;

    int dims = 1;
    float scale;
    CCMIOID mapID;
    CCMIOReadVerticesf
    (
        &err, vertices, &dims, &scale, &mapID, verts.begin(),
        offset, offsetPlusSize
    );
    CCMIOReadMap(&err, mapID, mapData.begin(), offset, offsetPlusSize);

    // CCMIOSize size;
    // CCMIOEntityDescription(&err, vertices, &size, nullptr);
    // char *desc = new char[size + 1];
    // CCMIOEntityDescription(&err, vertices, nullptr, desc);
    // Pout<< "label: '" << desc << "'" << endl;
    // delete [] desc;

    // Convert to foamPoints
    foamPoints.setSize(nVertices);
    foamPoints = Zero;
    foamPointMap.setSize(nVertices);

    forAll(foamPointMap, i)
    {
        foamPointMap[i] = mapData[i];
        for (direction cmpt = 0; cmpt < dims; cmpt++)
        {
            foamPoints[i][cmpt] = verts[dims*i + cmpt]*scale;
        }
    }
}


void ReadProblem
(
    CCMIOError &err,
    CCMIOID& problem,

    const Map<label>& prostarToFoamPatch,
    Map<word>& foamCellTypeNames,
    wordList& foamPatchTypes,
    wordList& foamPatchNames
)
{
    // ... walk through each cell type and print it...
    CCMIOID next;
    int i = 0;
    while
    (
        CCMIONextEntity(nullptr, problem, kCCMIOCellType, &i, &next)
     == kCCMIONoErr
    )
    {
        char *name;
        int size;
        int cellType;

        // ... if it has a material type.  (Note that we do not pass in
        // an array to get the name because we do not know how long the
        // string is yet.  Many parameters to CCMIO functions that
        // return
        // data can be nullptr if that data is not needed.)
        if
        (
            CCMIOReadOptstr(nullptr, next, "MaterialType", &size, nullptr)
         == kCCMIONoErr
        )
        {
            name = new char[size + 1];
            CCMIOReadOptstr(&err, next, "MaterialType", &size, name);
            CCMIOGetEntityIndex(&err, next, &cellType);

            foamCellTypeNames.insert(cellType, name);
            Pout<< "Celltype:" << cellType << " name:" << name << endl;

            delete [] name;
        }
    }

    // ... walk through each region description and print it...

    CCMIOID boundary;
    label regionI = 0;
    int k = 0;
    while
    (
        CCMIONextEntity(nullptr, problem, kCCMIOBoundaryRegion, &k, &boundary)
     == kCCMIONoErr
    )
    {
        // Index of foam patch
        label foamPatchi = -1;

        // Read prostar id

        int prostarI = -1;
        if
        (
            CCMIOReadOpti(nullptr, boundary, "ProstarRegionNumber", &prostarI)
         == kCCMIONoErr
        )
        {
            Pout<< "For region:" << regionI
                << " found ProstarRegionNumber:" << prostarI << endl;
        }
        else
        {
            prostarI = regionI;

            Pout<< "For region:" << regionI
                << "did not find ProstarRegionNumber entry. Assuming "
                << prostarI << endl;
        }


        if (prostarToFoamPatch.found(prostarI))
        {
            foamPatchi = prostarToFoamPatch[prostarI];

            // Read boundary type

            int size;
            if
            (
                CCMIOReadOptstr
                (
                    nullptr,
                    boundary,
                    "BoundaryType",
                    &size,
                    nullptr
                ) == kCCMIONoErr
            )
            {
                char* s = new char[size + 1];
                CCMIOReadOptstr(nullptr, boundary, "BoundaryType", &size, s);
                s[size] = '\0';
                foamPatchTypes[foamPatchi] = string::validate<word>(string(s));
                delete [] s;
            }


            // foamPatchMap.append(prostarI);


            // Read boundary name:
            // - from BoundaryName field (Prostar)
            // - from 'Label' field (ccm+?)
            // - make up one from prostar id.

            if
            (
                CCMIOReadOptstr
                (
                    nullptr,
                    boundary,
                    "BoundaryName",
                    &size,
                    nullptr
                ) == kCCMIONoErr
            )
            {
                char* name = new char[size + 1];
                CCMIOReadOptstr(nullptr, boundary, "BoundaryName", &size, name);
                name[size] = '\0';
                foamPatchNames[foamPatchi] =
                    string::validate<word>(string(name));
                delete [] name;
            }
            else if
            (
                CCMIOReadOptstr(nullptr, boundary, "Label", &size, nullptr)
             == kCCMIONoErr
            )
            {
                char* name = new char[size + 1];
                CCMIOReadOptstr(nullptr, boundary, "Label", &size, name);
                name[size] = '\0';
                foamPatchNames[foamPatchi] =
                    string::validate<word>(string(name));
                delete [] name;
            }
            else
            {
                foamPatchNames[foamPatchi] =
                    foamPatchTypes[foamPatchi]
                  + Foam::name(foamPatchi);
                Pout<< "Made up name:" << foamPatchNames[foamPatchi]
                    << endl;
            }

            Pout<< "Read patch:" << foamPatchi
                << " name:" << foamPatchNames[foamPatchi]
                << " foamPatchTypes:" << foamPatchTypes[foamPatchi]
                << endl;
        }

        regionI++;
    }
}


void ReadCells
(
    CCMIOError &err,
    CCMIOID& topology,
    labelList& foamCellMap,
    labelList& foamCellType,
    Map<label>& prostarToFoamPatch,
    DynamicList<label>& foamPatchSizes,
    DynamicList<label>& foamPatchStarts,
    labelList& foamFaceMap,
    labelList& foamOwner,
    labelList& foamNeighbour,
    faceList& foamFaces
)
{
    // Read the cells.
    // ~~~~~~~~~~~~~~~

    //  Store cell IDs so that we know what cells are in
    // this post data.
    CCMIOID id;
    CCMIOGetEntity(&err, topology, kCCMIOCells, 0, &id);
    CCMIOSize nCells;
    CCMIOEntitySize(&err, id, &nCells, nullptr);

    std::vector<int> mapData(nCells);
    std::vector<int> cellType(nCells);

    CCMIOID mapID;
    CCMIOReadCells(&err, id, &mapID, &cellType[0], 0, nCells);
    CCMIOReadMap(&err, mapID, &mapData[0], 0, nCells);
    CheckError(err, "Error reading cells");

    foamCellMap.setSize(nCells);
    foamCellType.setSize(nCells);
    forAll(foamCellMap, i)
    {
        foamCellMap[i] = mapData[i];
        foamCellType[i] = cellType[i];
    }


    // Read the internal faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    CCMIOGetEntity(&err, topology, kCCMIOInternalFaces, 0, &id);
    CCMIOSize nInternalFaces;
    CCMIOEntitySize(&err, id, &nInternalFaces, nullptr);
    Pout<< "nInternalFaces:" << nInternalFaces << endl;

    // Determine patch sizes before reading internal faces
    label foamNFaces = nInternalFaces;
    int index = 0;
    while
    (
        CCMIONextEntity(nullptr, topology, kCCMIOBoundaryFaces, &index, &id)
     == kCCMIONoErr
    )
    {
        CCMIOSize size;
        CCMIOEntitySize(&err, id, &size, nullptr);

        Pout<< "Read kCCMIOBoundaryFaces entry with " << size
            << " faces." << endl;

        foamPatchStarts.append(foamNFaces);
        foamPatchSizes.append(size);
        foamNFaces += size;
    }
    foamPatchStarts.shrink();
    foamPatchSizes.shrink();

    Pout<< "patchSizes:" << foamPatchSizes << endl;
    Pout<< "patchStarts:" << foamPatchStarts << endl;
    Pout<< "nFaces:" << foamNFaces << endl;

    mapData.resize(nInternalFaces);
    CCMIOGetEntity(&err, topology, kCCMIOInternalFaces, 0, &id);
    CCMIOSize size;
    CCMIOReadFaces(&err, id, kCCMIOInternalFaces, nullptr, &size, nullptr,
                   kCCMIOStart, kCCMIOEnd);
    std::vector<int> faces(size);
    CCMIOReadFaces(&err, id, kCCMIOInternalFaces, &mapID, nullptr, &faces[0],
                   kCCMIOStart, kCCMIOEnd);
    std::vector<int> faceCells(2*nInternalFaces);
    CCMIOReadFaceCells(&err, id, kCCMIOInternalFaces, &faceCells[0],
                       kCCMIOStart, kCCMIOEnd);
    CCMIOReadMap(&err, mapID, &mapData[0], kCCMIOStart, kCCMIOEnd);
    CheckError(err, "Error reading internal faces");

    // Copy into Foam lists
    foamFaceMap.setSize(foamNFaces);
    foamFaces.setSize(foamNFaces);
    foamOwner.setSize(foamNFaces);
    foamNeighbour.setSize(foamNFaces);

    unsigned int pos = 0;

    for (unsigned int facei = 0; facei < nInternalFaces; facei++)
    {
        foamFaceMap[facei] = mapData[facei];
        foamOwner[facei] = faceCells[2*facei];
        foamNeighbour[facei] = faceCells[2*facei+1];
        face& f = foamFaces[facei];

        f.setSize(faces[pos++]);
        forAll(f, fp)
        {
            f[fp] = faces[pos++];
        }
    }

    // Read the boundary faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    index = 0;
    label regionI = 0;
    while
    (
        CCMIONextEntity(nullptr, topology, kCCMIOBoundaryFaces, &index, &id)
     == kCCMIONoErr
    )
    {
        CCMIOSize nFaces;
        CCMIOEntitySize(&err, id, &nFaces, nullptr);

        mapData.resize(nFaces);
        faceCells.resize(nFaces);
        CCMIOReadFaces(&err, id, kCCMIOBoundaryFaces, nullptr, &size, nullptr,
                       kCCMIOStart, kCCMIOEnd);
        faces.resize(size);
        CCMIOReadFaces
        (
            &err,
            id,
            kCCMIOBoundaryFaces,
            &mapID,
            nullptr,
            &faces[0],
            kCCMIOStart,
            kCCMIOEnd
        );
        CCMIOReadFaceCells
        (
            &err,
            id,
            kCCMIOBoundaryFaces,
            &faceCells[0],
            kCCMIOStart,
            kCCMIOEnd
        );
        CCMIOReadMap(&err, mapID, &mapData[0], kCCMIOStart, kCCMIOEnd);
        CheckError(err, "Error reading boundary faces");

        // Read prostar id
        int prostarI;
        if
        (
            CCMIOReadOpti(nullptr, id, "ProstarRegionNumber", &prostarI)
         == kCCMIONoErr
        )
        {
            Pout<< "For region:" << regionI
                << " found ProstarRegionNumber:" << prostarI << endl;
        }
        else
        {
            prostarI = regionI;

            Pout<< "For region:" << regionI
                << " did not find ProstarRegionNumber entry. Assuming "
                << prostarI << endl;
        }
        prostarToFoamPatch.insert(prostarI, regionI);


        Pout<< "region:" << regionI
            << " ProstarRegionNumber:" << prostarI
            << " foamPatchStart:"
            << foamPatchStarts[regionI]
            << " size:"
            << foamPatchSizes[regionI]
            << endl;

        // Copy into Foam list.
        unsigned int pos = 0;

        for (unsigned int i = 0; i < nFaces; i++)
        {
            label foamFacei = foamPatchStarts[regionI] + i;

            foamFaceMap[foamFacei] = mapData[i];
            foamOwner[foamFacei] = faceCells[i];
            foamNeighbour[foamFacei] = -1;
            face& f = foamFaces[foamFacei];

            f.setSize(faces[pos++]);
            forAll(f, fp)
            {
                f[fp] = faces[pos++];
            }
        }
        Pout<< endl;

        regionI++;
    }
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "read CCM files as written by proSTAR/ccm\n"
        " - does not handle 'interfaces' (couples), cyclics or data\n"
        " - does not handle mesh regions (porosity, solids, ...)\n"
    );
    argList::noParallel();
    argList::validArgs.append("ccmFile");

    #include "setRootCase.H"
    #include "createTime.H"

    // Foam mesh data
    // ~~~~~~~~~~~~~~

    // Coordinates
    pointField foamPoints;
    labelList foamPointMap;

    // Cell type
    labelList foamCellType;
    labelList foamCellMap;

    // Patching info
    Map<label> prostarToFoamPatch;
    DynamicList<label> foamPatchSizes;
    DynamicList<label> foamPatchStarts;
    // Face connectivity
    labelList foamFaceMap;
    labelList foamOwner;
    labelList foamNeighbour;
    faceList foamFaces;

    // Celltypes (not the names of the zones; just the material type)
    // and patch type names
    Map<word> foamCellTypeNames;
    wordList foamPatchTypes;
    wordList foamPatchNames;

    {
        const fileName ccmFile = args[1];
        const word ccmExt = ccmFile.ext();

        if (!isFile(ccmFile))
        {
            FatalErrorInFunction
                << "Cannot read file " << ccmFile
                << exit(FatalError);
        }

        if (ccmExt != "ccm" && ccmExt != "ccmg")
        {
            FatalErrorInFunction
                << "Illegal extension " << ccmExt << " for file " << ccmFile
                << nl << "Allowed extensions are '.ccm', '.ccmg'"
                << exit(FatalError);
        }

        // Open the file.  Because we did not initialize 'err' we need to pass
        // in nullptr (which always means kCCMIONoErr) and then assign the
        // return value to 'err'.).
        CCMIOID root;
        CCMIOError err = CCMIOOpenFile
        (
            nullptr,
            ccmFile.c_str(),
            kCCMIORead,
            &root
        );

        // We are going to assume that we have a state with a known name.
        // We could instead use CCMIONextEntity() to walk through all the
        // states in the file and present the list to the user for selection.
        CCMIOID state;
        int stateI = 0;
        CCMIONextEntity(&err, root, kCCMIOState, &stateI, &state);
        CheckError(err, "Error opening state");

        unsigned int size;
        CCMIOEntityDescription(&err, state, &size, nullptr);
        char *desc = new char[size + 1];
        CCMIOEntityDescription(&err, state, nullptr, desc);
        Pout<< "Reading state '" << kDefaultState << "' (" << desc << ")"
            << endl;
        delete [] desc;

        // Find the first processor (i has previously been initialized to 0) and
        // read the mesh and solution information.
        int i = 0;
        CCMIOID processor;
        CCMIONextEntity(&err, state, kCCMIOProcessor, &i, &processor);
        CCMIOID solution, vertices, topology;
        CCMIOReadProcessor
        (
            &err,
            processor,
            &vertices,
            &topology,
            nullptr,
            &solution
        );

        if (err != kCCMIONoErr)
        {
            // Maybe no solution;  try again
            err = kCCMIONoErr;
            CCMIOReadProcessor
            (
                &err,
                processor,
                &vertices,
                &topology,
                nullptr,
                nullptr
            );
            if (err != kCCMIONoErr)
            {
                FatalErrorInFunction
                    << "Could not read the file."
                    << exit(FatalError);
            }
        }

        ReadVertices(err, vertices, foamPointMap, foamPoints);

        Pout<< "nPoints:" << foamPoints.size() << endl
            << "bounding box:" << boundBox(foamPoints) << endl
            << endl;

        ReadCells
        (
            err,
            topology,
            foamCellMap,
            foamCellType,
            prostarToFoamPatch,
            foamPatchSizes,
            foamPatchStarts,
            foamFaceMap,
            foamOwner,
            foamNeighbour,
            foamFaces
        );

        Pout<< "nCells:" << foamCellMap.size() << endl
            << "nFaces:" << foamOwner.size() << endl
            << "nPatches:" << foamPatchStarts.size() << endl
            << "nInternalFaces:" << foamPatchStarts[0] << endl
            << endl;

        // Create some default patch names/types. These will be overwritten
        // by any problem description (if it is there)
        foamPatchTypes.setSize(foamPatchStarts.size());
        foamPatchNames.setSize(foamPatchStarts.size());

        forAll(foamPatchNames, i)
        {
            foamPatchNames[i] = word("patch") + Foam::name(i);
            foamPatchTypes[i] = "patch";
        }

        // Get problem description

        CCMIOID problem;
        int problemI = 0;
        CCMIONextEntity
        (
            &err,
            root,
            kCCMIOProblemDescription,
            &problemI,
            &problem
        );
        CheckError(err, "Error stepping to first problem description");

        if (CCMIOIsValidEntity(problem))   // if we have a problem description
        {
            ReadProblem
            (
                err,
                problem,
                prostarToFoamPatch,

                foamCellTypeNames,
                foamPatchTypes,
                foamPatchNames
            );
        }


        CCMIOCloseFile(&err, vertices);
        CCMIOCloseFile(&err, topology);
        CCMIOCloseFile(&err, solution);
        CCMIOCloseFile(&err, root);
    }


    Pout<< "foamPatchNames:" << foamPatchNames << endl;


    Pout<< "foamOwner : min:" << min(foamOwner)
        << " max:" << max(foamOwner)
        << nl
        << "foamNeighbour : min:" << min(foamNeighbour)
        << " max:" << max(foamNeighbour)
        << nl
        << "foamCellType : min:" << min(foamCellType)
        << " max:" << max(foamCellType)
        << nl << endl;



    // We now have extracted all info from CCMIO:
    // - coordinates (points)
    // - face to point addressing (faces)
    // - face to cell addressing (owner, neighbour)
    // - cell based data (cellId)


    // Renumber vertex labels to Foam point labels
    {
        label maxCCMPointi = max(foamPointMap);
        labelList toFoamPoints(invert(maxCCMPointi+1, foamPointMap));

        forAll(foamFaces, facei)
        {
            inplaceRenumber(toFoamPoints, foamFaces[facei]);
        }
    }

    // Renumber cell labels
    {
        label maxCCMCelli = max(foamCellMap);
        labelList toFoamCells(invert(maxCCMCelli+1, foamCellMap));

        inplaceRenumber(toFoamCells, foamOwner);
        inplaceRenumber(toFoamCells, foamNeighbour);
    }


    //
    // Basic mesh info complete. Now convert to Foam convention:
    // - owner < neighbour
    // - face vertices such that normal points away from owner
    // - order faces: upper-triangular for internal faces; boundary faces after
    //   internal faces
    //

    // Set owner/neighbour so owner < neighbour
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(foamNeighbour, facei)
    {
        label nbr = foamNeighbour[facei];
        label own = foamOwner[facei];

        if (nbr >= foamCellType.size() || own >= foamCellType.size())
        {
            FatalErrorInFunction
                << "face:" << facei
                << " nbr:" << nbr
                << " own:" << own
                << " nCells:" << foamCellType.size()
                << exit(FatalError);
        }

        if (nbr >= 0)
        {
            if (nbr < own)
            {
                foamOwner[facei] = foamNeighbour[facei];
                foamNeighbour[facei] = own;
                foamFaces[facei].flip();
            }
        }


        // And check the face
        const face& f = foamFaces[facei];

        forAll(f, fp)
        {
            if (f[fp] < 0 || f[fp] >= foamPoints.size())
            {
                FatalErrorInFunction
                    << " f:" << f << abort(FatalError);
            }
        }
    }


    // Do upper-triangular ordering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Create cells (inverse of face-to-cell addressing)
        cellList foamCells;
        primitiveMesh::calcCells
        (
            foamCells,
            foamOwner,
            foamNeighbour,
            foamCellType.size()
        );

        // Determine face order for upper-triangular ordering
        labelList oldToNew
        (
            getInternalFaceOrder
            (
                foamCells,
                foamOwner,
                foamNeighbour
            )
        );

        // Reorder faces accordingly
        forAll(foamCells, celli)
        {
            inplaceRenumber(oldToNew, foamCells[celli]);
        }

        // Reorder faces.
        inplaceReorder(oldToNew, foamFaces);
        inplaceReorder(oldToNew, foamOwner);
        inplaceReorder(oldToNew, foamNeighbour);
    }


    // Construct fvMesh
    // ~~~~~~~~~~~~~~~~

    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        xferMove<pointField>(foamPoints),
        xferMove<faceList>(foamFaces),
        xferCopy<labelList>(foamOwner),
        xferMove<labelList>(foamNeighbour)
    );

    // Create patches. Use patch types to determine what Foam types to generate.
    List<polyPatch*> newPatches(foamPatchNames.size());

    label meshFacei = foamPatchStarts[0];

    forAll(newPatches, patchi)
    {
        const word& patchName = foamPatchNames[patchi];
        const word& patchType = foamPatchTypes[patchi];

        Pout<< "Patch:" << patchName << " start at:" << meshFacei
            << " size:" << foamPatchSizes[patchi]
            << " end at:" << meshFacei+foamPatchSizes[patchi]
            << endl;

        if (patchType == "wall")
        {
            newPatches[patchi] =
                new wallPolyPatch
                (
                    patchName,
                    foamPatchSizes[patchi],
                    meshFacei,
                    patchi,
                    mesh.boundaryMesh(),
                    patchType
                );
        }
        else if (patchType == "symmetryplane")
        {
            newPatches[patchi] =
                new symmetryPolyPatch
                (
                    patchName,
                    foamPatchSizes[patchi],
                    meshFacei,
                    patchi,
                    mesh.boundaryMesh(),
                    patchType
                );
        }
        else if (patchType == "empty")
        {
            // Note: not ccm name, introduced by us above.
            newPatches[patchi] =
                new emptyPolyPatch
                (
                    patchName,
                    foamPatchSizes[patchi],
                    meshFacei,
                    patchi,
                    mesh.boundaryMesh(),
                    patchType
                );
        }
        else
        {
            // All other ccm types become straight polyPatch:
            // 'inlet', 'outlet', ...
            newPatches[patchi] =
                new polyPatch
                (
                    patchName,
                    foamPatchSizes[patchi],
                    meshFacei,
                    patchi,
                    mesh.boundaryMesh(),
                    word::null
                );
        }

        meshFacei += foamPatchSizes[patchi];
    }

    if (meshFacei != foamOwner.size())
    {
        FatalErrorInFunction
            << "meshFacei:" << meshFacei
            << " nFaces:" << foamOwner.size()
            << abort(FatalError);
    }
    mesh.addFvPatches(newPatches);



    // Construct sets and zones from types
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label maxType = max(foamCellType);
    label minType = min(foamCellType);

    if (maxType > minType)
    {
        // From foamCellType physical region to Foam cellZone
        Map<label> typeToZone;
        // Storage for cell zones.
        List<DynamicList<label>> zoneCells(0);

        forAll(foamCellType, celli)
        {
            storeCellInZone
            (
                celli,
                foamCellType[celli],
                typeToZone,
                zoneCells
            );
        }

        mesh.cellZones().clear();
        mesh.cellZones().setSize(typeToZone.size());

        label nValidCellZones = 0;

        forAllConstIter(Map<label>, typeToZone, iter)
        {
            label type = iter.key();
            label zoneI = iter();

            zoneCells[zoneI].shrink();

            word zoneName = "cellZone_" + name(type);

            Info<< "Writing zone " << type
                << " size " << zoneCells[zoneI].size()
                << " to cellZone "
                << zoneName << " and cellSet " << zoneName
                << endl;

            cellSet cset(mesh, zoneName, zoneCells[zoneI]);
            cset.write();

            mesh.cellZones().set
            (
                nValidCellZones,
                new cellZone
                (
                    zoneName,
                    zoneCells[zoneI],
                    nValidCellZones,
                    mesh.cellZones()
                )
            );
            nValidCellZones++;
        }
        mesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;
    }


    Info<< "Writing mesh to " << mesh.objectRegistry::objectPath()
        << "..." << nl << endl;


    // Construct field with calculated bc to hold Star cell Id.
    volScalarField cellIdField
    (
        IOobject
        (
            "cellId",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("cellId", dimless, 0.0)
    );

    forAll(foamCellMap, celli)
    {
        cellIdField[celli] = foamCellMap[celli];
    }

    // Construct field with calculated bc to hold cell type.
    volScalarField cellTypeField
    (
        IOobject
        (
            "cellType",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("cellType", dimless, 0.0)
    );

    forAll(foamCellType, celli)
    {
        cellTypeField[celli] = foamCellType[celli];
    }

    Info<< "Writing cellIds as volScalarField to " << cellIdField.objectPath()
        << "..." << nl << endl;
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
