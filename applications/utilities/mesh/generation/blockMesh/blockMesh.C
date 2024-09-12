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

Application
    blockMesh

Description
    A multi-block mesh generator.

    Uses the block mesh description found in
      - \c system/blockMeshDict
      - \c system/\<region\>/blockMeshDict
      - \c constant/polyMesh/blockMeshDict
      - \c constant/\<region\>/polyMesh/blockMeshDict

Usage
    \b blockMesh [OPTION]

    Options:
      - \par -blockTopology
        Write the topology as a set of edges in OBJ format.

      - \par -region \<name\>
        Specify an alternative mesh region.

      - \par -dict \<filename\>
        Specify alternative dictionary for the block mesh description.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "systemDict.H"
#include "blockMesh.H"
#include "mergePatchPairs.H"
#include "polyTopoChange.H"
#include "emptyPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    Foam::argList::removeOption("noFunctionObjects");
    #include "addDictOption.H"
    argList::addBoolOption
    (
        "blockTopology",
        "write block edges and centres as .obj files"
    );
    argList::addBoolOption
    (
        "noClean",
        "keep the existing files in the polyMesh"
    );

    argList::addNote
    (
        "Block description\n"
        "\n"
        "  For a given block, the correspondence between the ordering of\n"
        "  vertex labels and face labels is shown below.\n"
        "  For vertex numbering in the sequence 0 to 7 (block, centre):\n"
        "    faces 0 (f0) and 1 are left and right, respectively;\n"
        "    faces 2 and 3 are front and back; \n"
        "    and faces 4 and 5 are bottom and top::\n"
        "\n"
        "                 7 ---- 6\n"
        "            f5   |\\     |\\   f3\n"
        "            |    | 4 ---- 5   \\\n"
        "            |    3 |--- 2 |    \\\n"
        "            |     \\|     \\|    f2\n"
        "            f4     0 ---- 1\n"
        "\n"
        "       Z         f0 ----- f1\n"
        "       |  Y\n"
        "       | /\n"
        "       O --- X\n"
    );

    #include "addMeshOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "setMeshPath.H"
    #include "createTime.H"

    const word dictName("blockMeshDict");


    // Check if the mesh is specified otherwise mesh the default mesh
    if (!meshPath.empty())
    {
        Info<< nl << "Generating mesh " << meshPath << endl;
    }

    word regionName;
    word regionPath;

    // Check if the region is specified otherwise mesh the default region
    if (args.optionReadIfPresent("region", regionName, polyMesh::defaultRegion))
    {
        Info<< nl << "Generating mesh for region " << regionName << endl;
        regionPath = regionName;
    }

    if (!args.optionFound("noClean"))
    {
        const fileName polyMeshPath
        (
            runTime.path()/runTime.constant()
           /meshPath/regionPath/polyMesh::meshSubDir
        );

        if (exists(polyMeshPath))
        {
            if (exists(polyMeshPath/dictName))
            {
                Info<< "Not deleting polyMesh directory " << nl
                    << "    " << polyMeshPath << nl
                    << "    because it contains " << dictName << endl;
            }
            else
            {
                Info<< "Deleting polyMesh directory" << nl
                    << "    " << polyMeshPath << endl;
                rmDir(polyMeshPath);
            }
        }
    }

    typeIOobject<IOdictionary> meshDictIO
    (
        systemDictIO(dictName, args, runTime, regionName, meshPath)
    );

    if (!meshDictIO.headerOk())
    {
        FatalErrorInFunction
            << "Cannot find file " << meshDictIO.relativeObjectPath()
            << nl
            << exit(FatalError);
    }

    Info<< "Creating block mesh from\n    "
        << meshDictIO.relativeObjectPath() << endl;

    IOdictionary meshDict(meshDictIO);
    blockMesh blocks(meshDict, runTime.constant()/meshPath, regionName);


    if (args.optionFound("blockTopology"))
    {
        // Write mesh as edges.
        {
            fileName objMeshFile("blockTopology.obj");

            OFstream str(runTime.path()/objMeshFile);

            Info<< nl << "Dumping block structure as Lightwave obj format"
                << " to " << objMeshFile << endl;

            blocks.writeTopology(str);
        }

        // Write centres of blocks
        {
            fileName objCcFile("blockCentres.obj");

            OFstream str(runTime.path()/objCcFile);

            Info<< nl << "Dumping block centres as Lightwave obj format"
                << " to " << objCcFile << endl;

            const polyMesh& topo = blocks.topology();

            const pointField& cellCentres = topo.cellCentres();

            forAll(cellCentres, celli)
            {
                // point cc = b.blockShape().centre(b.points());
                const point& cc = cellCentres[celli];

                str << "v " << cc.x() << ' ' << cc.y() << ' ' << cc.z() << nl;
            }
        }

        Info<< nl << "end" << endl;

        return 0;
    }


    Info<< nl << "Creating polyMesh from blockMesh" << endl;

    word defaultFacesName = "defaultFaces";
    word defaultFacesType = emptyPolyPatch::typeName;
    polyMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            meshPath,
            runTime
        ),
        clone(blocks.points()),           // could we reuse space?
        blocks.cells(),
        blocks.patches(),
        blocks.patchNames(),
        blocks.patchDicts(),
        defaultFacesName,
        defaultFacesType
    );


    // Read in a list of dictionaries for the merge patch pairs
    if (meshDict.found("mergePatchPairs"))
    {
        List<Pair<word>> patchPairNames
        (
            meshDict.lookup("mergePatchPairs")
        );

        if (patchPairNames.size())
        {
            const word oldInstance = mesh.pointsInstance();
            mergePatchPairs(mesh, patchPairNames, 1e-4);
            mesh.setInstance(oldInstance);
        }
    }
    else
    {
        Info<< nl << "No patch pairs to merge" << endl;
    }

    label nZones = blocks.numZonedBlocks();

    if (nZones > 0)
    {
        Info<< nl << "Adding cell zones" << endl;

        // Map from zoneName to cellZone index
        HashTable<label> zoneMap(nZones);

        // Cells per zone.
        List<DynamicList<label>> zoneCells(nZones);

        // Running cell counter
        label celli = 0;

        // Largest zone so far
        label freeZoneI = 0;

        forAll(blocks, blockI)
        {
            const block& b = blocks[blockI];
            const List<FixedList<label, 8>> blockCells = b.cells();
            const word& zoneName = b.zoneName();

            if (zoneName.size())
            {
                HashTable<label>::const_iterator iter = zoneMap.find(zoneName);

                label zoneI;

                if (iter == zoneMap.end())
                {
                    zoneI = freeZoneI++;

                    Info<< "    " << zoneI << '\t' << zoneName << endl;

                    zoneMap.insert(zoneName, zoneI);
                }
                else
                {
                    zoneI = iter();
                }

                forAll(blockCells, i)
                {
                    zoneCells[zoneI].append(celli++);
                }
            }
            else
            {
                celli += blockCells.size();
            }
        }


        List<cellZone*> cz(zoneMap.size());

        forAllConstIter(HashTable<label>, zoneMap, iter)
        {
            label zoneI = iter();

            cz[zoneI] = new cellZone
            (
                iter.key(),
                zoneCells[zoneI].shrink(),
                mesh.cellZones()
            );
        }

        mesh.pointZones().setSize(0);
        mesh.faceZones().setSize(0);
        mesh.cellZones().setSize(0);
        mesh.addZones(List<pointZone*>(0), List<faceZone*>(0), cz);
    }


    // Detect any cyclic patches and force re-ordering of the faces
    {
        const polyPatchList& patches = mesh.boundaryMesh();
        bool hasCyclic = false;
        forAll(patches, patchi)
        {
            if (isA<cyclicPolyPatch>(patches[patchi]))
            {
                hasCyclic = true;
                break;
            }
        }

        if (hasCyclic)
        {
            Info<< nl << "Detected cyclic patches; ordering boundary faces"
                << endl;
            const word oldInstance = mesh.instance();
            polyTopoChange meshMod(mesh);
            meshMod.changeMesh(mesh);
            mesh.setInstance(oldInstance);
        }
    }


    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< nl << "Writing polyMesh" << endl;
    mesh.removeFiles();
    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed writing polyMesh."
            << exit(FatalError);
    }


    // Write summary
    {
        const polyPatchList& patches = mesh.boundaryMesh();

        Info<< "----------------" << nl
            << "Mesh Information" << nl
            << "----------------" << nl
            << "  " << "boundingBox: " << boundBox(mesh.points()) << nl
            << "  " << "nPoints: " << mesh.nPoints() << nl
            << "  " << "nCells: " << mesh.nCells() << nl
            << "  " << "nFaces: " << mesh.nFaces() << nl
            << "  " << "nInternalFaces: " << mesh.nInternalFaces() << nl;

        Info<< "----------------" << nl
            << "Patches" << nl
            << "----------------" << nl;

        forAll(patches, patchi)
        {
            const polyPatch& p = patches[patchi];

            Info<< "  " << "patch " << patchi
                << " (start: " << p.start()
                << " size: " << p.size()
                << ") name: " << p.name()
                << nl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
