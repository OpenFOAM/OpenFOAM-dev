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

Application
    blockMesh

Description
    A multi-block mesh generator.

    Uses the block mesh description found in
    \a constant/polyMesh/blockMeshDict
    (or \a constant/\<region\>/polyMesh/blockMeshDict).

Usage

    - blockMesh [OPTION]

    \param -blockTopology \n
    Write the topology as a set of edges in OBJ format.

    \param -region \<name\> \n
    Specify an alternative mesh region.

    \param -dict \<filename\> \n
    Specify alternative dictionary for the block mesh description.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "IOdictionary.H"
#include "IOPtrList.H"

#include "blockMesh.H"
#include "attachPolyTopoChanger.H"
#include "emptyPolyPatch.H"
#include "cellSet.H"

#include "argList.H"
#include "OSspecific.H"
#include "OFstream.H"

#include "Pair.H"
#include "slidingInterface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption
    (
        "blockTopology",
        "write block edges and centres as .obj files"
    );
    argList::addOption
    (
        "dict",
        "file",
        "specify alternative dictionary for the blockMesh description"
    );

#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"

    const word dictName("blockMeshDict");

    word regionName;
    fileName polyMeshDir;

    if (args.optionReadIfPresent("region", regionName, polyMesh::defaultRegion))
    {
        // constant/<region>/polyMesh/blockMeshDict
        polyMeshDir = regionName/polyMesh::meshSubDir;

        Info<< nl << "Generating mesh for region " << regionName << endl;
    }
    else
    {
        // constant/polyMesh/blockMeshDict
        polyMeshDir = polyMesh::meshSubDir;
    }

    IOobject meshDictIO
    (
        dictName,
        runTime.constant(),
        polyMeshDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (args.optionFound("dict"))
    {
        const fileName dictPath = args["dict"];

        meshDictIO = IOobject
        (
            (
                isDir(dictPath)
              ? dictPath/dictName
              : dictPath
            ),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );
    }

    if (!meshDictIO.headerOk())
    {
        FatalErrorIn(args.executable())
            << "Cannot open mesh description file\n    "
            << meshDictIO.objectPath()
            << nl
            << exit(FatalError);
    }

    Info<< "Creating block mesh from\n    "
        << meshDictIO.objectPath() << endl;

    blockMesh::verbose(true);

    IOdictionary meshDict(meshDictIO);
    blockMesh blocks(meshDict, regionName);


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

            forAll(cellCentres, cellI)
            {
                //point cc = b.blockShape().centre(b.points());
                const point& cc = cellCentres[cellI];

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
            runTime
        ),
        xferCopy(blocks.points()),           // could we re-use space?
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
        List<Pair<word> > mergePatchPairs
        (
            meshDict.lookup("mergePatchPairs")
        );

#       include "mergePatchPairs.H"
    }
    else
    {
        Info<< nl << "There are no merge patch pairs edges" << endl;
    }


    // Set any cellZones (note: cell labelling unaffected by above
    // mergePatchPairs)

    label nZones = blocks.numZonedBlocks();

    if (nZones > 0)
    {
        Info<< nl << "Adding cell zones" << endl;

        // Map from zoneName to cellZone index
        HashTable<label> zoneMap(nZones);

        // Cells per zone.
        List<DynamicList<label> > zoneCells(nZones);

        // Running cell counter
        label cellI = 0;

        // Largest zone so far
        label freeZoneI = 0;

        forAll(blocks, blockI)
        {
            const block& b = blocks[blockI];
            const labelListList& blockCells = b.cells();
            const word& zoneName = b.blockDef().zoneName();

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
                    zoneCells[zoneI].append(cellI++);
                }
            }
            else
            {
                cellI += b.cells().size();
            }
        }


        List<cellZone*> cz(zoneMap.size());

        Info<< nl << "Writing cell zones as cellSets" << endl;

        forAllConstIter(HashTable<label>, zoneMap, iter)
        {
            label zoneI = iter();

            cz[zoneI] = new cellZone
            (
                iter.key(),
                zoneCells[zoneI].shrink(),
                zoneI,
                mesh.cellZones()
            );

            // Write as cellSet for ease of processing
            cellSet cset(mesh, iter.key(), zoneCells[zoneI].shrink());
            cset.write();
        }

        mesh.pointZones().setSize(0);
        mesh.faceZones().setSize(0);
        mesh.cellZones().setSize(0);
        mesh.addZones(List<pointZone*>(0), List<faceZone*>(0), cz);
    }

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< nl << "Writing polyMesh" << endl;
    mesh.removeFiles();
    if (!mesh.write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing polyMesh."
            << exit(FatalError);
    }


    //
    // write some information
    //
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

        forAll(patches, patchI)
        {
            const polyPatch& p = patches[patchI];

            Info<< "  " << "patch " << patchI
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
