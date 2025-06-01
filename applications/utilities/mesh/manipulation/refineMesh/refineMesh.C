/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    refineMesh

Description
    Utility to refine cells in multiple directions.

    Command-line option handling:
    - If -all specified or no refineMeshDict exists or, refine all cells
    - If -dict \<file\> specified refine according to \<file\>
    - If refineMeshDict exists refine according to refineMeshDict

    When the refinement or all cells is selected apply 3D refinement for 3D
    cases and 2D refinement for 2D cases.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "zoneGenerator.H"
#include "cellSet.H"
#include "hexRef8.H"
#include "multiDirRefinement.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "systemDict.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Max cos angle for edges to be considered aligned with axis.
static const scalar edgeTol = 1e-3;

// Print edge statistics on mesh.
void printEdgeStats(const polyMesh& mesh)
{
    label nX = 0;
    label nY = 0;
    label nZ = 0;

    scalar minX = great;
    scalar maxX = -great;
    static const vector x(1, 0, 0);

    scalar minY = great;
    scalar maxY = -great;
    static const vector y(0, 1, 0);

    scalar minZ = great;
    scalar maxZ = -great;
    static const vector z(0, 0, 1);

    scalar minOther = great;
    scalar maxOther = -great;

    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));

    const edgeList& edges = mesh.edges();

    forAll(edges, edgeI)
    {
        if (isMasterEdge[edgeI])
        {
            const edge& e = edges[edgeI];

            vector eVec(e.vec(mesh.points()));

            scalar eMag = mag(eVec);

            eVec /= eMag;

            if (mag(eVec & x) > 1-edgeTol)
            {
                minX = min(minX, eMag);
                maxX = max(maxX, eMag);
                nX++;
            }
            else if (mag(eVec & y) > 1-edgeTol)
            {
                minY = min(minY, eMag);
                maxY = max(maxY, eMag);
                nY++;
            }
            else if (mag(eVec & z) > 1-edgeTol)
            {
                minZ = min(minZ, eMag);
                maxZ = max(maxZ, eMag);
                nZ++;
            }
            else
            {
                minOther = min(minOther, eMag);
                maxOther = max(maxOther, eMag);
            }
        }
    }

    label nEdges = mesh.nEdges();
    reduce(nEdges, sumOp<label>());
    reduce(nX, sumOp<label>());
    reduce(nY, sumOp<label>());
    reduce(nZ, sumOp<label>());

    reduce(minX, minOp<scalar>());
    reduce(maxX, maxOp<scalar>());

    reduce(minY, minOp<scalar>());
    reduce(maxY, maxOp<scalar>());

    reduce(minZ, minOp<scalar>());
    reduce(maxZ, maxOp<scalar>());

    reduce(minOther, minOp<scalar>());
    reduce(maxOther, maxOp<scalar>());

    Info<< "Mesh edge statistics:" << nl;
    if (minX > great || maxX > -great)
    {
        Info<< "    x aligned :  number:" << nX << "\tminLen:" << minX
            << "\tmaxLen:" << maxX << nl;
    }
    if (minY > great || maxY > -great)
    {
        Info<< "    y aligned :  number:" << nY << "\tminLen:" << minY
            << "\tmayLen:" << maxY << nl;
    }
    if (minZ > great || maxZ > -great)
    {
        Info<< "    z aligned :  number:" << nZ << "\tminLen:" << minZ
            << "\tmazLen:" << maxZ << nl;
    }
    if (minOther > great || maxOther > -great)
    {
        Info<< "    other     :  number:" << nEdges - nX - nY - nZ
            << "\tminLen:" << minOther
            << "\tmaxLen:" << maxOther << nl;
    }

    Info<< endl;
}


void refineZone
(
    polyMesh& mesh,
    const labelList& refCells,
    autoPtr<hexRef8>& meshCutterPtr,
    const dictionary& refineDict,
    const dictionary& zoneDict
)
{
    const label nCells0 = mesh.globalData().nTotalCells();

    if (meshCutterPtr.valid())
    {
        hexRef8& meshCutter = meshCutterPtr();

        // Maintain 2:1 ratio
        const labelList newCellsToRefine
        (
            meshCutter.consistentRefinement
            (
                refCells,
                true                 // extend set
            )
        );

        // Mesh changing engine.
        polyTopoChange meshMod(mesh);

        // Play refinement commands into mesh changer.
        meshCutter.setRefinement(newCellsToRefine, meshMod);

        // Create mesh, return map from old to new mesh.
        autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh);

        // Update mesh objects
        mesh.topoChange(map);

        // Update hexRef8 cells/vertices
        meshCutter.topoChange(map);
    }
    else
    {
        multiDirRefinement
        (
            mesh,
            refCells,
            refineDict,
            zoneDict.isDict("coordinates")
          ? zoneDict.subDict("coordinates")
          : refineDict.subDict("coordinates")
        );
    }

    Info<< "    total number of cell increased from " << nCells0
        << " to " << mesh.globalData().nTotalCells()
        << nl << endl;
}


int main(int argc, char *argv[])
{
    argList::addNote("Refine cells in multiple directions");

    #include "addOverwriteOption.H"
    #include "addMeshOption.H"
    #include "addRegionOption.H"
    #include "addDictOption.H"

    Foam::argList::addBoolOption
    (
        "all",
        "Refine all cells"
    );

    #include "setRootCase.H"
    #include "createTimeNoFunctionObjects.H"
    #include "createSpecifiedPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    // Some stats
    Info<< "Read mesh:" << nl
        << "    cells:" << mesh.globalData().nTotalCells() << nl
        << "    faces:" << mesh.globalData().nTotalFaces() << nl
        << "    points:" << mesh.globalData().nTotalPoints() << nl << endl;
    printEdgeStats(mesh);

    const bool refineAllCells = args.optionFound("all");
    const bool overwrite = args.optionFound("overwrite");

    autoPtr<hexRef8> meshCutterPtr;

    if (refineAllCells)
    {
        Info<< "Refining all cells" << nl << endl;

        // Select all cells
        labelList refCells(identityMap(mesh.nCells()));

        dictionary refineDict;

        if (mesh.nGeometricD() == 3)
        {
            Info<< "3D case; refining all directions" << nl << endl;

            wordList directions(3);
            directions[0] = "e1";
            directions[1] = "e2";
            directions[2] = "e3";
            refineDict.add("directions", directions);

            // Use hex cutter
            refineDict.add("useHexTopology", "true");
        }
        else
        {
            const Vector<label> dirs(mesh.geometricD());

            wordList directions(2);

            if (dirs.x() == -1)
            {
                Info<< "2D case; refining in directions y,z\n" << endl;
                directions[0] = "e2";
                directions[1] = "e3";
            }
            else if (dirs.y() == -1)
            {
                Info<< "2D case; refining in directions x,z\n" << endl;
                directions[0] = "e1";
                directions[1] = "e3";
            }
            else
            {
                Info<< "2D case; refining in directions x,y\n" << endl;
                directions[0] = "e1";
                directions[1] = "e2";
            }

            refineDict.add("directions", directions);

            // Use standard cutter
            refineDict.add("useHexTopology", "false");
        }

        refineDict.add("coordinateSystem", "global");

        dictionary coeffDict;
        coeffDict.add("e1", vector(1, 0, 0));
        coeffDict.add("e2", vector(0, 1, 0));
        refineDict.add("globalCoeffs", coeffDict);

        refineDict.add("geometricCut", "false");
        refineDict.add("writeMesh", "false");

        multiDirRefinement multiRef(mesh, refCells, refineDict, refineDict);
    }
    else
    {
        const dictionary refineDict
        (
            systemDict("refineMeshDict", args, mesh)
        );

        const bool hexRef8Refine(refineDict.lookupOrDefault("hexRef8", false));

        // Construct refiner without unrefinement.
        // Read existing point/cell level.
        if (hexRef8Refine)
        {
            meshCutterPtr = new hexRef8(mesh);
            hexRef8& meshCutter = meshCutterPtr();

            Info<< "Refining using hexRef8" << nl
                << "    cellLevel :"
                << " min:" << gMin(meshCutter.cellLevel())
                << " max:" << gMax(meshCutter.cellLevel()) << nl
                << "    pointLevel :"
                << " min:" << gMin(meshCutter.pointLevel())
                << " max:" << gMax(meshCutter.pointLevel()) << nl
                << endl;
        }

        if (refineDict.found("set"))
        {
            const word setName(refineDict.lookup("set"));
            const cellSet cells(mesh, setName);

            Info<< "Refining "
                << returnReduce(cells.size(), sumOp<label>())
                << " cells in set "
                << cells.instance()/cells.local()/cells.name()
                << endl;

            refineZone
            (
                mesh,
                cells.toc(),
                meshCutterPtr,
                refineDict,
                refineDict
            );
        }
        else if (refineDict.found("zone"))
        {
            labelList refCells;

            if (refineDict.isDict("zone"))
            {
                autoPtr<zoneGenerator> zg
                (
                    zoneGenerator::New
                    (
                        "zone",
                        zoneGenerator::zoneTypes::cell,
                        mesh,
                        refineDict.subDict("zone")
                    )
                );

                refCells = zg->generate().cZone();

                Info<< "Refining "
                    << returnReduce(refCells.size(), sumOp<label>())
                    << " cells in zone " << zg->zoneName()
                    << " of type " << zg->type() << endl;
            }
            else
            {
                const word cellZoneName(refineDict.lookup("zone"));
                refCells = mesh.cellZones()[cellZoneName];

                Info<< "Refining "
                    << returnReduce(refCells.size(), sumOp<label>())
                    << " cells in zone " << cellZoneName << endl;
            }

            refineZone
            (
                mesh,
                refCells,
                meshCutterPtr,
                refineDict,
                refineDict.isDict("zone")
                  ? refineDict.subDict("zone")
                  : refineDict
            );
        }
        else if (refineDict.found("zones"))
        {
            const dictionary& zones = refineDict.subDict("zones");

            forAllConstIter(dictionary, zones, iter)
            {
                const word& name = iter().keyword();
                const dictionary& zoneDict = iter().dict();

                autoPtr<zoneGenerator> zg
                (
                    zoneGenerator::New
                    (
                        name,
                        zoneGenerator::zoneTypes::cell,
                        mesh,
                        zoneDict
                    )
                );

                const labelList refCells(zg->generate().cZone());

                Info<< "Refining "
                    << returnReduce(refCells.size(), sumOp<label>())
                    << " cells in zone " << zg->zoneName()
                    << " of type " << zg->type() << endl;

                refineZone
                (
                    mesh,
                    refCells,
                    meshCutterPtr,
                    refineDict,
                    zoneDict
                );
            }
        }
    }

    printEdgeStats(mesh);

    if (overwrite)
    {
        mesh.setInstance(oldInstance);

        if (meshCutterPtr.valid())
        {
            meshCutterPtr->setInstance(oldInstance);
        }
    }
    else
    {
        runTime++;
    }

    Info<< "Writing mesh to ";
    mesh.write();
    Info<< mesh.facesInstance()/mesh.meshDir() << endl;

    if (meshCutterPtr.valid())
    {
        Info<< "Writing hexRef8 refinement level files to "
            << mesh.facesInstance()/mesh.meshDir() << endl;
        meshCutterPtr->write();
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
