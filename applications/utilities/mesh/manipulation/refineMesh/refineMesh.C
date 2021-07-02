/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
#include "polyMesh.H"
#include "Time.H"
#include "cellSet.H"
#include "multiDirRefinement.H"
#include "labelIOList.H"
#include "IOdictionary.H"
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


    Info<< "Mesh edge statistics:" << nl
        << "    x aligned :  number:" << nX << "\tminLen:" << minX
        << "\tmaxLen:" << maxX << nl
        << "    y aligned :  number:" << nY << "\tminLen:" << minY
        << "\tmaxLen:" << maxY << nl
        << "    z aligned :  number:" << nZ << "\tminLen:" << minZ
        << "\tmaxLen:" << maxZ << nl
        << "    other     :  number:" << nEdges - nX - nY - nZ
        << "\tminLen:" << minOther
        << "\tmaxLen:" << maxOther << nl << endl;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "refine cells in multiple directions"
    );

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"
    #include "addDictOption.H"

    Foam::argList::addBoolOption
    (
        "all",
        "Refine all cells"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createNamedPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    printEdgeStats(mesh);

    //
    // Read/construct control dictionary
    //

    const bool readDict = args.optionFound("dict");
    const bool refineAllCells = args.optionFound("all");
    const bool overwrite = args.optionFound("overwrite");

    // List of cells to refine
    labelList refCells;

    // Dictionary to control refinement
    const word dictName("refineMeshDict");
    IOobject dictIO(systemDictIO(dictName, args, runTime));
    dictionary refineDict;
    if (readDict)
    {
        if (dictIO.typeHeaderOk<IOdictionary>(true))
        {
            Info<< "Refining according to " << dictIO.path() << nl << endl;
            refineDict = IOdictionary(dictIO);
        }
        else
        {
            FatalErrorInFunction
                << "Cannot open specified refinement dictionary "
                << dictIO.path() << exit(FatalError);
        }
    }
    else if (!refineAllCells)
    {
        if (dictIO.typeHeaderOk<IOdictionary>(true))
        {
            Info<< "Refining according to " << dictIO.path() << nl << endl;
            refineDict = IOdictionary(dictIO);
        }
        else
        {
            Info<< "Refinement dictionary " << dictIO.path() << " not found"
                << nl << endl;
        }
    }

    if (refineDict.size())
    {
        const word setName(refineDict.lookup("set"));

        cellSet cells(mesh, setName);

        Info<< "Read " << returnReduce(cells.size(), sumOp<label>())
            << " cells from cellSet "
            << cells.instance()/cells.local()/cells.name()
            << endl << endl;

        refCells = cells.toc();
    }
    else
    {
        Info<< "Refining all cells" << nl << endl;

        // Select all cells
        refCells = identity(mesh.nCells());

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

        dictionary coeffsDict;
        coeffsDict.add("e1", vector(1, 0, 0));
        coeffsDict.add("e2", vector(0, 1, 0));
        refineDict.add("globalCoeffs", coeffsDict);

        refineDict.add("geometricCut", "false");
        refineDict.add("writeMesh", "false");
    }


    string oldTimeName(runTime.timeName());

    if (!overwrite)
    {
        runTime++;
    }


    // Multi-directional refinement (does multiple iterations)
    multiDirRefinement multiRef(mesh, refCells, refineDict);


    // Write resulting mesh
    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    mesh.write();


    // Get list of cell splits.
    // (is for every cell in old mesh the cells they have been split into)
    const labelListList& oldToNew = multiRef.addedCells();


    // Create cellSet with added cells for easy inspection
    cellSet newCells(mesh, "refinedCells", refCells.size());

    forAll(oldToNew, oldCelli)
    {
        const labelList& added = oldToNew[oldCelli];

        forAll(added, i)
        {
            newCells.insert(added[i]);
        }
    }

    Info<< "Writing refined cells ("
        << returnReduce(newCells.size(), sumOp<label>())
        << ") to cellSet "
        << newCells.instance()/newCells.local()/newCells.name()
        << endl << endl;

    newCells.write();



    //
    // Invert cell split to construct map from new to old
    //

    labelIOList newToOld
    (
        IOobject
        (
            "cellMap",
            runTime.timeName(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.nCells()
    );
    newToOld.note() =
        "From cells in mesh at "
      + runTime.timeName()
      + " to cells in mesh at "
      + oldTimeName;


    forAll(oldToNew, oldCelli)
    {
        const labelList& added = oldToNew[oldCelli];

        if (added.size())
        {
            forAll(added, i)
            {
                newToOld[added[i]] = oldCelli;
            }
        }
        else
        {
            // Unrefined cell
            newToOld[oldCelli] = oldCelli;
        }
    }

    Info<< "Writing map from new to old cell to "
        << newToOld.localObjectPath() << nl << endl;

    newToOld.write();

    printEdgeStats(mesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
