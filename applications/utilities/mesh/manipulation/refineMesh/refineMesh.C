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
    refineMesh

Description
    Utility to refine cells in multiple directions.

    Either supply -all option to refine all cells (3D refinement for 3D
    cases; 2D for 2D cases) or reads a refineMeshDict with
    - cellSet to refine
    - directions to refine

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "undoableMeshCutter.H"
#include "hexCellLooper.H"
#include "cellSet.H"
#include "twoDPointCorrector.H"
#include "directions.H"
#include "OFstream.H"
#include "multiDirRefinement.H"
#include "labelIOList.H"
#include "wedgePolyPatch.H"
#include "plane.H"
#include "SubField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Max cos angle for edges to be considered aligned with axis.
static const scalar edgeTol = 1e-3;


// Calculate some edge statistics on mesh.
void printEdgeStats(const primitiveMesh& mesh)
{
    label nX = 0;
    label nY = 0;
    label nZ = 0;

    scalar minX = GREAT;
    scalar maxX = -GREAT;
    vector x(1, 0, 0);

    scalar minY = GREAT;
    scalar maxY = -GREAT;
    vector y(0, 1, 0);

    scalar minZ = GREAT;
    scalar maxZ = -GREAT;
    vector z(0, 0, 1);

    scalar minOther = GREAT;
    scalar maxOther = -GREAT;

    const edgeList& edges = mesh.edges();

    forAll(edges, edgeI)
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

    Pout<< "Mesh edge statistics:" << endl
        << "    x aligned :  number:" << nX << "\tminLen:" << minX
        << "\tmaxLen:" << maxX << endl
        << "    y aligned :  number:" << nY << "\tminLen:" << minY
        << "\tmaxLen:" << maxY << endl
        << "    z aligned :  number:" << nZ << "\tminLen:" << minZ
        << "\tmaxLen:" << maxZ << endl
        << "    other     :  number:" << mesh.nEdges() - nX - nY - nZ
        << "\tminLen:" << minOther
        << "\tmaxLen:" << maxOther << endl << endl;
}


// Return index of coordinate axis.
label axis(const vector& normal)
{
    label axisIndex = -1;

    if (mag(normal & point(1, 0, 0)) > (1-edgeTol))
    {
        axisIndex = 0;
    }
    else if (mag(normal & point(0, 1, 0)) > (1-edgeTol))
    {
        axisIndex = 1;
    }
    else if (mag(normal & point(0, 0, 1)) > (1-edgeTol))
    {
        axisIndex = 2;
    }

    return axisIndex;
}


//- Returns -1 or cartesian coordinate component (0=x, 1=y, 2=z) of normal
//  in case of 2D mesh
label twoDNess(const polyMesh& mesh)
{
    const pointField& ctrs = mesh.cellCentres();

    if (ctrs.size() < 2)
    {
        return -1;
    }


    //
    // 1. All cell centres on single plane aligned with x, y or z
    //

    // Determine 3 points to base plane on.
    vector vec10 = ctrs[1] - ctrs[0];
    vec10 /= mag(vec10);

    label otherCellI = -1;

    for (label cellI = 2; cellI < ctrs.size(); cellI++)
    {
        vector vec(ctrs[cellI] - ctrs[0]);
        vec /= mag(vec);

        if (mag(vec & vec10) < 0.9)
        {
            // ctrs[cellI] not in line with n
            otherCellI = cellI;

            break;
        }
    }

    if (otherCellI == -1)
    {
        // Cannot find cell to make decent angle with cell0-cell1 vector.
        // Note: what to do here? All cells (almost) in one line. Maybe 1D case?
        return -1;
    }

    plane cellPlane(ctrs[0], ctrs[1], ctrs[otherCellI]);


    forAll(ctrs, cellI)
    {
        const labelList& cEdges = mesh.cellEdges()[cellI];

        scalar minLen = GREAT;

        forAll(cEdges, i)
        {
            minLen = min(minLen, mesh.edges()[cEdges[i]].mag(mesh.points()));
        }

        if (cellPlane.distance(ctrs[cellI]) > 1e-6*minLen)
        {
            // Centres not in plane
            return  -1;
        }
    }

    label axisIndex = axis(cellPlane.normal());

    if (axisIndex == -1)
    {
        return axisIndex;
    }


    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    //
    // 2. No edges without points on boundary
    //

    // Mark boundary points
    boolList boundaryPoint(mesh.points().size(), false);

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        forAll(patch, patchFaceI)
        {
            const face& f = patch[patchFaceI];

            forAll(f, fp)
            {
                boundaryPoint[f[fp]] = true;
            }
        }
    }


    const edgeList& edges = mesh.edges();

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        if (!boundaryPoint[e.start()] && !boundaryPoint[e.end()])
        {
            // Edge has no point on boundary.
            return -1;
        }
    }


    // 3. For all non-wedge patches: all faces either perp or aligned with
    //    cell-plane normal. (wedge patches already checked upon construction)

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        if (!isA<wedgePolyPatch>(patch))
        {
            const vectorField& n = patch.faceAreas();

            const scalarField cosAngle(mag(n/mag(n) & cellPlane.normal()));

            if (mag(min(cosAngle) - max(cosAngle)) > 1e-6)
            {
                // cosAngle should be either ~1 over all faces (2D front and
                // back) or ~0 (all other patches perp to 2D)
                return -1;
            }
        }
    }

    return axisIndex;
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
    const bool overwrite = args.optionFound("overwrite");

    // List of cells to refine
    labelList refCells;

    // Dictionary to control refinement
    dictionary refineDict;

    if (readDict)
    {
        const word dictName("refineMeshDict");
        #include "setSystemMeshDictionaryIO.H"

        Info<< "Refining according to " << dictName << nl << endl;

        refineDict = IOdictionary(dictIO);

        const word setName(refineDict.lookup("set"));

        cellSet cells(mesh, setName);

        Pout<< "Read " << cells.size() << " cells from cellSet "
            << cells.instance()/cells.local()/cells.name()
            << endl << endl;

        refCells = cells.toc();
    }
    else
    {
        Info<< "Refining all cells" << nl << endl;

        // Select all cells
        refCells.setSize(mesh.nCells());

        forAll(mesh.cells(), cellI)
        {
            refCells[cellI] = cellI;
        }


        // Set refinement directions based on 2D/3D
        label axisIndex = twoDNess(mesh);

        if (axisIndex == -1)
        {
            Info<< "3D case; refining all directions" << nl << endl;

            wordList directions(3);
            directions[0] = "tan1";
            directions[1] = "tan2";
            directions[2] = "normal";
            refineDict.add("directions", directions);

            // Use hex cutter
            refineDict.add("useHexTopology", "true");
        }
        else
        {
            wordList directions(2);

            if (axisIndex == 0)
            {
                Info<< "2D case; refining in directions y,z\n" << endl;
                directions[0] = "tan2";
                directions[1] = "normal";
            }
            else if (axisIndex == 1)
            {
                Info<< "2D case; refining in directions x,z\n" << endl;
                directions[0] = "tan1";
                directions[1] = "normal";
            }
            else
            {
                Info<< "2D case; refining in directions x,y\n" << endl;
                directions[0] = "tan1";
                directions[1] = "tan2";
            }

            refineDict.add("directions", directions);

            // Use standard cutter
            refineDict.add("useHexTopology", "false");
        }

        refineDict.add("coordinateSystem", "global");

        dictionary coeffsDict;
        coeffsDict.add("tan1", vector(1, 0, 0));
        coeffsDict.add("tan2", vector(0, 1, 0));
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

    forAll(oldToNew, oldCellI)
    {
        const labelList& added = oldToNew[oldCellI];

        forAll(added, i)
        {
            newCells.insert(added[i]);
        }
    }

    Pout<< "Writing refined cells (" << newCells.size() << ") to cellSet "
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


    forAll(oldToNew, oldCellI)
    {
        const labelList& added = oldToNew[oldCellI];

        if (added.size())
        {
            forAll(added, i)
            {
                newToOld[added[i]] = oldCellI;
            }
        }
        else
        {
            // Unrefined cell
            newToOld[oldCellI] = oldCellI;
        }
    }

    Info<< "Writing map from new to old cell to "
        << newToOld.objectPath() << nl << endl;

    newToOld.write();


    // Some statistics.

    printEdgeStats(mesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
