/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "checkTools.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "hexMatcher.H"
#include "wedgeMatcher.H"
#include "prismMatcher.H"
#include "pyrMatcher.H"
#include "tetWedgeMatcher.H"
#include "tetMatcher.H"
#include "IOmanip.H"
#include "pointSet.H"
#include "faceSet.H"
#include "cellSet.H"
#include "Time.H"
#include "surfaceWriter.H"
#include "syncTools.H"
#include "globalIndex.H"
#include "PatchTools.H"


void Foam::printMeshStats(const polyMesh& mesh, const bool allTopology)
{
    Info<< "Mesh stats" << nl
        << "    points:           "
        << returnReduce(mesh.points().size(), sumOp<label>()) << nl;

    label nInternalPoints = returnReduce
    (
        mesh.nInternalPoints(),
        sumOp<label>()
    );

    if (nInternalPoints != -Pstream::nProcs())
    {
        Info<< "    internal points:  " << nInternalPoints << nl;

        if (returnReduce(mesh.nInternalPoints(), minOp<label>()) == -1)
        {
            WarningInFunction
                << "Some processors have their points sorted into internal"
                << " and external and some do not." << endl
                << "This can cause problems later on." << endl;
        }
    }

    if (allTopology && nInternalPoints != -Pstream::nProcs())
    {
        label nEdges = returnReduce(mesh.nEdges(), sumOp<label>());
        label nInternalEdges = returnReduce
        (
            mesh.nInternalEdges(),
            sumOp<label>()
        );
        label nInternal1Edges = returnReduce
        (
            mesh.nInternal1Edges(),
            sumOp<label>()
        );
        label nInternal0Edges = returnReduce
        (
            mesh.nInternal0Edges(),
            sumOp<label>()
        );

        Info<< "    edges:            " << nEdges << nl
            << "    internal edges:   " << nInternalEdges << nl
            << "    internal edges using one boundary point:   "
            << nInternal1Edges-nInternal0Edges << nl
            << "    internal edges using two boundary points:  "
            << nInternalEdges-nInternal1Edges << nl;
    }

    label nFaces = returnReduce(mesh.faces().size(), sumOp<label>());
    label nIntFaces = returnReduce(mesh.faceNeighbour().size(), sumOp<label>());
    label nCells = returnReduce(mesh.cells().size(), sumOp<label>());

    Info<< "    faces:            " << nFaces << nl
        << "    internal faces:   " << nIntFaces << nl
        << "    cells:            " << nCells << nl
        << "    faces per cell:   "
        << scalar(nFaces + nIntFaces)/max(1, nCells) << nl
        << "    boundary patches: " << mesh.boundaryMesh().size() << nl
        << "    point zones:      " << mesh.pointZones().size() << nl
        << "    face zones:       " << mesh.faceZones().size() << nl
        << "    cell zones:       " << mesh.cellZones().size() << nl
        << endl;

    // Construct shape recognizers
    hexMatcher hex;
    prismMatcher prism;
    wedgeMatcher wedge;
    pyrMatcher pyr;
    tetWedgeMatcher tetWedge;
    tetMatcher tet;

    // Counters for different cell types
    label nHex = 0;
    label nWedge = 0;
    label nPrism = 0;
    label nPyr = 0;
    label nTet = 0;
    label nTetWedge = 0;
    label nUnknown = 0;

    Map<label> polyhedralFaces;

    for (label celli = 0; celli < mesh.nCells(); celli++)
    {
        if (hex.isA(mesh, celli))
        {
            nHex++;
        }
        else if (tet.isA(mesh, celli))
        {
            nTet++;
        }
        else if (pyr.isA(mesh, celli))
        {
            nPyr++;
        }
        else if (prism.isA(mesh, celli))
        {
            nPrism++;
        }
        else if (wedge.isA(mesh, celli))
        {
            nWedge++;
        }
        else if (tetWedge.isA(mesh, celli))
        {
            nTetWedge++;
        }
        else
        {
            nUnknown++;
            polyhedralFaces(mesh.cells()[celli].size())++;
        }
    }

    reduce(nHex,sumOp<label>());
    reduce(nPrism,sumOp<label>());
    reduce(nWedge,sumOp<label>());
    reduce(nPyr,sumOp<label>());
    reduce(nTetWedge,sumOp<label>());
    reduce(nTet,sumOp<label>());
    reduce(nUnknown,sumOp<label>());

    Info<< "Overall number of cells of each type:" << nl
        << "    hexahedra:     " << nHex << nl
        << "    prisms:        " << nPrism << nl
        << "    wedges:        " << nWedge << nl
        << "    pyramids:      " << nPyr << nl
        << "    tet wedges:    " << nTetWedge << nl
        << "    tetrahedra:    " << nTet << nl
        << "    polyhedra:     " << nUnknown
        << endl;

    if (nUnknown > 0)
    {
        Pstream::mapCombineGather(polyhedralFaces, plusEqOp<label>());

        Info<< "    Breakdown of polyhedra by number of faces:" << nl
            << "        faces" << "   number of cells" << endl;

        const labelList sortedKeys = polyhedralFaces.sortedToc();

        forAll(sortedKeys, keyI)
        {
            const label nFaces = sortedKeys[keyI];

            Info<< setf(std::ios::right) << setw(13)
                << nFaces << "   " << polyhedralFaces[nFaces] << nl;
        }
    }

    Info<< endl;
}


void Foam::mergeAndWrite
(
    const polyMesh& mesh,
    const surfaceWriter& writer,
    const word& name,
    const indirectPrimitivePatch setPatch,
    const fileName& outputDir
)
{
    if (Pstream::parRun())
    {
        labelList pointToGlobal;
        labelList uniqueMeshPointLabels;
        autoPtr<globalIndex> globalPoints;
        autoPtr<globalIndex> globalFaces;
        faceList mergedFaces;
        pointField mergedPoints;
        Foam::PatchTools::gatherAndMerge
        (
            mesh,
            setPatch.localFaces(),
            setPatch.meshPoints(),
            setPatch.meshPointMap(),

            pointToGlobal,
            uniqueMeshPointLabels,
            globalPoints,
            globalFaces,

            mergedFaces,
            mergedPoints
        );

        // Write
        if (Pstream::master())
        {
            writer.write(outputDir, name, mergedPoints, mergedFaces);
        }
    }
    else
    {
        writer.write
        (
            outputDir,
            name,
            setPatch.localPoints(),
            setPatch.localFaces()
        );
    }
}


void Foam::mergeAndWrite
(
    const surfaceWriter& writer,
    const faceSet& set
)
{
    const polyMesh& mesh = refCast<const polyMesh>(set.db());

    const indirectPrimitivePatch setPatch
    (
        IndirectList<face>(mesh.faces(), set.sortedToc()),
        mesh.points()
    );

    fileName outputDir
    (
        set.time().path()
      / (Pstream::parRun() ? ".." : "")
      / "postProcessing"
      / mesh.pointsInstance()
      / set.name()
    );
    outputDir.clean();

    mergeAndWrite(mesh, writer, set.name(), setPatch, outputDir);
}


void Foam::mergeAndWrite
(
    const surfaceWriter& writer,
    const cellSet& set
)
{
    const polyMesh& mesh = refCast<const polyMesh>(set.db());
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();


    // Determine faces on outside of cellSet
    PackedBoolList isInSet(mesh.nCells());
    forAllConstIter(cellSet, set, iter)
    {
        isInSet[iter.key()] = true;
    }


    boolList bndInSet(mesh.nFaces()-mesh.nInternalFaces());
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const labelList& fc = pp.faceCells();
        forAll(fc, i)
        {
            bndInSet[pp.start()+i-mesh.nInternalFaces()] = isInSet[fc[i]];
        }
    }
    syncTools::swapBoundaryFaceList(mesh, bndInSet);


    DynamicList<label> outsideFaces(3*set.size());
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        bool ownVal = isInSet[mesh.faceOwner()[facei]];
        bool neiVal = isInSet[mesh.faceNeighbour()[facei]];

        if (ownVal != neiVal)
        {
            outsideFaces.append(facei);
        }
    }


    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const labelList& fc = pp.faceCells();
        if (pp.coupled())
        {
            forAll(fc, i)
            {
                label facei = pp.start()+i;

                bool neiVal = bndInSet[facei-mesh.nInternalFaces()];
                if (isInSet[fc[i]] && !neiVal)
                {
                    outsideFaces.append(facei);
                }
            }
        }
        else
        {
            forAll(fc, i)
            {
                if (isInSet[fc[i]])
                {
                    outsideFaces.append(pp.start()+i);
                }
            }
        }
    }


    const indirectPrimitivePatch setPatch
    (
        IndirectList<face>(mesh.faces(), outsideFaces),
        mesh.points()
    );

    fileName outputDir
    (
        set.time().path()
      / (Pstream::parRun() ? ".." : "")
      / "postProcessing"
      / mesh.pointsInstance()
      / set.name()
    );
    outputDir.clean();

    mergeAndWrite(mesh, writer, set.name(), setPatch, outputDir);
}


void Foam::mergeAndWrite
(
    const writer<scalar>& writer,
    const pointSet& set
)
{
    const polyMesh& mesh = refCast<const polyMesh>(set.db());

    pointField mergedPts;
    labelList mergedIDs;

    if (Pstream::parRun())
    {
        // Note: we explicitly do not merge the points
        // (mesh.globalData().mergePoints etc) since this might
        // hide any synchronisation problem

        globalIndex globalNumbering(mesh.nPoints());

        mergedPts.setSize(returnReduce(set.size(), sumOp<label>()));
        mergedIDs.setSize(mergedPts.size());

        labelList setPointIDs(set.sortedToc());

        // Get renumbered local data
        pointField myPoints(mesh.points(), setPointIDs);
        labelList myIDs(setPointIDs.size());
        forAll(setPointIDs, i)
        {
            myIDs[i] = globalNumbering.toGlobal(setPointIDs[i]);
        }

        if (Pstream::master())
        {
            // Insert master data first
            label pOffset = 0;
            SubList<point>(mergedPts, myPoints.size(), pOffset) = myPoints;
            SubList<label>(mergedIDs, myIDs.size(), pOffset) = myIDs;
            pOffset += myPoints.size();

            // Receive slave ones
            for (int slave=1; slave<Pstream::nProcs(); slave++)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);

                pointField slavePts(fromSlave);
                labelList slaveIDs(fromSlave);

                SubList<point>(mergedPts, slavePts.size(), pOffset) = slavePts;
                SubList<label>(mergedIDs, slaveIDs.size(), pOffset) = slaveIDs;
                pOffset += slaveIDs.size();
            }
        }
        else
        {
            // Construct processor stream with estimate of size. Could
            // be improved.
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo(),
                myPoints.byteSize() + myIDs.byteSize()
            );
            toMaster << myPoints << myIDs;
        }
    }
    else
    {
        mergedIDs = set.sortedToc();
        mergedPts = pointField(mesh.points(), mergedIDs);
    }


    // Write with scalar pointID
    if (Pstream::master())
    {
        scalarField scalarPointIDs(mergedIDs.size());
        forAll(mergedIDs, i)
        {
            scalarPointIDs[i] = 1.0*mergedIDs[i];
        }

        coordSet points(set.name(), "distance", mergedPts, mag(mergedPts));

        List<const scalarField*> flds(1, &scalarPointIDs);

        wordList fldNames(1, "pointID");

        // Output e.g. pointSet p0 to
        // postProcessing/<time>/p0.vtk
        fileName outputDir
        (
            set.time().path()
          / (Pstream::parRun() ? ".." : "")
          / "postProcessing"
          / mesh.pointsInstance()
          // set.name()
        );
        outputDir.clean();
        mkDir(outputDir);

        fileName outputFile(outputDir/writer.getFileName(points, wordList()));
        // fileName outputFile(outputDir/set.name());

        OFstream os(outputFile);

        writer.write(points, fldNames, flds, os);
    }
}


// ************************************************************************* //
