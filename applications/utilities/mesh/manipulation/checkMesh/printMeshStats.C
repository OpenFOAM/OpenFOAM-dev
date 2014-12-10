#include "printMeshStats.H"
#include "polyMesh.H"
#include "globalMeshData.H"

#include "hexMatcher.H"
#include "wedgeMatcher.H"
#include "prismMatcher.H"
#include "pyrMatcher.H"
#include "tetWedgeMatcher.H"
#include "tetMatcher.H"
#include "IOmanip.H"


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
            WarningIn("Foam::printMeshStats(const polyMesh&, const bool)")
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

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        if (hex.isA(mesh, cellI))
        {
            nHex++;
        }
        else if (tet.isA(mesh, cellI))
        {
            nTet++;
        }
        else if (pyr.isA(mesh, cellI))
        {
            nPyr++;
        }
        else if (prism.isA(mesh, cellI))
        {
            nPrism++;
        }
        else if (wedge.isA(mesh, cellI))
        {
            nWedge++;
        }
        else if (tetWedge.isA(mesh, cellI))
        {
            nTetWedge++;
        }
        else
        {
            nUnknown++;
            polyhedralFaces(mesh.cells()[cellI].size())++;
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
