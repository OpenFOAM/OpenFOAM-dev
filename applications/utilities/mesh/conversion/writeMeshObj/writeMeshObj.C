/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

Application
    writeMeshObj

Description
    For mesh debugging: writes mesh as three separate OBJ files which can
    be viewed with e.g. javaview.

    meshPoints_XXX.obj : all points and edges as lines.
    meshFaceCentres_XXX.obj : all face centres.
    meshCellCentres_XXX.obj : all cell centres.

    patch_YYY_XXX.obj : all face centres of patch YYY

    Optional: - patch faces (as polygons) : patchFaces_YYY_XXX.obj
              - non-manifold edges : patchEdges_YYY_XXX.obj

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "polyMesh.H"
#include "OFstream.H"
#include "meshTools.H"
#include "cellSet.H"
#include "faceSet.H"
#include "SubField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeOBJ(const point& pt, Ostream& os)
{
    os  << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
}

// All edges of mesh
void writePoints(const polyMesh& mesh, const fileName& timeName)
{
    label vertI = 0;

    fileName pointFile(mesh.time().path()/"meshPoints_" + timeName + ".obj");

    Info<< "Writing mesh points and edges to " << pointFile << endl;

    OFstream pointStream(pointFile);

    forAll(mesh.points(), pointi)
    {
        writeOBJ(mesh.points()[pointi], pointStream);
        vertI++;
    }

    forAll(mesh.edges(), edgeI)
    {
        const edge& e = mesh.edges()[edgeI];

        pointStream << "l " << e.start() + 1 << ' ' << e.end() + 1 << nl;
    }
}


// Edges for subset of cells
void writePoints
(
    const polyMesh& mesh,
    const labelList& cellLabels,
    const fileName& timeName
)
{
    fileName fName(mesh.time().path()/"meshPoints_" + timeName + ".obj");

    Info<< "Writing mesh points and edges to " << fName << endl;

    OFstream str(fName);

    // OBJ file vertex
    label vertI = 0;

    // From point to OBJ file vertex
    Map<label> pointToObj(6*cellLabels.size());

    forAll(cellLabels, i)
    {
        label celli = cellLabels[i];

        const labelList& cEdges = mesh.cellEdges()[celli];

        forAll(cEdges, cEdgeI)
        {
            const edge& e = mesh.edges()[cEdges[cEdgeI]];

            label v0;

            Map<label>::iterator e0Fnd = pointToObj.find(e[0]);

            if (e0Fnd == pointToObj.end())
            {
                meshTools::writeOBJ(str, mesh.points()[e[0]]);
                v0 = vertI++;
                pointToObj.insert(e[0], v0);
            }
            else
            {
                v0 = e0Fnd();
            }

            label v1;

            Map<label>::iterator e1Fnd = pointToObj.find(e[1]);

            if (e1Fnd == pointToObj.end())
            {
                meshTools::writeOBJ(str, mesh.points()[e[1]]);
                v1 = vertI++;
                pointToObj.insert(e[1], v1);
            }
            else
            {
                v1 = e1Fnd();
            }


            str << "l " << v0+1 << ' ' << v1+1 << nl;
        }
    }
}


// Edges of single cell
void writePoints
(
    const polyMesh& mesh,
    const label celli,
    const fileName& timeName
)
{
    fileName fName
    (
        mesh.time().path()
      / "meshPoints_" + timeName + '_' + name(celli) + ".obj"
    );

    Info<< "Writing mesh points and edges to " << fName << endl;

    OFstream pointStream(fName);

    const cell& cFaces = mesh.cells()[celli];

    meshTools::writeOBJ(pointStream, mesh.faces(), mesh.points(), cFaces);
}



// All face centres
void writeFaceCentres(const polyMesh& mesh,const fileName& timeName)
{
    fileName faceFile
    (
        mesh.time().path()
      / "meshFaceCentres_" + timeName + ".obj"
    );

    Info<< "Writing mesh face centres to " << faceFile << endl;

    OFstream faceStream(faceFile);

    forAll(mesh.faceCentres(), facei)
    {
        writeOBJ(mesh.faceCentres()[facei], faceStream);
    }
}


void writeCellCentres(const polyMesh& mesh, const fileName& timeName)
{
    fileName cellFile
    (
        mesh.time().path()/"meshCellCentres_" + timeName + ".obj"
    );

    Info<< "Writing mesh cell centres to " << cellFile << endl;

    OFstream cellStream(cellFile);

    forAll(mesh.cellCentres(), celli)
    {
        writeOBJ(mesh.cellCentres()[celli], cellStream);
    }
}


void writePatchCentres
(
    const polyMesh& mesh,
    const fileName& timeName
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        fileName faceFile
        (
            mesh.time().path()/"patch_" + pp.name() + '_' + timeName + ".obj"
        );

        Info<< "Writing patch face centres to " << faceFile << endl;

        OFstream patchFaceStream(faceFile);

        forAll(pp.faceCentres(), facei)
        {
            writeOBJ(pp.faceCentres()[facei], patchFaceStream);
        }
    }
}


void writePatchFaces
(
    const polyMesh& mesh,
    const fileName& timeName
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        fileName faceFile
        (
            mesh.time().path()
          / "patchFaces_" + pp.name() + '_' + timeName + ".obj"
        );

        Info<< "Writing patch faces to " << faceFile << endl;

        OFstream patchFaceStream(faceFile);

        forAll(pp.localPoints(), pointi)
        {
            writeOBJ(pp.localPoints()[pointi], patchFaceStream);
        }

        forAll(pp.localFaces(), facei)
        {
            const face& f = pp.localFaces()[facei];

            patchFaceStream<< 'f';

            forAll(f, fp)
            {
                patchFaceStream << ' ' << f[fp]+1;
            }
            patchFaceStream << nl;
        }
    }
}


void writePatchBoundaryEdges
(
    const polyMesh& mesh,
    const fileName& timeName
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        fileName edgeFile
        (
            mesh.time().path()
          / "patchEdges_" + pp.name() + '_' + timeName + ".obj"
        );

        Info<< "Writing patch edges to " << edgeFile << endl;

        OFstream patchEdgeStream(edgeFile);

        forAll(pp.localPoints(), pointi)
        {
            writeOBJ(pp.localPoints()[pointi], patchEdgeStream);
        }

        for (label edgeI = pp.nInternalEdges(); edgeI < pp.nEdges(); edgeI++)
        {
            if (pp.edgeFaces()[edgeI].size() == 1)
            {
                const edge& e = pp.edges()[edgeI];

                patchEdgeStream<< "l " << e[0]+1 << ' ' << e[1]+1 << nl;
            }
        }
    }
}


void writePointCells
(
    const polyMesh& mesh,
    const label pointi,
    const fileName& timeName
)
{
    const labelList& pCells = mesh.pointCells()[pointi];

    labelHashSet allEdges(6*pCells.size());

    forAll(pCells, i)
    {
        const labelList& cEdges = mesh.cellEdges()[pCells[i]];

        forAll(cEdges, i)
        {
            allEdges.insert(cEdges[i]);
        }
    }


    fileName pFile
    (
        mesh.time().path()
      / "pointEdges_" + timeName + '_' + name(pointi) + ".obj"
    );

    Info<< "Writing pointEdges to " << pFile << endl;

    OFstream pointStream(pFile);

    label vertI = 0;

    forAllConstIter(labelHashSet, allEdges, iter)
    {
        const edge& e = mesh.edges()[iter.key()];

        meshTools::writeOBJ(pointStream, mesh.points()[e[0]]); vertI++;
        meshTools::writeOBJ(pointStream, mesh.points()[e[1]]); vertI++;
        pointStream<< "l " << vertI-1 << ' ' << vertI << nl;
    }
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "for mesh debugging: write mesh as separate OBJ files"
    );

    timeSelector::addOptions();
    argList::addBoolOption
    (
        "patchFaces",
        "write patch faces edges"
    );
    argList::addBoolOption
    (
        "patchEdges",
        "write patch boundary edges"
    );
    argList::addOption
    (
        "cell",
        "int",
        "write points for the specified cell"
    );
    argList::addOption
    (
        "face",
        "int",
        "write specified face"
    );
    argList::addOption
    (
        "point",
        "int",
        "write specified point"
    );
    argList::addOption
    (
        "cellSet",
        "name",
        "write points for specified cellSet"
    );
    argList::addOption
    (
        "faceSet",
        "name",
        "write points for specified faceSet"
    );
    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    const bool patchFaces = args.optionFound("patchFaces");
    const bool patchEdges = args.optionFound("patchEdges");
    const bool doCell     = args.optionFound("cell");
    const bool doPoint    = args.optionFound("point");
    const bool doFace     = args.optionFound("face");
    const bool doCellSet  = args.optionFound("cellSet");
    const bool doFaceSet  = args.optionFound("faceSet");


    Info<< "Writing mesh objects as .obj files such that the object"
        << " numbering" << endl
        << "(for points, faces, cells) is consistent with"
        << " Foam numbering (starting from 0)." << endl << endl;

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedPolyMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        polyMesh::readUpdateState state = mesh.readUpdate();

        if (!timeI || state != polyMesh::UNCHANGED)
        {
            if (patchFaces)
            {
                writePatchFaces(mesh, runTime.timeName());
            }
            if (patchEdges)
            {
                writePatchBoundaryEdges(mesh, runTime.timeName());
            }
            if (doCell)
            {
                label celli = args.optionRead<label>("cell");

                writePoints(mesh, celli, runTime.timeName());
            }
            if (doPoint)
            {
                label pointi = args.optionRead<label>("point");

                writePointCells(mesh, pointi, runTime.timeName());
            }
            if (doFace)
            {
                label facei = args.optionRead<label>("face");

                fileName fName
                (
                    mesh.time().path()
                  / "meshPoints_"
                  + runTime.timeName()
                  + '_'
                  + name(facei)
                  + ".obj"
                );

                Info<< "Writing mesh points and edges to " << fName << endl;

                OFstream str(fName);

                const face& f = mesh.faces()[facei];

                meshTools::writeOBJ(str, faceList(1, f), mesh.points());
            }
            if (doCellSet)
            {
                const word setName = args["cellSet"];

                cellSet cells(mesh, setName);

                Info<< "Read " << cells.size() << " cells from set " << setName
                    << endl;

                writePoints(mesh, cells.toc(), runTime.timeName());
            }
            if (doFaceSet)
            {
                const word setName = args["faceSet"];

                faceSet faces(mesh, setName);

                Info<< "Read " << faces.size() << " faces from set " << setName
                    << endl;

                fileName fName
                (
                    mesh.time().path()
                  / "meshPoints_"
                  + runTime.timeName()
                  + '_'
                  + setName
                  + ".obj"
                );

                Info<< "Writing mesh points and edges to " << fName << endl;

                OFstream str(fName);

                meshTools::writeOBJ
                (
                    str,
                    mesh.faces(),
                    mesh.points(),
                    faces.toc()
                );
            }
            else if
            (
                !patchFaces
             && !patchEdges
             && !doCell
             && !doPoint
             && !doFace
             && !doCellSet
             && !doFaceSet
            )
            {
                // points & edges
                writePoints(mesh, runTime.timeName());

                // face centres
                writeFaceCentres(mesh, runTime.timeName());

                // cell centres
                writeCellCentres(mesh, runTime.timeName());

                // Patch face centres
                writePatchCentres(mesh, runTime.timeName());
            }
        }
        else
        {
            Info<< "No mesh." << endl;
        }

        Info<< nl << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
