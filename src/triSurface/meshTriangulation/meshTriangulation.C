/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "meshTriangulation.H"
#include "polyMesh.H"
#include "faceTriangulation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::meshTriangulation::isInternalFace
(
    const primitiveMesh& mesh,
    const boolList& includedCell,
    const label faceI
)
{
    if (mesh.isInternalFace(faceI))
    {
        label own = mesh.faceOwner()[faceI];
        label nei = mesh.faceNeighbour()[faceI];

        if (includedCell[own] && includedCell[nei])
        {
            // Neighbouring cell will get included in subset
            // as well so face is internal.
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


void Foam::meshTriangulation::getFaces
(
    const primitiveMesh& mesh,
    const boolList& includedCell,
    boolList& faceIsCut,
    label& nFaces,
    label& nInternalFaces
)
{
    // All faces to be triangulated.
    faceIsCut.setSize(mesh.nFaces());
    faceIsCut = false;

    nFaces = 0;
    nInternalFaces = 0;

    forAll(includedCell, cellI)
    {
        // Include faces of cut cells only.
        if (includedCell[cellI])
        {
            const labelList& cFaces = mesh.cells()[cellI];

            forAll(cFaces, i)
            {
                label faceI = cFaces[i];

                if (!faceIsCut[faceI])
                {
                    // First visit of face.
                    nFaces++;
                    faceIsCut[faceI] = true;

                    // See if would become internal or external face
                    if (isInternalFace(mesh, includedCell, faceI))
                    {
                        nInternalFaces++;
                    }
                }
            }
        }
    }

    Pout<< "Subset consists of " << nFaces << " faces out of " << mesh.nFaces()
        << " of which " << nInternalFaces << " are internal" << endl;
}


void Foam::meshTriangulation::insertTriangles
(
    const triFaceList& faceTris,
    const label faceI,
    const label regionI,
    const bool reverse,

    List<labelledTri>& triangles,
    label& triI
)
{
    // Copy triangles. Optionally reverse them
    forAll(faceTris, i)
    {
        const triFace& f = faceTris[i];

        labelledTri& tri = triangles[triI];

        if (reverse)
        {
            tri[0] = f[0];
            tri[2] = f[1];
            tri[1] = f[2];
        }
        else
        {
            tri[0] = f[0];
            tri[1] = f[1];
            tri[2] = f[2];
        }

        tri.region() = regionI;

        faceMap_[triI] = faceI;

        triI++;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
Foam::meshTriangulation::meshTriangulation()
:
    triSurface(),
    nInternalFaces_(0),
    faceMap_()
{}


// Construct from faces of cells
Foam::meshTriangulation::meshTriangulation
(
    const polyMesh& mesh,
    const label internalFacesPatch,
    const boolList& includedCell,
    const bool faceCentreDecomposition
)
:
    triSurface(),
    nInternalFaces_(0),
    faceMap_()
{
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // All faces to be triangulated.
    boolList faceIsCut;
    label nFaces, nInternalFaces;

    getFaces
    (
        mesh,
        includedCell,
        faceIsCut,
        nFaces,
        nInternalFaces
    );


    // Find upper limit for number of triangles
    // (can be less if triangulation fails)
    label nTotTri = 0;

    if (faceCentreDecomposition)
    {
        forAll(faceIsCut, faceI)
        {
            if (faceIsCut[faceI])
            {
                nTotTri += faces[faceI].size();
            }
        }
    }
    else
    {
        forAll(faceIsCut, faceI)
        {
            if (faceIsCut[faceI])
            {
                nTotTri += faces[faceI].nTriangles(points);
            }
        }
    }
    Pout<< "nTotTri : " << nTotTri << endl;


    // Storage for new and old points (only for faceCentre decomposition;
    // for triangulation uses only existing points)
    pointField newPoints;

    if (faceCentreDecomposition)
    {
        newPoints.setSize(mesh.nPoints() + faces.size());
        forAll(mesh.points(), pointI)
        {
            newPoints[pointI] = mesh.points()[pointI];
        }
        // Face centres
        forAll(faces, faceI)
        {
            newPoints[mesh.nPoints() + faceI] = mesh.faceCentres()[faceI];
        }
    }

    // Storage for all triangles
    List<labelledTri> triangles(nTotTri);
    faceMap_.setSize(nTotTri);
    label triI = 0;


    if (faceCentreDecomposition)
    {
        // Decomposition around face centre
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Triangulate internal faces
        forAll(faceIsCut, faceI)
        {
            if (faceIsCut[faceI] && isInternalFace(mesh, includedCell, faceI))
            {
                // Face was internal to the mesh and will be 'internal' to
                // the surface.

                // Triangulate face
                const face& f = faces[faceI];

                forAll(f, fp)
                {
                    faceMap_[triI] = faceI;

                    triangles[triI++] =
                        labelledTri
                        (
                            f[fp],
                            f.nextLabel(fp),
                            mesh.nPoints() + faceI,     // face centre
                            internalFacesPatch
                        );
                }
            }
        }
        nInternalFaces_ = triI;


        // Triangulate external faces
        forAll(faceIsCut, faceI)
        {
            if (faceIsCut[faceI] && !isInternalFace(mesh, includedCell, faceI))
            {
                // Face will become outside of the surface.

                label patchI = -1;
                bool reverse = false;

                if (mesh.isInternalFace(faceI))
                {
                    patchI = internalFacesPatch;

                    // Check orientation. Check which side of the face gets
                    // included (note: only one side is).
                    if (includedCell[mesh.faceOwner()[faceI]])
                    {
                        reverse = false;
                    }
                    else
                    {
                        reverse = true;
                    }
                }
                else
                {
                    // Face was already outside so orientation ok.

                    patchI = patches.whichPatch(faceI);

                    reverse = false;
                }


                // Triangulate face
                const face& f = faces[faceI];

                if (reverse)
                {
                    forAll(f, fp)
                    {
                        faceMap_[triI] = faceI;

                        triangles[triI++] =
                            labelledTri
                            (
                                f.nextLabel(fp),
                                f[fp],
                                mesh.nPoints() + faceI,     // face centre
                                patchI
                            );
                    }
                }
                else
                {
                    forAll(f, fp)
                    {
                        faceMap_[triI] = faceI;

                        triangles[triI++] =
                            labelledTri
                            (
                                f[fp],
                                f.nextLabel(fp),
                                mesh.nPoints() + faceI,     // face centre
                                patchI
                            );
                    }
                }
            }
        }
    }
    else
    {
        // Triangulation using existing vertices
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Triangulate internal faces
        forAll(faceIsCut, faceI)
        {
            if (faceIsCut[faceI] && isInternalFace(mesh, includedCell, faceI))
            {
                // Face was internal to the mesh and will be 'internal' to
                // the surface.

                // Triangulate face. Fall back to naive triangulation if failed.
                faceTriangulation faceTris(points, faces[faceI], true);

                if (faceTris.empty())
                {
                    WarningIn("meshTriangulation::meshTriangulation")
                        << "Could not find triangulation for face " << faceI
                        << " vertices " << faces[faceI] << " coords "
                        << IndirectList<point>(points, faces[faceI])() << endl;
                }
                else
                {
                    // Copy triangles. Make them internalFacesPatch
                    insertTriangles
                    (
                        faceTris,
                        faceI,
                        internalFacesPatch,
                        false,                  // no reverse

                        triangles,
                        triI
                    );
                }
            }
        }
        nInternalFaces_ = triI;


        // Triangulate external faces
        forAll(faceIsCut, faceI)
        {
            if (faceIsCut[faceI] && !isInternalFace(mesh, includedCell, faceI))
            {
                // Face will become outside of the surface.

                label patchI = -1;
                bool reverse = false;

                if (mesh.isInternalFace(faceI))
                {
                    patchI = internalFacesPatch;

                    // Check orientation. Check which side of the face gets
                    // included (note: only one side is).
                    if (includedCell[mesh.faceOwner()[faceI]])
                    {
                        reverse = false;
                    }
                    else
                    {
                        reverse = true;
                    }
                }
                else
                {
                    // Face was already outside so orientation ok.

                    patchI = patches.whichPatch(faceI);

                    reverse = false;
                }

                // Triangulate face
                faceTriangulation faceTris(points, faces[faceI], true);

                if (faceTris.empty())
                {
                    WarningIn("meshTriangulation::meshTriangulation")
                        << "Could not find triangulation for face " << faceI
                        << " vertices " << faces[faceI] << " coords "
                        << IndirectList<point>(points, faces[faceI])() << endl;
                }
                else
                {
                    // Copy triangles. Optionally reverse them
                    insertTriangles
                    (
                        faceTris,
                        faceI,
                        patchI,
                        reverse,    // whether to reverse

                        triangles,
                        triI
                    );
                }
            }
        }
    }

    // Shrink if necessary (because of invalid triangulations)
    triangles.setSize(triI);
    faceMap_.setSize(triI);

    Pout<< "nInternalFaces_:" << nInternalFaces_ << endl;
    Pout<< "triangles:" << triangles.size() << endl;


    geometricSurfacePatchList surfPatches(patches.size());

    forAll(patches, patchI)
    {
        surfPatches[patchI] =
            geometricSurfacePatch
            (
                patches[patchI].physicalType(),
                patches[patchI].name(),
                patchI
            );
    }

    // Create globally numbered tri surface
    if (faceCentreDecomposition)
    {
        // Use newPoints (mesh points + face centres)
        triSurface globalSurf(triangles, surfPatches, newPoints);

        // Create locally numbered tri surface
        triSurface::operator=
        (
            triSurface
            (
                globalSurf.localFaces(),
                surfPatches,
                globalSurf.localPoints()
            )
        );
    }
    else
    {
        // Use mesh points
        triSurface globalSurf(triangles, surfPatches, mesh.points());

        // Create locally numbered tri surface
        triSurface::operator=
        (
            triSurface
            (
                globalSurf.localFaces(),
                surfPatches,
                globalSurf.localPoints()
            )
        );
    }
}


// ************************************************************************* //
