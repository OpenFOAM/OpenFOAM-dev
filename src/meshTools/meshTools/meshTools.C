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

\*---------------------------------------------------------------------------*/

#include "meshTools.H"
#include "polyMesh.H"
#include "hexMatcher.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Check if n is in same direction as normals of all faceLabels
bool Foam::meshTools::visNormal
(
    const vector& n,
    const vectorField& faceNormals,
    const labelList& faceLabels
)
{
    forAll(faceLabels, i)
    {
        if ((faceNormals[faceLabels[i]] & n) < SMALL)
        {
            // Found normal in different direction from n.
            return false;
        }
    }

    return true;
}


Foam::vectorField Foam::meshTools::calcBoxPointNormals(const primitivePatch& pp)
{
    vectorField octantNormal(8);
    octantNormal[mXmYmZ] = vector(-1, -1, -1);
    octantNormal[pXmYmZ] = vector(1, -1, -1);
    octantNormal[mXpYmZ] = vector(-1, 1, -1);
    octantNormal[pXpYmZ] = vector(1, 1, -1);

    octantNormal[mXmYpZ] = vector(-1, -1, 1);
    octantNormal[pXmYpZ] = vector(1, -1, 1);
    octantNormal[mXpYpZ] = vector(-1, 1, 1);
    octantNormal[pXpYpZ] = vector(1, 1, 1);

    octantNormal /= mag(octantNormal);


    vectorField pn(pp.nPoints());

    const vectorField& faceNormals = pp.faceNormals();
    const vectorField& pointNormals = pp.pointNormals();
    const labelListList& pointFaces = pp.pointFaces();

    forAll(pointFaces, pointI)
    {
        const labelList& pFaces = pointFaces[pointI];

        if (visNormal(pointNormals[pointI], faceNormals, pFaces))
        {
            pn[pointI] = pointNormals[pointI];
        }
        else
        {
            WarningIn
            (
                "Foam::meshTools::calcBoxPointNormals(const primitivePatch& pp)"
            )   << "Average point normal not visible for point:"
                << pp.meshPoints()[pointI] << endl;

            label visOctant =
                mXmYmZMask
              | pXmYmZMask
              | mXpYmZMask
              | pXpYmZMask
              | mXmYpZMask
              | pXmYpZMask
              | mXpYpZMask
              | pXpYpZMask;

            forAll(pFaces, i)
            {
                const vector& n = faceNormals[pFaces[i]];

                if (n.x() > SMALL)
                {
                    // All -x octants become invisible
                    visOctant &= ~mXmYmZMask;
                    visOctant &= ~mXmYpZMask;
                    visOctant &= ~mXpYmZMask;
                    visOctant &= ~mXpYpZMask;
                }
                else if (n.x() < -SMALL)
                {
                    // All +x octants become invisible
                    visOctant &= ~pXmYmZMask;
                    visOctant &= ~pXmYpZMask;
                    visOctant &= ~pXpYmZMask;
                    visOctant &= ~pXpYpZMask;
                }

                if (n.y() > SMALL)
                {
                    visOctant &= ~mXmYmZMask;
                    visOctant &= ~mXmYpZMask;
                    visOctant &= ~pXmYmZMask;
                    visOctant &= ~pXmYpZMask;
                }
                else if (n.y() < -SMALL)
                {
                    visOctant &= ~mXpYmZMask;
                    visOctant &= ~mXpYpZMask;
                    visOctant &= ~pXpYmZMask;
                    visOctant &= ~pXpYpZMask;
                }

                if (n.z() > SMALL)
                {
                    visOctant &= ~mXmYmZMask;
                    visOctant &= ~mXpYmZMask;
                    visOctant &= ~pXmYmZMask;
                    visOctant &= ~pXpYmZMask;
                }
                else if (n.z() < -SMALL)
                {
                    visOctant &= ~mXmYpZMask;
                    visOctant &= ~mXpYpZMask;
                    visOctant &= ~pXmYpZMask;
                    visOctant &= ~pXpYpZMask;
                }
            }

            label visI = -1;

            label mask = 1;

            for (label octant = 0; octant < 8; octant++)
            {
                if (visOctant & mask)
                {
                    // Take first visible octant found. Note: could use
                    // combination of visible octants here.
                    visI = octant;

                    break;
                }
                mask <<= 1;
            }

            if (visI != -1)
            {
                // Take a vector in this octant.
                pn[pointI] = octantNormal[visI];
            }
            else
            {
                pn[pointI] = vector::zero;

                WarningIn
                (
                    "Foam::meshTools::calcBoxPointNormals"
                    "(const primitivePatch& pp)"
                )   << "No visible octant for point:" << pp.meshPoints()[pointI]
                    << " cooord:" << pp.points()[pp.meshPoints()[pointI]] << nl
                    << "Normal set to " << pn[pointI] << endl;
            }
        }
    }

    return pn;
}


Foam::vector Foam::meshTools::normEdgeVec
(
    const primitiveMesh& mesh,
    const label edgeI
)
{
    vector eVec = mesh.edges()[edgeI].vec(mesh.points());

    eVec /= mag(eVec);

    return eVec;
}


void Foam::meshTools::writeOBJ
(
    Ostream& os,
    const point& pt
)
{
    os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
}


void Foam::meshTools::writeOBJ
(
    Ostream& os,
    const triad& t,
    const point& pt
)
{
    forAll(t, dirI)
    {
        writeOBJ(os, pt, pt + t[dirI]);
    }
}


void Foam::meshTools::writeOBJ
(
    Ostream& os,
    const point& p1,
    const point& p2,
    label& count
)
{
    os << "v" << ' ' << p1.x() << ' ' << p1.y() << ' ' << p1.z() << endl;
    os << "v" << ' ' << p2.x() << ' ' << p2.y() << ' ' << p2.z() << endl;

    os << "l" << " " << (count + 1) << " " << (count + 2) << endl;

    count += 2;
}


void Foam::meshTools::writeOBJ
(
    Ostream& os,
    const point& p1,
    const point& p2
)
{
    os << "v" << ' ' << p1.x() << ' ' << p1.y() << ' ' << p1.z() << endl;

    os << "vn"
        << ' ' << p2.x() - p1.x()
        << ' ' << p2.y() - p1.y()
        << ' ' << p2.z() - p1.z() << endl;
}


void Foam::meshTools::writeOBJ
(
    Ostream& os,
    const cellList& cells,
    const faceList& faces,
    const pointField& points,
    const labelList& cellLabels
)
{
    labelHashSet usedFaces(4*cellLabels.size());

    forAll(cellLabels, i)
    {
        const cell& cFaces = cells[cellLabels[i]];

        forAll(cFaces, j)
        {
            usedFaces.insert(cFaces[j]);
        }
    }

    writeOBJ(os, faces, points, usedFaces.toc());
}


bool Foam::meshTools::edgeOnCell
(
    const primitiveMesh& mesh,
    const label cellI,
    const label edgeI
)
{
    return findIndex(mesh.edgeCells(edgeI), cellI) != -1;
}


bool Foam::meshTools::edgeOnFace
(
    const primitiveMesh& mesh,
    const label faceI,
    const label edgeI
)
{
    return findIndex(mesh.faceEdges(faceI), edgeI) != -1;
}


// Return true if faceI part of cellI
bool Foam::meshTools::faceOnCell
(
    const primitiveMesh& mesh,
    const label cellI,
    const label faceI
)
{
    if (mesh.isInternalFace(faceI))
    {
        if
        (
            (mesh.faceOwner()[faceI] == cellI)
         || (mesh.faceNeighbour()[faceI] == cellI)
        )
        {
            return true;
        }
    }
    else
    {
        if (mesh.faceOwner()[faceI] == cellI)
        {
            return true;
        }
    }
    return false;
}


Foam::label Foam::meshTools::findEdge
(
    const edgeList& edges,
    const labelList& candidates,
    const label v0,
    const label v1
)
{
    forAll(candidates, i)
    {
        label edgeI = candidates[i];

        const edge& e = edges[edgeI];

        if ((e[0] == v0 && e[1] == v1) || (e[0] == v1 && e[1] == v0))
        {
            return edgeI;
        }
    }
    return -1;
}


Foam::label Foam::meshTools::findEdge
(
    const primitiveMesh& mesh,
    const label v0,
    const label v1
)
{
    const edgeList& edges = mesh.edges();

    const labelList& v0Edges = mesh.pointEdges()[v0];

    forAll(v0Edges, i)
    {
        label edgeI = v0Edges[i];

        const edge& e = edges[edgeI];

        if ((e.start() == v1) || (e.end() == v1))
        {
            return edgeI;
        }
    }
    return -1;
}


Foam::label Foam::meshTools::getSharedEdge
(
    const primitiveMesh& mesh,
    const label f0,
    const label f1
)
{
    const labelList& f0Edges = mesh.faceEdges()[f0];
    const labelList& f1Edges = mesh.faceEdges()[f1];

    forAll(f0Edges, f0EdgeI)
    {
        label edge0 = f0Edges[f0EdgeI];

        forAll(f1Edges, f1EdgeI)
        {
            label edge1 = f1Edges[f1EdgeI];

            if (edge0 == edge1)
            {
                return edge0;
            }
        }
    }
    FatalErrorIn
    (
        "meshTools::getSharedEdge(const primitiveMesh&, const label"
        ", const label)"
    )   << "Faces " << f0 << " and " << f1 << " do not share an edge"
        << abort(FatalError);

    return -1;

}


Foam::label Foam::meshTools::getSharedFace
(
    const primitiveMesh& mesh,
    const label cell0I,
    const label cell1I
)
{
    const cell& cFaces = mesh.cells()[cell0I];

    forAll(cFaces, cFaceI)
    {
        label faceI = cFaces[cFaceI];

        if
        (
            mesh.isInternalFace(faceI)
         && (
                mesh.faceOwner()[faceI] == cell1I
             || mesh.faceNeighbour()[faceI] == cell1I
            )
        )
        {
            return faceI;
        }
    }


    FatalErrorIn
    (
        "meshTools::getSharedFace(const primitiveMesh&, const label"
        ", const label)"
    )   << "No common face for"
        << "  cell0I:" << cell0I << "  faces:" << cFaces
        << "  cell1I:" << cell1I << "  faces:"
        << mesh.cells()[cell1I]
        << abort(FatalError);

    return -1;
}


// Get the two faces on cellI using edgeI.
void Foam::meshTools::getEdgeFaces
(
    const primitiveMesh& mesh,
    const label cellI,
    const label edgeI,
    label& face0,
    label& face1
)
{
    const labelList& eFaces = mesh.edgeFaces(edgeI);

    face0 = -1;
    face1 = -1;

    forAll(eFaces, eFaceI)
    {
        label faceI = eFaces[eFaceI];

        if (faceOnCell(mesh, cellI, faceI))
        {
            if (face0 == -1)
            {
                face0 = faceI;
            }
            else
            {
                face1 = faceI;

                return;
            }
        }
    }

    if ((face0 == -1) || (face1 == -1))
    {
        FatalErrorIn
        (
            "meshTools::getEdgeFaces(const primitiveMesh&, const label"
            ", const label, label&, label&"
        )   << "Can not find faces using edge " << mesh.edges()[edgeI]
            << " on cell " << cellI << abort(FatalError);
    }
}


// Return label of other edge connected to vertex
Foam::label Foam::meshTools::otherEdge
(
    const primitiveMesh& mesh,
    const labelList& edgeLabels,
    const label thisEdgeI,
    const label thisVertI
)
{
    forAll(edgeLabels, edgeLabelI)
    {
        label edgeI = edgeLabels[edgeLabelI];

        if (edgeI != thisEdgeI)
        {
            const edge& e = mesh.edges()[edgeI];

            if ((e.start() == thisVertI) || (e.end() == thisVertI))
            {
                return edgeI;
            }
        }
    }

    FatalErrorIn
    (
        "meshTools::otherEdge(const primitiveMesh&, const labelList&"
        ", const label, const label)"
    )   << "Can not find edge in "
        << UIndirectList<edge>(mesh.edges(), edgeLabels)()
        << " connected to edge "
        << thisEdgeI << " with vertices " << mesh.edges()[thisEdgeI]
        << " on side " << thisVertI << abort(FatalError);

    return -1;
}


// Return face on other side of edgeI
Foam::label Foam::meshTools::otherFace
(
    const primitiveMesh& mesh,
    const label cellI,
    const label faceI,
    const label edgeI
)
{
    label face0;
    label face1;

    getEdgeFaces(mesh, cellI, edgeI, face0, face1);

    if (face0 == faceI)
    {
        return face1;
    }
    else
    {
        return face0;
    }
}


// Return face on other side of edgeI
Foam::label Foam::meshTools::otherCell
(
    const primitiveMesh& mesh,
    const label otherCellI,
    const label faceI
)
{
    if (!mesh.isInternalFace(faceI))
    {
        FatalErrorIn
        (
            "meshTools::otherCell(const primitiveMesh&, const label"
            ", const label)"
        )   << "Face " << faceI << " is not internal"
            << abort(FatalError);
    }

    label newCellI = mesh.faceOwner()[faceI];

    if (newCellI == otherCellI)
    {
        newCellI = mesh.faceNeighbour()[faceI];
    }
    return newCellI;
}


// Returns label of edge nEdges away from startEdge (in the direction of
// startVertI)
Foam::label Foam::meshTools::walkFace
(
    const primitiveMesh& mesh,
    const label faceI,
    const label startEdgeI,
    const label startVertI,
    const label nEdges
)
{
    const labelList& fEdges = mesh.faceEdges(faceI);

    label edgeI = startEdgeI;

    label vertI = startVertI;

    for (label iter = 0; iter < nEdges; iter++)
    {
        edgeI = otherEdge(mesh, fEdges, edgeI, vertI);

        vertI = mesh.edges()[edgeI].otherVertex(vertI);
    }

    return edgeI;
}


void Foam::meshTools::constrainToMeshCentre
(
    const polyMesh& mesh,
    point& pt
)
{
    const Vector<label>& dirs = mesh.geometricD();

    const point& min = mesh.bounds().min();
    const point& max = mesh.bounds().max();

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        if (dirs[cmpt] == -1)
        {
            pt[cmpt] = 0.5*(min[cmpt] + max[cmpt]);
        }
    }
}


void Foam::meshTools::constrainToMeshCentre
(
    const polyMesh& mesh,
    pointField& pts
)
{
    const Vector<label>& dirs = mesh.geometricD();

    const point& min = mesh.bounds().min();
    const point& max = mesh.bounds().max();

    bool isConstrained = false;
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        if (dirs[cmpt] == -1)
        {
            isConstrained = true;
            break;
        }
    }

    if (isConstrained)
    {
        forAll(pts, i)
        {
            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
            {
                if (dirs[cmpt] == -1)
                {
                    pts[i][cmpt] = 0.5*(min[cmpt] + max[cmpt]);
                }
            }
        }
    }
}


//- Set the constrained components of directions/velocity to zero
void Foam::meshTools::constrainDirection
(
    const polyMesh& mesh,
    const Vector<label>& dirs,
    vector& d
)
{
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        if (dirs[cmpt] == -1)
        {
            d[cmpt] = 0.0;
        }
    }
}


void Foam::meshTools::constrainDirection
(
    const polyMesh& mesh,
    const Vector<label>& dirs,
    vectorField& d
)
{
    bool isConstrained = false;
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        if (dirs[cmpt] == -1)
        {
            isConstrained = true;
            break;
        }
    }

    if (isConstrained)
    {
        forAll(d, i)
        {
            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
            {
                if (dirs[cmpt] == -1)
                {
                    d[i][cmpt] = 0.0;
                }
            }
        }
    }
}


void Foam::meshTools::getParallelEdges
(
    const primitiveMesh& mesh,
    const label cellI,
    const label e0,
    label& e1,
    label& e2,
    label& e3
)
{
    // Go to any face using e0
    label faceI = meshTools::otherFace(mesh, cellI, -1, e0);

    // Opposite edge on face
    e1 = meshTools::walkFace(mesh, faceI, e0, mesh.edges()[e0].end(), 2);

    faceI = meshTools::otherFace(mesh, cellI, faceI, e1);

    e2 = meshTools::walkFace(mesh, faceI, e1, mesh.edges()[e1].end(), 2);

    faceI = meshTools::otherFace(mesh, cellI, faceI, e2);

    e3 = meshTools::walkFace(mesh, faceI, e2, mesh.edges()[e2].end(), 2);
}


Foam::vector Foam::meshTools::edgeToCutDir
(
    const primitiveMesh& mesh,
    const label cellI,
    const label startEdgeI
)
{
    if (!hexMatcher().isA(mesh, cellI))
    {
        FatalErrorIn
        (
            "Foam::meshTools::getCutDir(const label, const label)"
        )   << "Not a hex : cell:" << cellI << abort(FatalError);
    }


    vector avgVec(normEdgeVec(mesh, startEdgeI));

    label edgeI = startEdgeI;

    label faceI = -1;

    for (label i = 0; i < 3; i++)
    {
        // Step to next face, next edge
        faceI = meshTools::otherFace(mesh, cellI, faceI, edgeI);

        vector eVec(normEdgeVec(mesh, edgeI));

        if ((eVec & avgVec) > 0)
        {
            avgVec += eVec;
        }
        else
        {
            avgVec -= eVec;
        }

        label vertI = mesh.edges()[edgeI].end();

        edgeI = meshTools::walkFace(mesh, faceI, edgeI, vertI, 2);
    }

    avgVec /= mag(avgVec) + VSMALL;

    return avgVec;
}


// Find edges most aligned with cutDir
Foam::label Foam::meshTools::cutDirToEdge
(
    const primitiveMesh& mesh,
    const label cellI,
    const vector& cutDir
)
{
    if (!hexMatcher().isA(mesh, cellI))
    {
        FatalErrorIn
        (
            "Foam::meshTools::getCutDir(const label, const vector&)"
        )   << "Not a hex : cell:" << cellI << abort(FatalError);
    }

    const labelList& cEdges = mesh.cellEdges()[cellI];

    labelHashSet doneEdges(2*cEdges.size());

    scalar maxCos = -GREAT;
    label maxEdgeI = -1;

    for (label i = 0; i < 4; i++)
    {
        forAll(cEdges, cEdgeI)
        {
            label e0 = cEdges[cEdgeI];

            if (!doneEdges.found(e0))
            {
                vector avgDir(edgeToCutDir(mesh, cellI, e0));

                scalar cosAngle = mag(avgDir & cutDir);

                if (cosAngle > maxCos)
                {
                    maxCos = cosAngle;
                    maxEdgeI = e0;
                }

                // Mark off edges in cEdges.
                label e1, e2, e3;
                getParallelEdges(mesh, cellI, e0, e1, e2, e3);

                doneEdges.insert(e0);
                doneEdges.insert(e1);
                doneEdges.insert(e2);
                doneEdges.insert(e3);
            }
        }
    }

    forAll(cEdges, cEdgeI)
    {
        if (!doneEdges.found(cEdges[cEdgeI]))
        {
            FatalErrorIn
            (
                "meshTools::cutDirToEdge(const label, const vector&)"
            )   << "Cell:" << cellI << " edges:" << cEdges << endl
                << "Edge:" << cEdges[cEdgeI] << " not yet handled"
                << abort(FatalError);
        }
    }

    if (maxEdgeI == -1)
    {
        FatalErrorIn
        (
            "meshTools::cutDirToEdge(const label, const vector&)"
        )   << "Problem : did not find edge aligned with " << cutDir
            << " on cell " << cellI << abort(FatalError);
    }

    return maxEdgeI;
}


// ************************************************************************* //
