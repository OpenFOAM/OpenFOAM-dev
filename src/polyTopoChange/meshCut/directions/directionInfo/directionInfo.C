/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "directionInfo.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Find edge among edgeLabels that uses v0 and v1
Foam::label Foam::directionInfo::findEdge
(
    const primitiveMesh& mesh,
    const labelList& edgeLabels,
    const label v1,
    const label v0
)
{
    forAll(edgeLabels, edgeLabelI)
    {
        label edgeI = edgeLabels[edgeLabelI];

        if (mesh.edges()[edgeI] == edge(v0, v1))
        {
            return edgeI;
        }
    }

    FatalErrorInFunction
        << "Cannot find an edge among " << edgeLabels << endl
        << "that uses vertices " << v0
        << " and " << v1
        << abort(FatalError);

    return -1;
}


Foam::label Foam::directionInfo::lowest
(
    const label size,
    const label a,
    const label b
)
{
    // Get next point
    label a1 = (a + 1) % size;

    if (a1 == b)
    {
        return a;
    }
    else
    {
        label b1 = (b + 1) % size;

        if (b1 != a)
        {
            FatalErrorInFunction
                << "Problem : a:" << a << " b:" << b << " size:" << size
                << abort(FatalError);
        }

        return b;
    }
}


// Have edge on hex cell. Find corresponding edge on face. Return -1 if none
// found.
Foam::label Foam::directionInfo::edgeToFaceIndex
(
    const primitiveMesh& mesh,
    const label celli,
    const label facei,
    const label edgeI
)
{
    if ((edgeI < 0) || (edgeI >= mesh.nEdges()))
    {
        FatalErrorInFunction
            << "Illegal edge label:" << edgeI
            << " when projecting cut edge from cell " << celli
            << " to face " << facei
            << abort(FatalError);
    }

    const edge& e = mesh.edges()[edgeI];

    const face& f = mesh.faces()[facei];

    // edgeI is either
    // - in facei. Convert into index in face.
    // - connected (but not in) to face. Return -1.
    // - in face opposite facei. Convert into index in face.

    label fpA = findIndex(f, e.start());
    label fpB = findIndex(f, e.end());

    if (fpA != -1)
    {
        if (fpB != -1)
        {
            return lowest(f.size(), fpA, fpB);
        }
        else
        {
            // e.start() in face, e.end() not
            return -1;
        }
    }
    else
    {
        if (fpB != -1)
        {
            // e.end() in face, e.start() not
            return -1;
        }
        else
        {
            // Both not in face.
            // e is on opposite face. Determine corresponding edge on this face:
            // - determine two faces using edge (one is the opposite face,
            //   one is 'side' face
            // - walk on both these faces to opposite edge
            // - check if this opposite edge is on facei

            label f0I, f1I;

            meshTools::getEdgeFaces(mesh, celli, edgeI, f0I, f1I);

            // Walk to opposite edge on face f0
            label edge0I =
                meshTools::walkFace(mesh, f0I, edgeI, e.start(), 2);

            // Check if edge on facei.

            const edge& e0 = mesh.edges()[edge0I];

            fpA = findIndex(f, e0.start());
            fpB = findIndex(f, e0.end());

            if ((fpA != -1) && (fpB != -1))
            {
                return lowest(f.size(), fpA, fpB);
            }

            // Face0 is doesn't have an edge on facei (so must be the opposite
            // face) so try face1.

            // Walk to opposite edge on face f1
            label edge1I =
                meshTools::walkFace(mesh, f1I, edgeI, e.start(), 2);

            // Check if edge on facei.
            const edge& e1 = mesh.edges()[edge1I];

            fpA = findIndex(f, e1.start());
            fpB = findIndex(f, e1.end());

            if ((fpA != -1) && (fpB != -1))
            {
                return lowest(f.size(), fpA, fpB);
            }

            FatalErrorInFunction
                << "Found connected faces " << mesh.faces()[f0I] << " and "
                << mesh.faces()[f1I] << " sharing edge " << edgeI << endl
                << "But none seems to be connected to face " << facei
                << " vertices:" << f
                << abort(FatalError);

            return -1;
        }
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::directionInfo& wDist
)
{
    if (os.format() == IOstream::ASCII)
    {
        os << wDist.index_ << wDist.n_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&wDist.index_),
            sizeof(directionInfo)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const directionInfo&)");
    return os;

}


Foam::Istream& Foam::operator>>(Foam::Istream& is, Foam::directionInfo& wDist)
{
    if (is.format() == IOstream::ASCII)
    {
        is >> wDist.index_ >> wDist.n_;
    }
    else
    {
        is.read
        (
            reinterpret_cast<char*>(&wDist.index_),
            sizeof(directionInfo)
        );
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, directionInfo&)");
    return is;
}

// ************************************************************************* //
