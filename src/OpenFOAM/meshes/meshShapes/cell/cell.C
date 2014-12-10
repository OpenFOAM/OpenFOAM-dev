/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "cell.H"
#include "pyramidPointFaceRef.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::cell::typeName = "cell";

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::cell::labels(const faceUList& f) const
{
    // return the unordered list of vertex labels supporting the cell

    // count the maximum size of all vertices
    label maxVert = 0;

    const labelList& faces = *this;

    forAll(faces, faceI)
    {
        maxVert += f[faces[faceI]].size();
    }

    // set the fill-in list
    labelList p(maxVert);

    // in the first face there is no duplicates
    const labelList& first = f[faces[0]];

    forAll(first, pointI)
    {
        p[pointI] = first[pointI];
    }

    // re-use maxVert to count the real vertices
    maxVert = first.size();

    // go through the rest of the faces. For each vertex, check if the point is
    // already inserted (up to maxVert, which now carries the number of real
    // points. If not, add it at the end of the list.
    for (label faceI = 1; faceI < faces.size(); faceI++)
    {
        const labelList& curFace = f[faces[faceI]];

        forAll(curFace, pointI)
        {
            const label curPoint = curFace[pointI];

            bool found = false;

            for (register label checkI = 0; checkI < maxVert; checkI++)
            {
                if (curPoint == p[checkI])
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                p[maxVert] = curPoint;

                maxVert++;
            }
        }
    }

    // reset the size of the list
    p.setSize(maxVert);

    return p;
}


Foam::pointField Foam::cell::points
(
    const faceUList& f,
    const pointField& meshPoints
) const
{
    labelList pointLabels = labels(f);

    pointField p(pointLabels.size());

    forAll(p, i)
    {
        p[i] = meshPoints[pointLabels[i]];
    }

    return p;
}


Foam::edgeList Foam::cell::edges(const faceUList& f) const
{
    // return the lisf of cell edges

    const labelList& curFaces = *this;

    // create a list of edges
    label maxNoEdges = 0;

    forAll(curFaces, faceI)
    {
        maxNoEdges += f[curFaces[faceI]].nEdges();
    }

    edgeList allEdges(maxNoEdges);
    label nEdges = 0;

    forAll(curFaces, faceI)
    {
        const edgeList curFaceEdges = f[curFaces[faceI]].edges();

        forAll(curFaceEdges, faceEdgeI)
        {
            const edge& curEdge = curFaceEdges[faceEdgeI];

            bool edgeFound = false;

            for (label addedEdgeI = 0; addedEdgeI < nEdges; addedEdgeI++)
            {
                if (allEdges[addedEdgeI] == curEdge)
                {
                    edgeFound = true;

                    break;
                }
            }

            if (!edgeFound)
            {
                // Add the new edge onto the list
                allEdges[nEdges] = curEdge;
                nEdges++;
            }
        }
    }

    allEdges.setSize(nEdges);

    return allEdges;
}


Foam::point Foam::cell::centre
(
    const pointField& p,
    const faceUList& f
) const
{
    // When one wants to access the cell centre and magnitude, the
    // functionality on the mesh level should be used in preference to the
    // functions provided here. They do not rely to the functionality
    // implemented here, provide additional checking and are more efficient.
    // The cell::centre and cell::mag functions may be removed in the future.

    // WARNING!
    // In the old version of the code, it was possible to establish whether any
    // of the pyramids had a negative volume, caused either by the poor
    // estimate of the cell centre or by the fact that the cell was defined
    // inside out. In the new description of the cell, this can only be
    // established with the reference to the owner-neighbour face-cell
    // relationship, which is not available on this level. Thus, all the
    // pyramids are believed to be positive with no checking.

    // first calculate the aproximate cell centre as the average of all
    // face centres
    vector cEst = vector::zero;
    scalar sumArea = 0;

    const labelList& faces = *this;

    forAll(faces, faceI)
    {
        scalar a = f[faces[faceI]].mag(p);
        cEst += f[faces[faceI]].centre(p)*a;
        sumArea += a;
    }

    cEst /= sumArea + VSMALL;

    // Calculate the centre by breaking the cell into pyramids and
    // volume-weighted averaging their centres
    vector sumVc = vector::zero;

    scalar sumV = 0;

    forAll(faces, faceI)
    {
        // calculate pyramid volume. If it is greater than zero, OK.
        // If not, the pyramid is inside-out. Create a face with the opposite
        // order and recalculate pyramid centre!
        scalar pyrVol = pyramidPointFaceRef(f[faces[faceI]], cEst).mag(p);
        vector pyrCentre = pyramidPointFaceRef(f[faces[faceI]], cEst).centre(p);

        // if pyramid inside-out because face points inwards invert
        // N.B. pyramid remains unchanged
        if (pyrVol < 0)
        {
            pyrVol = -pyrVol;
        }

        sumVc += pyrVol*pyrCentre;
        sumV += pyrVol;
    }

    return sumVc/(sumV + VSMALL);
}


Foam::scalar Foam::cell::mag
(
    const pointField& p,
    const faceUList& f
) const
{
    // When one wants to access the cell centre and magnitude, the
    // functionality on the mesh level should be used in preference to the
    // functions provided here. They do not rely to the functionality
    // implemented here, provide additional checking and are more efficient.
    // The cell::centre and cell::mag functions may be removed in the future.

    // WARNING! See cell::centre

    // first calculate the aproximate cell centre as the average of all
    // face centres
    vector cEst = vector::zero;
    scalar nCellFaces = 0;

    const labelList& faces = *this;

    forAll(faces, faceI)
    {
        cEst += f[faces[faceI]].centre(p);
        nCellFaces += 1;
    }

    cEst /= nCellFaces;

    // Calculate the magnitude by summing the mags of the pyramids
    scalar v = 0;

    forAll(faces, faceI)
    {
        v += ::Foam::mag(pyramidPointFaceRef(f[faces[faceI]], cEst).mag(p));
    }

    return v;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool Foam::operator==(const cell& a, const cell& b)
{
    // Trivial reject: faces are different size
    if (a.size() != b.size())
    {
        return false;
    }

    List<bool> fnd(a.size(), false);

    forAll(b, bI)
    {
        label curLabel = b[bI];

        bool found = false;

        forAll(a, aI)
        {
            if (a[aI] == curLabel)
            {
                found = true;
                fnd[aI] = true;
                break;
            }
        }

        if (!found)
        {
            return false;
        }
    }

    // check if all faces on a were marked
    bool result = true;

    forAll(fnd, aI)
    {
        result = (result && fnd[aI]);
    }

    return result;
}


// ************************************************************************* //
