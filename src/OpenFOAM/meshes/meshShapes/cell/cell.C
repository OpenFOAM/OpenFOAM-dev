/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

    forAll(faces, facei)
    {
        maxVert += f[faces[facei]].size();
    }

    // set the fill-in list
    labelList p(maxVert);

    // in the first face there is no duplicates
    const labelList& first = f[faces[0]];

    forAll(first, pointi)
    {
        p[pointi] = first[pointi];
    }

    // re-use maxVert to count the real vertices
    maxVert = first.size();

    // go through the rest of the faces. For each vertex, check if the point is
    // already inserted (up to maxVert, which now carries the number of real
    // points. If not, add it at the end of the list.
    for (label facei = 1; facei < faces.size(); facei++)
    {
        const labelList& curFace = f[faces[facei]];

        forAll(curFace, pointi)
        {
            const label curPoint = curFace[pointi];

            bool found = false;

            for (label checkI = 0; checkI < maxVert; checkI++)
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

    forAll(curFaces, facei)
    {
        maxNoEdges += f[curFaces[facei]].nEdges();
    }

    edgeList allEdges(maxNoEdges);
    label nEdges = 0;

    forAll(curFaces, facei)
    {
        const edgeList curFaceEdges = f[curFaces[facei]].edges();

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

    // first calculate the approximate cell centre as the average of all
    // face centres
    vector cEst = Zero;
    scalar sumArea = 0;

    const labelList& faces = *this;

    forAll(faces, facei)
    {
        scalar a = f[faces[facei]].mag(p);
        cEst += f[faces[facei]].centre(p)*a;
        sumArea += a;
    }

    cEst /= sumArea + vSmall;

    // Calculate the centre by breaking the cell into pyramids and
    // volume-weighted averaging their centres
    vector sumVc = Zero;

    scalar sumV = 0;

    forAll(faces, facei)
    {
        // calculate pyramid volume. If it is greater than zero, OK.
        // If not, the pyramid is inside-out. Create a face with the opposite
        // order and recalculate pyramid centre!
        scalar pyrVol = pyramidPointFaceRef(f[faces[facei]], cEst).mag(p);
        vector pyrCentre = pyramidPointFaceRef(f[faces[facei]], cEst).centre(p);

        // if pyramid inside-out because face points inwards invert
        // N.B. pyramid remains unchanged
        if (pyrVol < 0)
        {
            pyrVol = -pyrVol;
        }

        sumVc += pyrVol*pyrCentre;
        sumV += pyrVol;
    }

    return sumVc/(sumV + vSmall);
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

    // first calculate the approximate cell centre as the average of all
    // face centres
    vector cEst = Zero;
    scalar nCellFaces = 0;

    const labelList& faces = *this;

    forAll(faces, facei)
    {
        cEst += f[faces[facei]].centre(p);
        nCellFaces += 1;
    }

    cEst /= nCellFaces;

    // Calculate the magnitude by summing the mags of the pyramids
    scalar v = 0;

    forAll(faces, facei)
    {
        v += ::Foam::mag(pyramidPointFaceRef(f[faces[facei]], cEst).mag(p));
    }

    return v;
}


Foam::boundBox Foam::cell::bb(const pointField& ps, const faceUList& fs) const
{
    boundBox result = boundBox::invertedBox;

    const cell& c = *this;

    forAll(c, cfi)
    {
        const face& f = fs[c[cfi]];

        forAll(f, fpi)
        {
            const point& p = ps[f[fpi]];

            result.min() = min(result.min(), p);
            result.max() = max(result.max(), p);
        }
    }

    return result;
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
