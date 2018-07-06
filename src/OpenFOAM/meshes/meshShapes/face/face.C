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

\*---------------------------------------------------------------------------*/

#include "face.H"
#include "triFace.H"
#include "triPointRef.H"
#include "mathematicalConstants.H"
#include "Swap.H"
#include "ConstCirculator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::face::typeName = "face";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::face::calcEdges(const pointField& points) const
{
    tmp<vectorField> tedges(new vectorField(size()));
    vectorField& edges = tedges.ref();

    forAll(*this, i)
    {
        label ni = fcIndex(i);

        point thisPt = points[operator[](i)];
        point nextPt = points[operator[](ni)];

        vector vec(nextPt - thisPt);
        vec /= Foam::mag(vec) + vSmall;

        edges[i] = vec;
    }

    return tedges;
}


Foam::scalar Foam::face::edgeCos
(
    const vectorField& edges,
    const label index
) const
{
    label leftEdgeI = left(index);
    label rightEdgeI = right(index);

    // Note negate on left edge to get correct left-pointing edge.
    return -(edges[leftEdgeI] & edges[rightEdgeI]);
}


Foam::label Foam::face::mostConcaveAngle
(
    const pointField& points,
    const vectorField& edges,
    scalar& maxAngle
) const
{
    vector a(area(points));

    label index = 0;
    maxAngle = -great;

    forAll(edges, i)
    {
        label leftEdgeI = left(i);
        label rightEdgeI = right(i);

        vector edgeNormal = edges[rightEdgeI] ^ edges[leftEdgeI];

        scalar edgeCos = edges[leftEdgeI] & edges[rightEdgeI];
        scalar edgeAngle = acos(max(-1.0, min(1.0, edgeCos)));

        scalar angle;

        if ((edgeNormal & a) > 0)
        {
            // Concave angle.
            angle = constant::mathematical::pi + edgeAngle;
        }
        else
        {
            // Convex angle. Note '-' to take into account that rightEdge
            // and leftEdge are head-to-tail connected.
            angle = constant::mathematical::pi - edgeAngle;
        }

        if (angle > maxAngle)
        {
            maxAngle = angle;
            index = i;
        }
    }

    return index;
}


Foam::label Foam::face::split
(
    const face::splitMode mode,
    const pointField& points,
    label& triI,
    label& quadI,
    faceList& triFaces,
    faceList& quadFaces
) const
{
    label oldIndices = (triI + quadI);

    if (size() <= 2)
    {
        FatalErrorInFunction
            << "Serious problem: asked to split a face with < 3 vertices"
            << abort(FatalError);
    }
    if (size() == 3)
    {
        // Triangle. Just copy.
        if (mode == COUNTTRIANGLE || mode == COUNTQUAD)
        {
            triI++;
        }
        else
        {
            triFaces[triI++] = *this;
        }
    }
    else if (size() == 4)
    {
        if (mode == COUNTTRIANGLE)
        {
            triI += 2;
        }
        if (mode == COUNTQUAD)
        {
            quadI++;
        }
        else if (mode == SPLITTRIANGLE)
        {
            //  Start at point with largest internal angle.
            const vectorField edges(calcEdges(points));

            scalar minAngle;
            label startIndex = mostConcaveAngle(points, edges, minAngle);

            label nextIndex = fcIndex(startIndex);
            label splitIndex = fcIndex(nextIndex);

            // Create triangles
            face triFace(3);
            triFace[0] = operator[](startIndex);
            triFace[1] = operator[](nextIndex);
            triFace[2] = operator[](splitIndex);

            triFaces[triI++] = triFace;

            triFace[0] = operator[](splitIndex);
            triFace[1] = operator[](fcIndex(splitIndex));
            triFace[2] = operator[](startIndex);

            triFaces[triI++] = triFace;
        }
        else
        {
            quadFaces[quadI++] = *this;
        }
    }
    else
    {
        // General case. Like quad: search for largest internal angle.

        const vectorField edges(calcEdges(points));

        scalar minAngle = 1;
        label startIndex = mostConcaveAngle(points, edges, minAngle);

        scalar bisectAngle = minAngle/2;
        vector rightEdge = edges[right(startIndex)];

        //
        // Look for opposite point which as close as possible bisects angle
        //

        // split candidate starts two points away.
        label index = fcIndex(fcIndex(startIndex));

        label minIndex = index;
        scalar minDiff = constant::mathematical::pi;

        for (label i = 0; i < size() - 3; i++)
        {
            vector splitEdge
            (
                points[operator[](index)]
              - points[operator[](startIndex)]
            );
            splitEdge /= Foam::mag(splitEdge) + vSmall;

            const scalar splitCos = splitEdge & rightEdge;
            const scalar splitAngle = acos(max(-1.0, min(1.0, splitCos)));
            const scalar angleDiff = fabs(splitAngle - bisectAngle);

            if (angleDiff < minDiff)
            {
                minDiff = angleDiff;
                minIndex = index;
            }

            // Go to next candidate
            index = fcIndex(index);
        }


        // Split into two subshapes.
        //     face1: startIndex to minIndex
        //     face2: minIndex to startIndex

        // Get sizes of the two subshapes
        label diff = 0;
        if (minIndex > startIndex)
        {
            diff = minIndex - startIndex;
        }
        else
        {
            // Folded around
            diff = minIndex + size() - startIndex;
        }

        label nPoints1 = diff + 1;
        label nPoints2 = size() - diff + 1;

        // Collect face1 points
        face face1(nPoints1);

        index = startIndex;
        for (label i = 0; i < nPoints1; i++)
        {
            face1[i] = operator[](index);
            index = fcIndex(index);
        }

        // Collect face2 points
        face face2(nPoints2);

        index = minIndex;
        for (label i = 0; i < nPoints2; i++)
        {
            face2[i] = operator[](index);
            index = fcIndex(index);
        }

        // Split faces
        face1.split(mode, points, triI, quadI, triFaces, quadFaces);
        face2.split(mode, points, triI, quadI, triFaces, quadFaces);
    }

    return (triI + quadI - oldIndices);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::face::face(const triFace& f)
:
    labelList(f)
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

int Foam::face::compare(const face& a, const face& b)
{
    // Basic rule: we assume that the sequence of labels in each list
    // will be circular in the same order (but not necessarily in the
    // same direction or from the same starting point).

    // Trivial reject: faces are different size
    label sizeA = a.size();
    label sizeB = b.size();

    if (sizeA != sizeB || sizeA == 0)
    {
        return 0;
    }
    else if (sizeA == 1)
    {
        if (a[0] == b[0])
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }

    ConstCirculator<face> aCirc(a);
    ConstCirculator<face> bCirc(b);

    // Rotate face b until its element matches the starting element of face a.
    do
    {
        if (aCirc() == bCirc())
        {
            // Set bCirc fulcrum to its iterator and increment the iterators
            bCirc.setFulcrumToIterator();
            ++aCirc;
            ++bCirc;

            break;
        }
    } while (bCirc.circulate(CirculatorBase::CLOCKWISE));

    // If the circulator has stopped then faces a and b do not share a matching
    // point. Doesn't work on matching, single element face.
    if (!bCirc.circulate())
    {
        return 0;
    }

    // Look forwards around the faces for a match
    do
    {
        if (aCirc() != bCirc())
        {
            break;
        }
    }
    while
    (
        aCirc.circulate(CirculatorBase::CLOCKWISE),
        bCirc.circulate(CirculatorBase::CLOCKWISE)
    );

    // If the circulator has stopped then faces a and b matched.
    if (!aCirc.circulate())
    {
        return 1;
    }
    else
    {
        // Reset the circulators back to their fulcrum
        aCirc.setIteratorToFulcrum();
        bCirc.setIteratorToFulcrum();
        ++aCirc;
        --bCirc;
    }

    // Look backwards around the faces for a match
    do
    {
        if (aCirc() != bCirc())
        {
            break;
        }
    }
    while
    (
        aCirc.circulate(CirculatorBase::CLOCKWISE),
        bCirc.circulate(CirculatorBase::ANTICLOCKWISE)
    );

    // If the circulator has stopped then faces a and b matched.
    if (!aCirc.circulate())
    {
        return -1;
    }

    return 0;
}


bool Foam::face::sameVertices(const face& a, const face& b)
{
    label sizeA = a.size();
    label sizeB = b.size();

    // Trivial reject: faces are different size
    if (sizeA != sizeB)
    {
        return false;
    }
    // Check faces with a single vertex
    else if (sizeA == 1)
    {
        if (a[0] == b[0])
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    forAll(a, i)
    {
        // Count occurrences of a[i] in a
        label aOcc = 0;
        forAll(a, j)
        {
            if (a[i] == a[j]) aOcc++;
        }

        // Count occurrences of a[i] in b
        label bOcc = 0;
        forAll(b, j)
        {
            if (a[i] == b[j]) bOcc++;
        }

        // Check if occurrences of a[i] in a and b are the same
        if (aOcc != bOcc) return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::face::collapse()
{
    if (size() > 1)
    {
        label ci = 0;
        for (label i=1; i<size(); i++)
        {
            if (operator[](i) != operator[](ci))
            {
                operator[](++ci) = operator[](i);
            }
        }

        if (operator[](ci) != operator[](0))
        {
            ci++;
        }

        setSize(ci);
    }

    return size();
}


void Foam::face::flip()
{
    const label n = size();

    if (n > 2)
    {
        for (label i=1; i < (n+1)/2; ++i)
        {
            Swap(operator[](i), operator[](n-i));
        }
    }
}


Foam::point Foam::face::centre(const pointField& points) const
{
    // Calculate the centre by breaking the face into triangles and
    // area-weighted averaging their centres

    const label nPoints = size();

    // If the face is a triangle, do a direct calculation
    if (nPoints == 3)
    {
        return
            (1.0/3.0)
           *(
               points[operator[](0)]
             + points[operator[](1)]
             + points[operator[](2)]
            );
    }


    point centrePoint = Zero;
    for (label pI=0; pI<nPoints; ++pI)
    {
        centrePoint += points[operator[](pI)];
    }
    centrePoint /= nPoints;

    scalar sumA = 0;
    vector sumAc = Zero;

    for (label pI=0; pI<nPoints; ++pI)
    {
        const point& nextPoint = points[operator[]((pI + 1) % nPoints)];

        // Calculate 3*triangle centre
        const vector ttc
        (
            points[operator[](pI)]
          + nextPoint
          + centrePoint
        );

        // Calculate 2*triangle area
        const scalar ta = Foam::mag
        (
            (points[operator[](pI)] - centrePoint)
          ^ (nextPoint - centrePoint)
        );

        sumA += ta;
        sumAc += ta*ttc;
    }

    if (sumA > vSmall)
    {
        return sumAc/(3.0*sumA);
    }
    else
    {
        return centrePoint;
    }
}


Foam::vector Foam::face::area(const pointField& p) const
{
    const label nPoints = size();

    // Calculate the area by summing the face triangle areas.
    // Changed to deal with small concavity by using a central decomposition
    //

    // If the face is a triangle, do a direct calculation to avoid round-off
    // error-related problems
    //
    if (nPoints == 3)
    {
        return triPointRef
        (
            p[operator[](0)],
            p[operator[](1)],
            p[operator[](2)]
        ).area();
    }

    label pI;

    point centrePoint = Zero;
    for (pI = 0; pI < nPoints; ++pI)
    {
        centrePoint += p[operator[](pI)];
    }
    centrePoint /= nPoints;

    vector a = Zero;

    point nextPoint = centrePoint;

    for (pI = 0; pI < nPoints; ++pI)
    {
        if (pI < nPoints - 1)
        {
            nextPoint = p[operator[](pI + 1)];
        }
        else
        {
            nextPoint = p[operator[](0)];
        }

        // Note: for best accuracy, centre point always comes last
        //
        a += triPointRef
        (
            p[operator[](pI)],
            nextPoint,
            centrePoint
        ).area();
    }

    return a;
}


Foam::vector Foam::face::normal(const pointField& points) const
{
    const vector a = area(points);
    const scalar maga = Foam::mag(a);
    return maga > 0 ? a/maga : Zero;
}


Foam::face Foam::face::reverseFace() const
{
    // Reverse the label list and return
    // The starting points of the original and reverse face are identical.

    const labelList& f = *this;
    labelList newList(size());

    newList[0] = f[0];

    for (label pointi = 1; pointi < newList.size(); pointi++)
    {
        newList[pointi] = f[size() - pointi];
    }

    return face(xferMove(newList));
}


Foam::label Foam::face::which(const label globalIndex) const
{
    const labelList& f = *this;

    forAll(f, localIdx)
    {
        if (f[localIdx] == globalIndex)
        {
            return localIdx;
        }
    }

    return -1;
}


Foam::scalar Foam::face::sweptVol
(
    const pointField& oldPoints,
    const pointField& newPoints
) const
{
    // This Optimization causes a small discrepancy between the swept-volume of
    // opposite faces of complex cells with triangular faces opposing polygons.
    // It could be used without problem for tetrahedral cells
    // if (size() == 3)
    // {
    //     return
    //     (
    //         triPointRef
    //         (
    //             oldPoints[operator[](0)],
    //             oldPoints[operator[](1)],
    //             oldPoints[operator[](2)]
    //         ).sweptVol
    //         (
    //             triPointRef
    //             (
    //                 newPoints[operator[](0)],
    //                 newPoints[operator[](1)],
    //                 newPoints[operator[](2)]
    //             )
    //         )
    //     );
    // }

    scalar sv = 0;

    // Calculate the swept volume by breaking the face into triangles and
    // summing their swept volumes.
    // Changed to deal with small concavity by using a central decomposition

    point centreOldPoint = centre(oldPoints);
    point centreNewPoint = centre(newPoints);

    label nPoints = size();

    for (label pi=0; pi<nPoints-1; ++pi)
    {
        // Note: for best accuracy, centre point always comes last
        sv += triPointRef
        (
            centreOldPoint,
            oldPoints[operator[](pi)],
            oldPoints[operator[](pi + 1)]
        ).sweptVol
        (
            triPointRef
            (
                centreNewPoint,
                newPoints[operator[](pi)],
                newPoints[operator[](pi + 1)]
            )
        );
    }

    sv += triPointRef
    (
        centreOldPoint,
        oldPoints[operator[](nPoints-1)],
        oldPoints[operator[](0)]
    ).sweptVol
    (
        triPointRef
        (
            centreNewPoint,
            newPoints[operator[](nPoints-1)],
            newPoints[operator[](0)]
        )
    );

    return sv;
}


Foam::tensor Foam::face::inertia
(
    const pointField& p,
    const point& refPt,
    scalar density
) const
{
    // If the face is a triangle, do a direct calculation
    if (size() == 3)
    {
        return triPointRef
        (
            p[operator[](0)],
            p[operator[](1)],
            p[operator[](2)]
        ).inertia(refPt, density);
    }

    const point ctr = centre(p);

    tensor J = Zero;

    forAll(*this, i)
    {
        J += triPointRef
        (
            p[operator[](i)],
            p[operator[](fcIndex(i))],
            ctr
        ).inertia(refPt, density);
    }

    return J;
}


Foam::edgeList Foam::face::edges() const
{
    const labelList& points = *this;

    edgeList e(points.size());

    for (label pointi = 0; pointi < points.size() - 1; ++pointi)
    {
        e[pointi] = edge(points[pointi], points[pointi + 1]);
    }

    // Add last edge
    e.last() = edge(points.last(), points[0]);

    return e;
}


int Foam::face::edgeDirection(const edge& e) const
{
    forAll(*this, i)
    {
        if (operator[](i) == e.start())
        {
            if (operator[](rcIndex(i)) == e.end())
            {
                // Reverse direction
                return -1;
            }
            else if (operator[](fcIndex(i)) == e.end())
            {
                // Forward direction
                return 1;
            }

            // No match
            return 0;
        }
        else if (operator[](i) == e.end())
        {
            if (operator[](rcIndex(i)) == e.start())
            {
                // Forward direction
                return 1;
            }
            else if (operator[](fcIndex(i)) == e.start())
            {
                // Reverse direction
                return -1;
            }

            // No match
            return 0;
        }
    }

    // Not found
    return 0;
}


Foam::label Foam::face::nTriangles(const pointField&) const
{
    return nTriangles();
}


Foam::label Foam::face::triangles
(
    const pointField& points,
    label& triI,
    faceList& triFaces
) const
{
    label quadI = 0;
    faceList quadFaces;

    return split(SPLITTRIANGLE, points, triI, quadI, triFaces, quadFaces);
}


Foam::label Foam::face::nTrianglesQuads
(
    const pointField& points,
    label& triI,
    label& quadI
) const
{
    faceList triFaces;
    faceList quadFaces;

    return split(COUNTQUAD, points, triI, quadI, triFaces, quadFaces);
}


Foam::label Foam::face::trianglesQuads
(
    const pointField& points,
    label& triI,
    label& quadI,
    faceList& triFaces,
    faceList& quadFaces
) const
{
    return split(SPLITQUAD, points, triI, quadI, triFaces, quadFaces);
}


Foam::label Foam::longestEdge(const face& f, const pointField& pts)
{
    const edgeList& eds = f.edges();

    label longestEdgeI = -1;
    scalar longestEdgeLength = -small;

    forAll(eds, edI)
    {
        scalar edgeLength = eds[edI].mag(pts);

        if (edgeLength > longestEdgeLength)
        {
            longestEdgeI = edI;
            longestEdgeLength = edgeLength;
        }
    }

    return longestEdgeI;
}


// ************************************************************************* //
