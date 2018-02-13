/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "orientedSurface.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(orientedSurface, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Return true if edge is used in opposite order in faces
bool Foam::orientedSurface::consistentEdge
(
    const edge& e,
    const triSurface::FaceType& f0,
    const triSurface::FaceType& f1
)
{
    return (f0.edgeDirection(e) > 0) ^ (f1.edgeDirection(e) > 0);
}


Foam::labelList Foam::orientedSurface::faceToEdge
(
    const triSurface& s,
    const labelList& changedFaces
)
{
    labelList changedEdges(3*changedFaces.size());
    label changedI = 0;

    forAll(changedFaces, i)
    {
        const labelList& fEdges = s.faceEdges()[changedFaces[i]];

        forAll(fEdges, j)
        {
            changedEdges[changedI++] = fEdges[j];
        }
    }
    changedEdges.setSize(changedI);

    return changedEdges;
}


Foam::labelList Foam::orientedSurface::edgeToFace
(
    const triSurface& s,
    const labelList& changedEdges,
    labelList& flip
)
{
    labelList changedFaces(2*changedEdges.size());
    label changedI = 0;

    forAll(changedEdges, i)
    {
        label edgeI = changedEdges[i];

        const labelList& eFaces = s.edgeFaces()[edgeI];

        if (eFaces.size() < 2)
        {
            // Do nothing, faces was already visited.
        }
        else if (eFaces.size() == 2)
        {
            label face0 = eFaces[0];
            label face1 = eFaces[1];

            const triSurface::FaceType& f0 = s.localFaces()[face0];
            const triSurface::FaceType& f1 = s.localFaces()[face1];

            if (flip[face0] == UNVISITED)
            {
                if (flip[face1] == UNVISITED)
                {
                    FatalErrorInFunction
                        << abort(FatalError);
                }
                else
                {
                    // Face1 has a flip state, face0 hasn't
                    if (consistentEdge(s.edges()[edgeI], f0, f1))
                    {
                        // Take over flip status
                        flip[face0] = (flip[face1] == FLIP ? FLIP : NOFLIP);
                    }
                    else
                    {
                        // Invert
                        flip[face0] = (flip[face1] == FLIP ? NOFLIP : FLIP);
                    }
                    changedFaces[changedI++] = face0;
                }
            }
            else
            {
                if (flip[face1] == UNVISITED)
                {
                    // Face0 has a flip state, face1 hasn't
                    if (consistentEdge(s.edges()[edgeI], f0, f1))
                    {
                        flip[face1] = (flip[face0] == FLIP ? FLIP : NOFLIP);
                    }
                    else
                    {
                        flip[face1] = (flip[face0] == FLIP ? NOFLIP : FLIP);
                    }
                    changedFaces[changedI++] = face1;
                }
            }
        }
        else
        {
            // Multiply connected. Do what?
        }
    }
    changedFaces.setSize(changedI);

    return changedFaces;
}


void Foam::orientedSurface::walkSurface
(
    const triSurface& s,
    const label startFacei,
    labelList& flipState
)
{
    // List of faces that were changed in the last iteration.
    labelList changedFaces(1, startFacei);
    // List of edges that were changed in the last iteration.
    labelList changedEdges;

    while (true)
    {
        changedEdges = faceToEdge(s, changedFaces);

        if (changedEdges.empty())
        {
            break;
        }

        changedFaces = edgeToFace(s, changedEdges, flipState);

        if (changedFaces.empty())
        {
            break;
        }
    }
}


void Foam::orientedSurface::propagateOrientation
(
    const triSurface& s,
    const point& samplePoint,
    const bool orientOutside,
    const label nearestFacei,
    const point& nearestPt,
    labelList& flipState
)
{
    //
    // Determine orientation to normal on nearest face
    //
    triSurfaceTools::sideType side = triSurfaceTools::surfaceSide
    (
        s,
        samplePoint,
        nearestFacei
    );

    if (side == triSurfaceTools::UNKNOWN)
    {
        // Non-closed surface. Do what? For now behave as if no flipping
        // necessary
        flipState[nearestFacei] = NOFLIP;
    }
    else if ((side == triSurfaceTools::OUTSIDE) == orientOutside)
    {
        // outside & orientOutside or inside & !orientOutside
        // Normals on surface pointing correctly. No need to flip normals
        flipState[nearestFacei] = NOFLIP;
    }
    else
    {
        // Need to flip normals.
        flipState[nearestFacei] = FLIP;
    }

    if (debug)
    {
        vector n = triSurfaceTools::surfaceNormal(s, nearestFacei, nearestPt);

        Pout<< "orientedSurface::propagateOrientation : starting face"
            << " orientation:" << nl
            << "     for samplePoint:" << samplePoint << nl
            << "     starting from point:" << nearestPt << nl
            << "     on face:" << nearestFacei << nl
            << "     with normal:" << n << nl
            << "     decided side:" << label(side)
            << endl;
    }

    // Walk the surface from nearestFacei, changing the flipstate.
    walkSurface(s, nearestFacei, flipState);
}


// Find side for zoneI only by counting the number of intersections. Determines
// if face is oriented consistent with outwards pointing normals.
void Foam::orientedSurface::findZoneSide
(
    const triSurfaceSearch& surfSearches,
    const labelList& faceZone,
    const label zoneI,
    const point& outsidePoint,
    label& zoneFacei,
    bool& isOutside
)
{
    const triSurface& s = surfSearches.surface();

    zoneFacei = -1;
    isOutside = false;

    pointField start(1, outsidePoint);
    List<List<pointIndexHit>> hits(1, List<pointIndexHit>());

    forAll(faceZone, facei)
    {
        if (faceZone[facei] == zoneI)
        {
            const point& fc = s.faceCentres()[facei];
            const vector& n = s.faceNormals()[facei];

            const vector d = fc - outsidePoint;
            const scalar magD = mag(d);

            // Check if normal different enough to decide upon
            if (magD > small && (mag(n & d/magD) > 1e-6))
            {
                pointField end(1, fc + d);

                //Info<< "Zone " << zoneI << " : Shooting ray"
                //    << " from " << outsidePoint
                //    << " to " << end
                //    << " to pierce triangle " << facei
                //    << " with centre " << fc << endl;


                surfSearches.findLineAll(start, end, hits);

                label zoneIndex = -1;
                forAll(hits[0], i)
                {
                    if (hits[0][i].index() == facei)
                    {
                        zoneIndex = i;
                        break;
                    }
                }

                if (zoneIndex != -1)
                {
                    zoneFacei = facei;

                    if ((zoneIndex%2) == 0)
                    {
                        // Odd number of intersections. Check if normal points
                        // in direction of ray
                        isOutside = ((n & d) < 0);
                    }
                    else
                    {
                        isOutside = ((n & d) > 0);
                    }
                    break;
                }
            }
        }
    }
}


bool Foam::orientedSurface::flipSurface
(
    triSurface& s,
    const labelList& flipState
)
{
    bool hasFlipped = false;

    // Flip tris in s
    forAll(flipState, facei)
    {
        if (flipState[facei] == UNVISITED)
        {
            FatalErrorInFunction
                << "unvisited face " << facei
                << abort(FatalError);
        }
        else if (flipState[facei] == FLIP)
        {
            labelledTri& tri = s[facei];

            label tmp = tri[0];

            tri[0] = tri[1];
            tri[1] = tmp;

            hasFlipped = true;
        }
    }
    // Recalculate normals
    if (hasFlipped)
    {
        s.clearOut();
    }
    return hasFlipped;
}


bool Foam::orientedSurface::orientConsistent(triSurface& s)
{
    bool anyFlipped = false;

    // Do initial flipping to make triangles consistent. Otherwise if the
    // nearest is e.g. on an edge in between inconsistent triangles it might
    // make the wrong decision.
    if (s.size() > 0)
    {
        // Whether face has to be flipped.
        //      UNVISITED: unvisited
        //      NOFLIP: no need to flip
        //      FLIP: need to flip
        labelList flipState(s.size(), UNVISITED);

        label facei = 0;
        while (true)
        {
            label startFacei = -1;
            while (facei < s.size())
            {
                if (flipState[facei] == UNVISITED)
                {
                    startFacei = facei;
                    break;
                }
                facei++;
            }

            if (startFacei == -1)
            {
                break;
            }

            flipState[startFacei] = NOFLIP;
            walkSurface(s, startFacei, flipState);
        }

        anyFlipped = flipSurface(s, flipState);
    }
    return anyFlipped;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
Foam::orientedSurface::orientedSurface()
:
    triSurface()
{}


// Construct from surface and point which defines outside
Foam::orientedSurface::orientedSurface
(
    const triSurface& surf,
    const point& samplePoint,
    const bool orientOutside
)
:
    triSurface(surf)
{
    orient(*this, samplePoint, orientOutside);
}


// Construct from surface. Calculate outside point.
Foam::orientedSurface::orientedSurface
(
    const triSurface& surf,
    const bool orientOutside
)
:
    triSurface(surf)
{
    // BoundBox calculation without localPoints
    treeBoundBox bb(surf.points(), surf.meshPoints());

    point outsidePoint = bb.max() + bb.span();

    orient(*this, outsidePoint, orientOutside);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::orientedSurface::orient
(
    triSurface& s,
    const point& samplePoint,
    const bool orientOutside
)
{
    // Do initial flipping to make triangles consistent. Otherwise if the
    // nearest is e.g. on an edge in between inconsistent triangles it might
    // make the wrong decision.
    bool topoFlipped = orientConsistent(s);


    // Whether face has to be flipped.
    //      UNVISITED: unvisited
    //      NOFLIP: no need to flip
    //      FLIP: need to flip
    labelList flipState(s.size(), UNVISITED);


    while (true)
    {
        // Linear search for nearest unvisited point on surface.

        scalar minDist = great;
        point minPoint;
        label minFacei = -1;

        forAll(s, facei)
        {
            if (flipState[facei] == UNVISITED)
            {
                pointHit curHit =
                    s[facei].nearestPoint(samplePoint, s.points());

                if (curHit.distance() < minDist)
                {
                    minDist = curHit.distance();
                    minPoint = curHit.rawPoint();
                    minFacei = facei;
                }
            }
        }

        // Did we find anything?
        if (minFacei == -1)
        {
            break;
        }

        // From this nearest face see if needs to be flipped and then
        // go outwards.
        propagateOrientation
        (
            s,
            samplePoint,
            orientOutside,
            minFacei,
            minPoint,
            flipState
        );
    }

    // Now finally flip triangles according to flipState.
    bool geomFlipped = flipSurface(s, flipState);

    return topoFlipped || geomFlipped;
}


bool Foam::orientedSurface::orient
(
    triSurface& s,
    const triSurfaceSearch& querySurf,
    const point& samplePoint,
    const bool orientOutside
)
{
    // Do initial flipping to make triangles consistent. Otherwise if the
    // nearest is e.g. on an edge in between inconsistent triangles it might
    // make the wrong decision.
    bool topoFlipped = orientConsistent(s);

    // Determine disconnected parts of surface
    boolList borderEdge(s.nEdges(), false);
    forAll(s.edgeFaces(), edgeI)
    {
        if (s.edgeFaces()[edgeI].size() != 2)
        {
            borderEdge[edgeI] = true;
        }
    }
    labelList faceZone;
    label nZones = s.markZones(borderEdge, faceZone);

    // Check intersection with one face per zone.

    labelList flipState(s.size(), UNVISITED);
    for (label zoneI = 0; zoneI < nZones; zoneI++)
    {
        label zoneFacei = -1;
        bool isOutside;
        findZoneSide
        (
            querySurf,
            faceZone,
            zoneI,
            samplePoint,

            zoneFacei,
            isOutside
        );

        if (isOutside == orientOutside)
        {
            flipState[zoneFacei] = NOFLIP;
        }
        else
        {
            flipState[zoneFacei] = FLIP;
        }
        walkSurface(s, zoneFacei, flipState);
    }

    // Now finally flip triangles according to flipState.
    bool geomFlipped = flipSurface(s, flipState);

    return topoFlipped || geomFlipped;
}


// ************************************************************************* //
