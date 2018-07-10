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

#include "isoSurfaceCell.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "mergePoints.H"
#include "tetMatcher.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isoSurfaceCell, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::isoSurfaceCell::isoFraction
(
    const scalar s0,
    const scalar s1
) const
{
    scalar d = s1-s0;

    if (mag(d) > vSmall)
    {
        return (iso_-s0)/d;
    }
    else
    {
        return -1.0;
    }
}


bool Foam::isoSurfaceCell::isTriCut
(
    const triFace& tri,
    const scalarField& pointValues
) const
{
    bool aLower = (pointValues[tri[0]] < iso_);
    bool bLower = (pointValues[tri[1]] < iso_);
    bool cLower = (pointValues[tri[2]] < iso_);

    return !(aLower == bLower && aLower == cLower);
}


Foam::isoSurfaceCell::cellCutType Foam::isoSurfaceCell::calcCutType
(
    const PackedBoolList& isTet,
    const scalarField& cellValues,
    const scalarField& pointValues,
    const label celli
) const
{
    const cell& cFaces = mesh_.cells()[celli];

    if (isTet.get(celli) == 1)
    {
        forAll(cFaces, cFacei)
        {
            const face& f = mesh_.faces()[cFaces[cFacei]];

            for (label fp = 1; fp < f.size() - 1; fp++)
            {
                triFace tri(f[0], f[fp], f[f.fcIndex(fp)]);

                if (isTriCut(tri, pointValues))
                {
                    return CUT;
                }
            }
        }
        return NOTCUT;
    }
    else
    {
        bool cellLower = (cellValues[celli] < iso_);

        // First check if there is any cut in cell
        bool edgeCut = false;

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];
            const face& f = mesh_.faces()[facei];

            // Check pyramids cut
            forAll(f, fp)
            {
                if ((pointValues[f[fp]] < iso_) != cellLower)
                {
                    edgeCut = true;
                    break;
                }
            }

            if (edgeCut)
            {
                break;
            }

            const label fp0 = mesh_.tetBasePtIs()[facei];
            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); i++)
            {
                label nextFp = f.fcIndex(fp);

                if (isTriCut(triFace(f[fp0], f[fp], f[nextFp]), pointValues))
                {
                    edgeCut = true;
                    break;
                }

                fp = nextFp;
            }

            if (edgeCut)
            {
                break;
            }
        }

        if (edgeCut)
        {
            // Count actual cuts (expensive since addressing needed)
            // Note: not needed if you don't want to preserve maxima/minima
            // centred around cellcentre. In that case just always return CUT

            const labelList& cPoints = mesh_.cellPoints(celli);

            label nPyrCuts = 0;

            forAll(cPoints, i)
            {
                if ((pointValues[cPoints[i]] < iso_) != cellLower)
                {
                    nPyrCuts++;
                }
            }

            if (nPyrCuts == cPoints.size())
            {
                return SPHERE;
            }
            else
            {
                return CUT;
            }
        }
        else
        {
            return NOTCUT;
        }
    }
}


void Foam::isoSurfaceCell::calcCutTypes
(
    const PackedBoolList& isTet,
    const scalarField& cVals,
    const scalarField& pVals
)
{
    cellCutType_.setSize(mesh_.nCells());
    nCutCells_ = 0;
    forAll(mesh_.cells(), celli)
    {
        cellCutType_[celli] = calcCutType(isTet, cVals, pVals, celli);

        if (cellCutType_[celli] == CUT)
        {
            nCutCells_++;
        }
    }

    if (debug)
    {
        Pout<< "isoSurfaceCell : detected " << nCutCells_
            << " candidate cut cells." << endl;
    }
}


Foam::labelPair Foam::isoSurfaceCell::findCommonPoints
(
    const labelledTri& tri0,
    const labelledTri& tri1
)
{
    labelPair common(-1, -1);

    label fp0 = 0;
    label fp1 = findIndex(tri1, tri0[fp0]);

    if (fp1 == -1)
    {
        fp0 = 1;
        fp1 = findIndex(tri1, tri0[fp0]);
    }

    if (fp1 != -1)
    {
        // So tri0[fp0] is tri1[fp1]

        // Find next common point
        label fp0p1 = tri0.fcIndex(fp0);
        label fp1p1 = tri1.fcIndex(fp1);
        label fp1m1 = tri1.rcIndex(fp1);

        if (tri0[fp0p1] == tri1[fp1p1] || tri0[fp0p1] == tri1[fp1m1])
        {
            common[0] = tri0[fp0];
            common[1] = tri0[fp0p1];
        }
    }
    return common;
}


Foam::point Foam::isoSurfaceCell::calcCentre(const triSurface& s)
{
    vector sum = Zero;

    forAll(s, i)
    {
        sum += s[i].centre(s.points());
    }
    return sum/s.size();
}


Foam::pointIndexHit Foam::isoSurfaceCell::collapseSurface
(
    const label celli,
    pointField& localPoints,
    DynamicList<labelledTri, 64>& localTris
) const
{
    pointIndexHit info(false, Zero, localTris.size());

    if (localTris.size() == 1)
    {
        const labelledTri& tri = localTris[0];
        info.setPoint(tri.centre(localPoints));
        info.setHit();
    }
    else if (localTris.size() == 2)
    {
        // Check if the two triangles share an edge.
        const labelledTri& tri0 = localTris[0];
        const labelledTri& tri1 = localTris[1];

        labelPair shared = findCommonPoints(tri0, tri1);

        if (shared[0] != -1)
        {
            const vector n0 = tri0.normal(localPoints);
            const vector n1 = tri1.normal(localPoints);

            if ((n0 & n1) < 0)
            {
                // Too big an angle. Do not simplify.
            }
            else
            {
                info.setPoint
                (
                    0.5
                  * (
                        tri0.centre(localPoints)
                      + tri1.centre(localPoints)
                    )
                );
                info.setHit();
            }
        }
    }
    else if (localTris.size())
    {
        // Check if single region. Rare situation.
        triSurface surf
        (
            localTris,
            geometricSurfacePatchList(0),
            localPoints,
            true
        );
        localTris.clearStorage();

        labelList faceZone;
        label nZones = surf.markZones
        (
            boolList(surf.nEdges(), false),
            faceZone
        );

        if (nZones == 1)
        {
            // Check that all normals make a decent angle
            scalar minCos = great;
            const vector& n0 = surf.faceNormals()[0];
            for (label i = 1; i < surf.size(); i++)
            {
                scalar cosAngle = (n0 & surf.faceNormals()[i]);
                if (cosAngle < minCos)
                {
                    minCos = cosAngle;
                }
            }

            if (minCos > 0)
            {
                info.setPoint(calcCentre(surf));
                info.setHit();
            }
        }
    }

    return info;
}


void Foam::isoSurfaceCell::calcSnappedCc
(
    const PackedBoolList& isTet,
    const scalarField& cVals,
    const scalarField& pVals,

    DynamicList<point>& snappedPoints,
    labelList& snappedCc
) const
{
    const pointField& cc = mesh_.cellCentres();
    const pointField& pts = mesh_.points();

    snappedCc.setSize(mesh_.nCells());
    snappedCc = -1;

    // Work arrays
    DynamicList<point, 64> localPoints(64);
    DynamicList<labelledTri, 64> localTris(64);
    Map<label> pointToLocal(64);

    forAll(mesh_.cells(), celli)
    {
        if (cellCutType_[celli] == CUT && isTet.get(celli) == 0)
        {
            scalar cVal = cVals[celli];

            const cell& cFaces = mesh_.cells()[celli];

            localPoints.clear();
            localTris.clear();
            pointToLocal.clear();

            // Create points for all intersections close to cell centre
            // (i.e. from pyramid edges)

            forAll(cFaces, cFacei)
            {
                const face& f = mesh_.faces()[cFaces[cFacei]];

                forAll(f, fp)
                {
                    label pointi = f[fp];

                    scalar s = isoFraction(cVal, pVals[pointi]);

                    if (s >= 0.0 && s <= 0.5)
                    {
                        if (pointToLocal.insert(pointi, localPoints.size()))
                        {
                            localPoints.append((1.0-s)*cc[celli]+s*pts[pointi]);
                        }
                    }
                }
            }

            if (localPoints.size() == 1)
            {
                // No need for any analysis.
                snappedCc[celli] = snappedPoints.size();
                snappedPoints.append(localPoints[0]);

                // Pout<< "cell:" << celli
                //    << " at " << mesh_.cellCentres()[celli]
                //    << " collapsing " << localPoints
                //    << " intersections down to "
                //    << snappedPoints[snappedCc[celli]] << endl;
            }
            else if (localPoints.size() == 2)
            {
                //? No need for any analysis.???
                snappedCc[celli] = snappedPoints.size();
                snappedPoints.append(0.5*(localPoints[0]+localPoints[1]));

                // Pout<< "cell:" << celli
                //    << " at " << mesh_.cellCentres()[celli]
                //    << " collapsing " << localPoints
                //    << " intersections down to "
                //    << snappedPoints[snappedCc[celli]] << endl;
            }
            else if (localPoints.size())
            {
                // Need to analyse
                forAll(cFaces, cFacei)
                {
                    label facei = cFaces[cFacei];
                    const face& f = mesh_.faces()[facei];

                    // Do a tetrahedralisation. Each face to cc becomes pyr.
                    // Each pyr gets split into tets by diagonalisation
                    // of face.

                    const label fp0 = mesh_.tetBasePtIs()[facei];
                    label fp = f.fcIndex(fp0);
                    for (label i = 2; i < f.size(); i++)
                    {
                        label nextFp = f.fcIndex(fp);
                        triFace tri(f[fp0], f[fp], f[nextFp]);

                        // Get fractions for the three edges to cell centre
                        FixedList<scalar, 3> s(3);
                        s[0] = isoFraction(cVal, pVals[tri[0]]);
                        s[1] = isoFraction(cVal, pVals[tri[1]]);
                        s[2] = isoFraction(cVal, pVals[tri[2]]);

                        if
                        (
                            (s[0] >= 0.0 && s[0] <= 0.5)
                         && (s[1] >= 0.0 && s[1] <= 0.5)
                         && (s[2] >= 0.0 && s[2] <= 0.5)
                        )
                        {
                            if
                            (
                                (mesh_.faceOwner()[facei] == celli)
                             == (cVal >= pVals[tri[0]])
                            )
                            {
                                localTris.append
                                (
                                    labelledTri
                                    (
                                        pointToLocal[tri[1]],
                                        pointToLocal[tri[0]],
                                        pointToLocal[tri[2]],
                                        0
                                    )
                                );
                            }
                            else
                            {
                                localTris.append
                                (
                                    labelledTri
                                    (
                                        pointToLocal[tri[0]],
                                        pointToLocal[tri[1]],
                                        pointToLocal[tri[2]],
                                        0
                                    )
                                );
                            }
                        }

                        fp = nextFp;
                    }
                }

                pointField surfPoints;
                surfPoints.transfer(localPoints);
                pointIndexHit info = collapseSurface
                (
                    celli,
                    surfPoints,
                    localTris
                );

                if (info.hit())
                {
                    snappedCc[celli] = snappedPoints.size();
                    snappedPoints.append(info.hitPoint());

                    // Pout<< "cell:" << celli
                    //    << " at " << mesh_.cellCentres()[celli]
                    //    << " collapsing " << surfPoints
                    //    << " intersections down to "
                    //    << snappedPoints[snappedCc[celli]] << endl;
                }
            }
        }
    }
}


void Foam::isoSurfaceCell::genPointTris
(
    const scalarField& cellValues,
    const scalarField& pointValues,
    const label pointi,
    const label facei,
    const label celli,
    DynamicList<point, 64>& localTriPoints
) const
{
    const pointField& cc = mesh_.cellCentres();
    const pointField& pts = mesh_.points();
    const face& f = mesh_.faces()[facei];

    const label fp0 = mesh_.tetBasePtIs()[facei];
    label fp = f.fcIndex(fp0);
    for (label i = 2; i < f.size(); i++)
    {
        label nextFp = f.fcIndex(fp);
        triFace tri(f[fp0], f[fp], f[nextFp]);

        label index = findIndex(tri, pointi);

        if (index == -1)
        {
            continue;
        }

        // Tet between index..index-1, index..index+1, index..cc
        label b = tri[tri.fcIndex(index)];
        label c = tri[tri.rcIndex(index)];

        // Get fractions for the three edges emanating from point
        FixedList<scalar, 3> s(3);
        s[0] = isoFraction(pointValues[pointi], pointValues[b]);
        s[1] = isoFraction(pointValues[pointi], pointValues[c]);
        s[2] = isoFraction(pointValues[pointi], cellValues[celli]);

        if
        (
            (s[0] >= 0.0 && s[0] <= 0.5)
         && (s[1] >= 0.0 && s[1] <= 0.5)
         && (s[2] >= 0.0 && s[2] <= 0.5)
        )
        {
            point p0 = (1.0-s[0])*pts[pointi] + s[0]*pts[b];
            point p1 = (1.0-s[1])*pts[pointi] + s[1]*pts[c];
            point p2 = (1.0-s[2])*pts[pointi] + s[2]*cc[celli];

            if
            (
                (mesh_.faceOwner()[facei] == celli)
             == (pointValues[pointi] > cellValues[celli])
            )
            {
                localTriPoints.append(p0);
                localTriPoints.append(p1);
                localTriPoints.append(p2);
            }
            else
            {
                localTriPoints.append(p1);
                localTriPoints.append(p0);
                localTriPoints.append(p2);
            }
        }

        fp = nextFp;
    }
}


void Foam::isoSurfaceCell::genPointTris
(
    const scalarField& pointValues,
    const label pointi,
    const label facei,
    const label celli,
    DynamicList<point, 64>& localTriPoints
) const
{
    const pointField& pts = mesh_.points();
    const cell& cFaces = mesh_.cells()[celli];

    FixedList<label, 4> tet;

    // Make tet from this face to the 4th point (same as cellcentre in
    // non-tet cells)
    const face& f = mesh_.faces()[facei];

    // Find 4th point
    label ccPointi = -1;
    forAll(cFaces, cFacei)
    {
        const face& f1 = mesh_.faces()[cFaces[cFacei]];
        forAll(f1, fp)
        {
            label p1 = f1[fp];

            if (findIndex(f, p1) == -1)
            {
                ccPointi = p1;
                break;
            }
        }
        if (ccPointi != -1)
        {
            break;
        }
    }


    // Tet between index..index-1, index..index+1, index..cc
    label index = findIndex(f, pointi);
    label b = f[f.fcIndex(index)];
    label c = f[f.rcIndex(index)];

    // Pout<< " p0:" << pointi << " b:" << b << " c:" << c
    //<< " d:" << ccPointi << endl;

    // Get fractions for the three edges emanating from point
    FixedList<scalar, 3> s(3);
    s[0] = isoFraction(pointValues[pointi], pointValues[b]);
    s[1] = isoFraction(pointValues[pointi], pointValues[c]);
    s[2] = isoFraction(pointValues[pointi], pointValues[ccPointi]);

    if
    (
        (s[0] >= 0.0 && s[0] <= 0.5)
     && (s[1] >= 0.0 && s[1] <= 0.5)
     && (s[2] >= 0.0 && s[2] <= 0.5)
    )
    {
        point p0 = (1.0-s[0])*pts[pointi] + s[0]*pts[b];
        point p1 = (1.0-s[1])*pts[pointi] + s[1]*pts[c];
        point p2 = (1.0-s[2])*pts[pointi] + s[2]*pts[ccPointi];

        if (mesh_.faceOwner()[facei] != celli)
        {
            localTriPoints.append(p0);
            localTriPoints.append(p1);
            localTriPoints.append(p2);
        }
        else
        {
            localTriPoints.append(p1);
            localTriPoints.append(p0);
            localTriPoints.append(p2);
        }
    }
}


void Foam::isoSurfaceCell::calcSnappedPoint
(
    const PackedBoolList& isTet,
    const scalarField& cVals,
    const scalarField& pVals,

    DynamicList<point>& snappedPoints,
    labelList& snappedPoint
) const
{
    // Determine if point is on boundary. Points on boundaries are never
    // snapped. Coupled boundaries are handled explicitly so not marked here.
    PackedBoolList isBoundaryPoint(mesh_.nPoints());
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (!pp.coupled())
        {
            label facei = pp.start();
            forAll(pp, i)
            {
                const face& f = mesh_.faces()[facei++];

                forAll(f, fp)
                {
                    isBoundaryPoint.set(f[fp], 1);
                }
            }
        }
    }


    const point greatPoint(great, great, great);

    pointField collapsedPoint(mesh_.nPoints(), greatPoint);


    // Work arrays
    DynamicList<point, 64> localTriPoints(100);
    labelHashSet localPointCells(100);

    forAll(mesh_.pointFaces(), pointi)
    {
        if (isBoundaryPoint.get(pointi) == 1)
        {
            continue;
        }

        const labelList& pFaces = mesh_.pointFaces()[pointi];

        bool anyCut = false;

        forAll(pFaces, i)
        {
            label facei = pFaces[i];

            if
            (
                cellCutType_[mesh_.faceOwner()[facei]] == CUT
             || (
                    mesh_.isInternalFace(facei)
                 && cellCutType_[mesh_.faceNeighbour()[facei]] == CUT
                )
            )
            {
                anyCut = true;
                break;
            }
        }

        if (!anyCut)
        {
            continue;
        }


        // Do a pointCells walk (by using pointFaces)

        localPointCells.clear();
        localTriPoints.clear();

        forAll(pFaces, pFacei)
        {
            label facei = pFaces[pFacei];
            label own = mesh_.faceOwner()[facei];

            if (isTet.get(own) == 1)
            {
                // Since tets have no cell centre to include make sure
                // we only generate a triangle once per point.
                if (localPointCells.insert(own))
                {
                    genPointTris(pVals, pointi, facei, own, localTriPoints);
                }
            }
            else
            {
                genPointTris
                (
                    cVals,
                    pVals,
                    pointi,
                    facei,
                    own,
                    localTriPoints
                );
            }

            if (mesh_.isInternalFace(facei))
            {
                label nei = mesh_.faceNeighbour()[facei];

                if (isTet.get(nei) == 1)
                {
                    if (localPointCells.insert(nei))
                    {
                        genPointTris(pVals, pointi, facei, nei, localTriPoints);
                    }
                }
                else
                {
                    genPointTris
                    (
                        cVals,
                        pVals,
                        pointi,
                        facei,
                        nei,
                        localTriPoints
                    );
                }
            }
        }

        if (localTriPoints.size() == 3)
        {
            // Single triangle. No need for any analysis. Average points.
            pointField points;
            points.transfer(localTriPoints);
            collapsedPoint[pointi] = sum(points)/points.size();

            // Pout<< "    point:" << pointi
            //    << " replacing coord:" << mesh_.points()[pointi]
            //    << " by average:" << collapsedPoint[pointi] << endl;
        }
        else if (localTriPoints.size())
        {
            // Convert points into triSurface.

            // Merge points and compact out non-valid triangles
            labelList triMap;               // merged to unmerged triangle
            labelList triPointReverseMap;   // unmerged to merged point
            triSurface surf
            (
                stitchTriPoints
                (
                    false,                  // do not check for duplicate tris
                    localTriPoints,
                    triPointReverseMap,
                    triMap
                )
            );

            labelList faceZone;
            label nZones = surf.markZones
            (
                boolList(surf.nEdges(), false),
                faceZone
            );

            if (nZones == 1)
            {
                // Check that all normals make a decent angle
                scalar minCos = great;
                const vector& n0 = surf.faceNormals()[0];
                for (label i = 1; i < surf.size(); i++)
                {
                    const vector& n = surf.faceNormals()[i];
                    scalar cosAngle = (n0 & n);
                    if (cosAngle < minCos)
                    {
                        minCos = cosAngle;
                    }
                }
                if (minCos > 0)
                {
                    collapsedPoint[pointi] = calcCentre(surf);
                }
            }
        }
    }

    syncTools::syncPointPositions
    (
        mesh_,
        collapsedPoint,
        minMagSqrEqOp<point>(),
        greatPoint
    );

    snappedPoint.setSize(mesh_.nPoints());
    snappedPoint = -1;

    forAll(collapsedPoint, pointi)
    {
        // Cannot do == comparison since might be transformed so have
        // truncation errors.
        if (magSqr(collapsedPoint[pointi]) < 0.5*magSqr(greatPoint))
        {
            snappedPoint[pointi] = snappedPoints.size();
            snappedPoints.append(collapsedPoint[pointi]);
        }
    }
}


Foam::triSurface Foam::isoSurfaceCell::stitchTriPoints
(
    const bool checkDuplicates,
    const List<point>& triPoints,
    labelList& triPointReverseMap,  // unmerged to merged point
    labelList& triMap               // merged to unmerged triangle
) const
{
    label nTris = triPoints.size()/3;

    if ((triPoints.size() % 3) != 0)
    {
        FatalErrorInFunction
            << "Problem: number of points " << triPoints.size()
            << " not a multiple of 3." << abort(FatalError);
    }

    pointField newPoints;
    mergePoints
    (
        triPoints,
        mergeDistance_,
        false,
        triPointReverseMap,
        newPoints
    );

    // Check that enough merged.
    if (debug)
    {
        Pout<< "isoSurfaceCell : merged from " << triPoints.size()
            << " points down to " << newPoints.size() << endl;

        pointField newNewPoints;
        labelList oldToNew;
        bool hasMerged = mergePoints
        (
            newPoints,
            mergeDistance_,
            true,
            oldToNew,
            newNewPoints
        );

        if (hasMerged)
        {
            FatalErrorInFunction
                << "Merged points contain duplicates"
                << " when merging with distance " << mergeDistance_ << endl
                << "merged:" << newPoints.size() << " re-merged:"
                << newNewPoints.size()
                << abort(FatalError);
        }
    }


    List<labelledTri> tris;
    {
        DynamicList<labelledTri> dynTris(nTris);
        label rawPointi = 0;
        DynamicList<label> newToOldTri(nTris);

        for (label oldTriI = 0; oldTriI < nTris; oldTriI++)
        {
            labelledTri tri
            (
                triPointReverseMap[rawPointi],
                triPointReverseMap[rawPointi+1],
                triPointReverseMap[rawPointi+2],
                0
            );
            if ((tri[0] != tri[1]) && (tri[0] != tri[2]) && (tri[1] != tri[2]))
            {
                newToOldTri.append(oldTriI);
                dynTris.append(tri);
            }

            rawPointi += 3;
        }

        triMap.transfer(newToOldTri);
        tris.transfer(dynTris);
    }


    // Use face centres to determine 'flat hole' situation (see RMT paper).
    // Two unconnected triangles get connected because (some of) the edges
    // separating them get collapsed. Below only checks for duplicate triangles,
    // not non-manifold edge connectivity.
    if (checkDuplicates)
    {
        if (debug)
        {
            Pout<< "isoSurfaceCell : merged from " << nTris
                << " down to " << tris.size() << " triangles." << endl;
        }

        pointField centres(tris.size());
        forAll(tris, triI)
        {
            centres[triI] = tris[triI].centre(newPoints);
        }

        pointField mergedCentres;
        labelList oldToMerged;
        bool hasMerged = mergePoints
        (
            centres,
            mergeDistance_,
            false,
            oldToMerged,
            mergedCentres
        );

        if (debug)
        {
            Pout<< "isoSurfaceCell : detected "
                << centres.size()-mergedCentres.size()
                << " duplicate triangles." << endl;
        }

        if (hasMerged)
        {
            // Filter out duplicates.
            label newTriI = 0;
            DynamicList<label> newToOldTri(tris.size());
            labelList newToMaster(mergedCentres.size(), -1);
            forAll(tris, triI)
            {
                label mergedI = oldToMerged[triI];

                if (newToMaster[mergedI] == -1)
                {
                    newToMaster[mergedI] = triI;
                    newToOldTri.append(triMap[triI]);
                    tris[newTriI++] = tris[triI];
                }
            }

            triMap.transfer(newToOldTri);
            tris.setSize(newTriI);
        }
    }

    return triSurface(tris, geometricSurfacePatchList(0), newPoints, true);
}


bool Foam::isoSurfaceCell::validTri(const triSurface& surf, const label facei)
{
    // Simple check on indices ok.

    const labelledTri& f = surf[facei];

    forAll(f, fp)
    {
        if (f[fp] < 0 || f[fp] >= surf.points().size())
        {
            WarningInFunction
                << "triangle " << facei << " vertices " << f
                << " uses point indices outside point range 0.."
                << surf.points().size()-1 << endl;

            return false;
        }
    }

    if ((f[0] == f[1]) || (f[0] == f[2]) || (f[1] == f[2]))
    {
        WarningInFunction
            << "triangle " << facei
            << " uses non-unique vertices " << f
            << endl;
        return false;
    }

    // duplicate triangle check

    const labelList& fFaces = surf.faceFaces()[facei];

    // Check if faceNeighbours use same points as this face.
    // Note: discards normal information - sides of baffle are merged.
    forAll(fFaces, i)
    {
        label nbrFacei = fFaces[i];

        if (nbrFacei <= facei)
        {
            // lower numbered faces already checked
            continue;
        }

        const labelledTri& nbrF = surf[nbrFacei];

        if
        (
            ((f[0] == nbrF[0]) || (f[0] == nbrF[1]) || (f[0] == nbrF[2]))
         && ((f[1] == nbrF[0]) || (f[1] == nbrF[1]) || (f[1] == nbrF[2]))
         && ((f[2] == nbrF[0]) || (f[2] == nbrF[1]) || (f[2] == nbrF[2]))
        )
        {
            WarningInFunction
                << "triangle " << facei << " vertices " << f
                << " coords:" << f.points(surf.points())
                << " has the same vertices as triangle " << nbrFacei
                << " vertices " << nbrF
                << endl;

            return false;
        }
    }
    return true;
}


void Foam::isoSurfaceCell::calcAddressing
(
    const triSurface& surf,
    List<FixedList<label, 3>>& faceEdges,
    labelList& edgeFace0,
    labelList& edgeFace1,
    Map<labelList>& edgeFacesRest
) const
{
    const pointField& points = surf.points();

    pointField edgeCentres(3*surf.size());
    label edgeI = 0;
    forAll(surf, triI)
    {
        const labelledTri& tri = surf[triI];
        edgeCentres[edgeI++] = 0.5*(points[tri[0]]+points[tri[1]]);
        edgeCentres[edgeI++] = 0.5*(points[tri[1]]+points[tri[2]]);
        edgeCentres[edgeI++] = 0.5*(points[tri[2]]+points[tri[0]]);
    }

    pointField mergedCentres;
    labelList oldToMerged;
    bool hasMerged = mergePoints
    (
        edgeCentres,
        mergeDistance_,
        false,
        oldToMerged,
        mergedCentres
    );

    if (debug)
    {
        Pout<< "isoSurfaceCell : detected "
            << mergedCentres.size()
            << " edges on " << surf.size() << " triangles." << endl;
    }

    if (!hasMerged)
    {
        return;
    }


    // Determine faceEdges
    faceEdges.setSize(surf.size());
    edgeI = 0;
    forAll(surf, triI)
    {
        faceEdges[triI][0] = oldToMerged[edgeI++];
        faceEdges[triI][1] = oldToMerged[edgeI++];
        faceEdges[triI][2] = oldToMerged[edgeI++];
    }


    // Determine edgeFaces
    edgeFace0.setSize(mergedCentres.size());
    edgeFace0 = -1;
    edgeFace1.setSize(mergedCentres.size());
    edgeFace1 = -1;
    edgeFacesRest.clear();

    forAll(oldToMerged, oldEdgeI)
    {
        label triI = oldEdgeI / 3;
        label edgeI = oldToMerged[oldEdgeI];

        if (edgeFace0[edgeI] == -1)
        {
            edgeFace0[edgeI] = triI;
        }
        else if (edgeFace1[edgeI] == -1)
        {
            edgeFace1[edgeI] = triI;
        }
        else
        {
            // WarningInFunction
            //    << "Edge " << edgeI << " with centre " << mergedCentres[edgeI]
            //    << " used by more than two triangles: " << edgeFace0[edgeI]
            //    << ", "
            //    << edgeFace1[edgeI] << " and " << triI << endl;
            Map<labelList>::iterator iter = edgeFacesRest.find(edgeI);

            if (iter != edgeFacesRest.end())
            {
                labelList& eFaces = iter();
                label sz = eFaces.size();
                eFaces.setSize(sz+1);
                eFaces[sz] = triI;
            }
            else
            {
                edgeFacesRest.insert(edgeI, labelList(1, triI));
            }
        }
    }
}


bool Foam::isoSurfaceCell::danglingTriangle
(
    const FixedList<label, 3>& fEdges,
    const labelList& edgeFace1
)
{
    label nOpen = 0;
    forAll(fEdges, i)
    {
        if (edgeFace1[fEdges[i]] == -1)
        {
            nOpen++;
        }
    }

    if (nOpen == 1 || nOpen == 2 || nOpen == 3)
    {
        return true;
    }
    else
    {
        return false;
    }
}


Foam::label Foam::isoSurfaceCell::markDanglingTriangles
(
    const List<FixedList<label, 3>>& faceEdges,
    const labelList& edgeFace0,
    const labelList& edgeFace1,
    const Map<labelList>& edgeFacesRest,
    boolList& keepTriangles
)
{
    keepTriangles.setSize(faceEdges.size());
    keepTriangles = true;

    label nDangling = 0;

    // Remove any dangling triangles
    forAllConstIter(Map<labelList>,  edgeFacesRest, iter)
    {
        // These are all the non-manifold edges. Filter out all triangles
        // with only one connected edge (= this edge)

        label edgeI = iter.key();
        const labelList& otherEdgeFaces = iter();

        // Remove all dangling triangles
        if (danglingTriangle(faceEdges[edgeFace0[edgeI]], edgeFace1))
        {
            keepTriangles[edgeFace0[edgeI]] = false;
            nDangling++;
        }
        if (danglingTriangle(faceEdges[edgeFace1[edgeI]], edgeFace1))
        {
            keepTriangles[edgeFace1[edgeI]] = false;
            nDangling++;
        }
        forAll(otherEdgeFaces, i)
        {
            label triI = otherEdgeFaces[i];
            if (danglingTriangle(faceEdges[triI], edgeFace1))
            {
                keepTriangles[triI] = false;
                nDangling++;
            }
        }
    }
    return nDangling;
}


Foam::triSurface Foam::isoSurfaceCell::subsetMesh
(
    const triSurface& s,
    const labelList& newToOldFaces,
    labelList& oldToNewPoints,
    labelList& newToOldPoints
)
{
    const boolList include
    (
        createWithValues<boolList>
        (
            s.size(),
            false,
            newToOldFaces,
            true
        )
    );

    newToOldPoints.setSize(s.points().size());
    oldToNewPoints.setSize(s.points().size());
    oldToNewPoints = -1;
    {
        label pointi = 0;

        forAll(include, oldFacei)
        {
            if (include[oldFacei])
            {
                // Renumber labels for face
                const triSurface::FaceType& f = s[oldFacei];

                forAll(f, fp)
                {
                    label oldPointi = f[fp];

                    if (oldToNewPoints[oldPointi] == -1)
                    {
                        oldToNewPoints[oldPointi] = pointi;
                        newToOldPoints[pointi++] = oldPointi;
                    }
                }
            }
        }
        newToOldPoints.setSize(pointi);
    }

    // Extract points
    pointField newPoints(newToOldPoints.size());
    forAll(newToOldPoints, i)
    {
        newPoints[i] = s.points()[newToOldPoints[i]];
    }
    // Extract faces
    List<labelledTri> newTriangles(newToOldFaces.size());

    forAll(newToOldFaces, i)
    {
        // Get old vertex labels
        const labelledTri& tri = s[newToOldFaces[i]];

        newTriangles[i][0] = oldToNewPoints[tri[0]];
        newTriangles[i][1] = oldToNewPoints[tri[1]];
        newTriangles[i][2] = oldToNewPoints[tri[2]];
        newTriangles[i].region() = tri.region();
    }

    // Reuse storage.
    return triSurface(newTriangles, s.patches(), newPoints, true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoSurfaceCell::isoSurfaceCell
(
    const polyMesh& mesh,
    const scalarField& cVals,
    const scalarField& pVals,
    const scalar iso,
    const bool regularise,
    const scalar mergeTol
)
:
    mesh_(mesh),
    cVals_(cVals),
    pVals_(pVals),
    iso_(iso),
    mergeDistance_(mergeTol*mesh.bounds().mag())
{
    if (debug)
    {
        Pout<< "isoSurfaceCell : mergeTol:" << mergeTol
            << " mesh span:" << mesh.bounds().mag()
            << " mergeDistance:" << mergeDistance_ << endl;
    }

    // Determine if cell is tet
    PackedBoolList isTet(mesh_.nCells());
    {
        tetMatcher tet;

        forAll(isTet, celli)
        {
            if (tet.isA(mesh_, celli))
            {
                isTet.set(celli, 1);
            }
        }
    }


    // Determine if any cut through cell
    calcCutTypes(isTet, cVals, pVals);

    DynamicList<point> snappedPoints(nCutCells_);

    // Per cc -1 or a point inside snappedPoints.
    labelList snappedCc;
    if (regularise)
    {
        calcSnappedCc
        (
            isTet,
            cVals,
            pVals,
            snappedPoints,
            snappedCc
        );
    }
    else
    {
        snappedCc.setSize(mesh_.nCells());
        snappedCc = -1;
    }

    if (debug)
    {
        Pout<< "isoSurfaceCell : shifted " << snappedPoints.size()
            << " cell centres to intersection." << endl;
    }

    snappedPoints.shrink();
    label nCellSnaps = snappedPoints.size();

    // Per point -1 or a point inside snappedPoints.
    labelList snappedPoint;
    if (regularise)
    {
        calcSnappedPoint
        (
            isTet,
            cVals,
            pVals,
            snappedPoints,
            snappedPoint
        );
    }
    else
    {
        snappedPoint.setSize(mesh_.nPoints());
        snappedPoint = -1;
    }

    if (debug)
    {
        Pout<< "isoSurfaceCell : shifted " << snappedPoints.size()-nCellSnaps
            << " vertices to intersection." << endl;
    }


    {
        DynamicList<point> triPoints(nCutCells_);
        DynamicList<label> triMeshCells(nCutCells_);

        generateTriPoints
        (
            cVals,
            pVals,

            mesh_.cellCentres(),
            mesh_.points(),

            snappedPoints,
            snappedCc,
            snappedPoint,

            triPoints,
            triMeshCells
        );

        if (debug)
        {
            Pout<< "isoSurfaceCell : generated " << triMeshCells.size()
                << " unmerged triangles." << endl;
        }

        // Merge points and compact out non-valid triangles
        labelList triMap;
        triSurface::operator=
        (
            stitchTriPoints
            (
                regularise,         // check for duplicate tris
                triPoints,
                triPointMergeMap_,  // unmerged to merged point
                triMap              // merged to unmerged triangle
            )
        );

        if (debug)
        {
            Pout<< "isoSurfaceCell : generated " << triMap.size()
                << " merged triangles." << endl;
        }

        meshCells_.setSize(triMap.size());
        forAll(triMap, i)
        {
            meshCells_[i] = triMeshCells[triMap[i]];
        }
    }


    if (debug)
    {
        Pout<< "isoSurfaceCell : checking " << size()
            << " triangles for validity." << endl;

        forAll(*this, triI)
        {
            // Copied from surfaceCheck
            validTri(*this, triI);
        }
    }


    if (regularise)
    {
        List<FixedList<label, 3>> faceEdges;
        labelList edgeFace0, edgeFace1;
        Map<labelList> edgeFacesRest;


        while (true)
        {
            // Calculate addressing
            calcAddressing
            (
                *this,
                faceEdges,
                edgeFace0,
                edgeFace1,
                edgeFacesRest
            );

            // See if any dangling triangles
            boolList keepTriangles;
            label nDangling = markDanglingTriangles
            (
                faceEdges,
                edgeFace0,
                edgeFace1,
                edgeFacesRest,
                keepTriangles
            );

            if (debug)
            {
                Pout<< "isoSurfaceCell : detected " << nDangling
                    << " dangling triangles." << endl;
            }

            if (nDangling == 0)
            {
                break;
            }

            // Create face map (new to old)
            labelList subsetTriMap(findIndices(keepTriangles, true));

            labelList subsetPointMap;
            labelList reversePointMap;
            triSurface::operator=
            (
                subsetMesh
                (
                    *this,
                    subsetTriMap,
                    reversePointMap,
                    subsetPointMap
                )
            );
            meshCells_ = labelField(meshCells_, subsetTriMap);
            inplaceRenumber(reversePointMap, triPointMergeMap_);
        }
    }
}


// ************************************************************************* //
