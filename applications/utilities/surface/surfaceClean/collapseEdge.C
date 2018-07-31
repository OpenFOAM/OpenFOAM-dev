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

#include "collapseEdge.H"

static void markPointNbrs
(
    const triSurface& surf,
    const label facei,
    const bool val,
    boolList& okToCollapse
)
{
    const triSurface::FaceType& f = surf.localFaces()[facei];

    forAll(f, fp)
    {
        const labelList& pFaces = surf.pointFaces()[f[fp]];

        forAll(pFaces, i)
        {
            okToCollapse[pFaces[i]] = false;
        }
    }
}


static triSurface pack
(
    const triSurface& surf,
    const pointField& localPoints,
    const labelList& pointMap
)
{
    List<labelledTri> newTriangles(surf.size());
    label newTriangleI = 0;

    forAll(surf, facei)
    {
        const labelledTri& f = surf.localFaces()[facei];

        label newA = pointMap[f[0]];
        label newB = pointMap[f[1]];
        label newC = pointMap[f[2]];

        if ((newA != newB) && (newA != newC) && (newB != newC))
        {
            newTriangles[newTriangleI++] =
                labelledTri(newA, newB, newC, f.region());
        }
    }
    newTriangles.setSize(newTriangleI);

    return triSurface(newTriangles, surf.patches(), localPoints);
}


// Collapses small edge to point, thus removing triangle.
label collapseEdge(triSurface& surf, const scalar minLen)
{
    label nTotalCollapsed = 0;

    while (true)
    {
        const pointField& localPoints = surf.localPoints();
        const List<labelledTri>& localFaces = surf.localFaces();


        // Mapping from old to new points
        labelList pointMap(surf.nPoints());
        forAll(pointMap, i)
        {
            pointMap[i] = i;
        }

        // Storage for new points.
        pointField newPoints(localPoints);

        // To protect neighbours of collapsed faces.
        boolList okToCollapse(surf.size(), true);
        label nCollapsed = 0;

        forAll(localFaces, facei)
        {
            if (okToCollapse[facei])
            {
                // Check edge lengths.
                const triSurface::FaceType& f = localFaces[facei];

                forAll(f, fp)
                {
                    label v = f[fp];
                    label v1 = f[f.fcIndex(fp)];

                    if (mag(localPoints[v1] - localPoints[v]) < minLen)
                    {
                        // Collapse f[fp1] onto f[fp].
                        pointMap[v1] = v;
                        newPoints[v] = 0.5*(localPoints[v1] + localPoints[v]);

                        // Pout<< "Collapsing triangle " << facei
                        //    << " to edge mid " << newPoints[v] << endl;

                        nCollapsed++;
                        okToCollapse[facei] = false;

                        // Protect point neighbours from collapsing.
                        markPointNbrs(surf, facei, false, okToCollapse);

                        break;
                    }
                }
            }
        }

        Info<< "collapseEdge : collapsing " << nCollapsed
            << " triangles to a single edge."
            << endl;

        nTotalCollapsed += nCollapsed;

        if (nCollapsed == 0)
        {
            break;
        }

        // Pack the triangles
        surf = pack(surf, newPoints, pointMap);
    }

    // Remove any unused vertices
    surf = triSurface(surf.localFaces(), surf.patches(), surf.localPoints());

    return nTotalCollapsed;
}


// ************************************************************************* //
