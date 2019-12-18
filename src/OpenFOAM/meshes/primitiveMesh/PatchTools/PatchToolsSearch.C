/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Description
    Searching and marking zones of the patch.

\*---------------------------------------------------------------------------*/

#include "PatchTools.H"
#include "PackedBoolList.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class BoolListType, class FaceList, class PointField>
void Foam::PatchTools::markZone
(
    const PrimitivePatch<FaceList, PointField>& p,
    const BoolListType& borderEdge,
    const label facei,
    const label currentZone,
    labelList&  faceZone
)
{
    const labelListList& faceEdges = p.faceEdges();
    const labelListList& edgeFaces = p.edgeFaces();

    // List of faces whose faceZone has been set.
    labelList changedFaces(1, facei);

    while (true)
    {
        // Pick up neighbours of changedFaces
        DynamicList<label> newChangedFaces(2*changedFaces.size());

        forAll(changedFaces, i)
        {
            label facei = changedFaces[i];

            const labelList& fEdges = faceEdges[facei];

            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                if (!borderEdge[edgeI])
                {
                    const labelList& eFaceLst = edgeFaces[edgeI];

                    forAll(eFaceLst, j)
                    {
                        label nbrFacei = eFaceLst[j];

                        if (faceZone[nbrFacei] == -1)
                        {
                            faceZone[nbrFacei] = currentZone;
                            newChangedFaces.append(nbrFacei);
                        }
                        else if (faceZone[nbrFacei] != currentZone)
                        {
                            FatalErrorInFunction
                                << "Zones " << faceZone[nbrFacei]
                                << " at face " << nbrFacei
                                << " connects to zone " << currentZone
                                << " at face " << facei
                                << abort(FatalError);
                        }
                    }
                }
            }
        }

        if (newChangedFaces.empty())
        {
            break;
        }

        // transfer from dynamic to normal list
        changedFaces.transfer(newChangedFaces);
    }
}


template<class BoolListType, class FaceList, class PointField>
Foam::label Foam::PatchTools::markZones
(
    const PrimitivePatch<FaceList, PointField>& p,
    const BoolListType& borderEdge,
    labelList& faceZone
)
{
    faceZone.setSize(p.size());
    faceZone = -1;

    label zoneI = 0;
    for (label startFacei = 0; startFacei < faceZone.size();)
    {
        // Find next non-visited face
        for (; startFacei < faceZone.size(); ++startFacei)
        {
            if (faceZone[startFacei] == -1)
            {
                faceZone[startFacei] = zoneI;
                markZone(p, borderEdge, startFacei, zoneI, faceZone);
                zoneI++;
                break;
            }
        }
    }

    return zoneI;
}


template<class BoolListType, class FaceList, class PointField>
void Foam::PatchTools::subsetMap
(
    const PrimitivePatch<FaceList, PointField>& p,
    const BoolListType& includeFaces,
    labelList& pointMap,
    labelList& faceMap
)
{
    typedef typename PrimitivePatch<FaceList, PointField>::FaceType FaceType;

    label facei  = 0;
    label pointi = 0;

    const List<FaceType>& localFaces = p.localFaces();

    faceMap.setSize(localFaces.size());
    pointMap.setSize(p.nPoints());

    boolList pointHad(pointMap.size(), false);

    forAll(p, oldFacei)
    {
        if (includeFaces[oldFacei])
        {
            // Store new faces compact
            faceMap[facei++] = oldFacei;

            // Renumber labels for face
            const FaceType& f = localFaces[oldFacei];

            forAll(f, fp)
            {
                const label ptLabel = f[fp];
                if (!pointHad[ptLabel])
                {
                    pointHad[ptLabel]  = true;
                    pointMap[pointi++] = ptLabel;
                }
            }
        }
    }

    // Trim
    faceMap.setSize(facei);
    pointMap.setSize(pointi);
}


template<class FaceList, class PointField>
void Foam::PatchTools::calcBounds
(
    const PrimitivePatch<FaceList, PointField>& p,
    boundBox& bb,
    label& nPoints
)
{
    typedef typename PrimitivePatch<FaceList, PointField>::FaceType FaceType;

    // Unfortunately nPoints constructs meshPoints() so do compact version
    // ourselves
    const PointField& points = p.points();

    PackedBoolList pointIsUsed(points.size());

    nPoints = 0;
    bb = boundBox::invertedBox;

    forAll(p, facei)
    {
        const FaceType& f = p[facei];

        forAll(f, fp)
        {
            label pointi = f[fp];
            if (pointIsUsed.set(pointi, 1u))
            {
                bb.min() = ::Foam::min(bb.min(), points[pointi]);
                bb.max() = ::Foam::max(bb.max(), points[pointi]);
                nPoints++;
            }
        }
    }
}


// ************************************************************************* //
