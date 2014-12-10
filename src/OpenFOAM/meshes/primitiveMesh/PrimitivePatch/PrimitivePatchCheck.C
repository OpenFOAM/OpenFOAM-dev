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

Description
    Checks topology of the patch.

\*---------------------------------------------------------------------------*/

#include "PrimitivePatch.H"
#include "Map.H"
#include "ListOps.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
visitPointRegion
(
    const label pointI,
    const labelList& pFaces,
    const label startFaceI,
    const label startEdgeI,
    boolList& pFacesHad
) const
{
    label index = findIndex(pFaces, startFaceI);

    if (!pFacesHad[index])
    {
        // Mark face as been visited.
        pFacesHad[index] = true;

        // Step to next edge on face which is still using pointI
        const labelList& fEdges = faceEdges()[startFaceI];

        label nextEdgeI = -1;

        forAll(fEdges, i)
        {
            label edgeI = fEdges[i];

            const edge& e = edges()[edgeI];

            if (edgeI != startEdgeI && (e[0] == pointI || e[1] == pointI))
            {
                nextEdgeI = edgeI;

                break;
            }
        }

        if (nextEdgeI == -1)
        {
            FatalErrorIn
            (
                "PrimitivePatch<Face, FaceList, PointField, PointType>::"
                "visitPointRegion"
            )   << "Problem: cannot find edge out of " << fEdges
                << "on face " << startFaceI << " that uses point " << pointI
                << " and is not edge " << startEdgeI << abort(FatalError);
        }

        // Walk to next face(s) across edge.
        const labelList& eFaces = edgeFaces()[nextEdgeI];

        forAll(eFaces, i)
        {
            if (eFaces[i] != startFaceI)
            {
                visitPointRegion
                (
                    pointI,
                    pFaces,
                    eFaces[i],
                    nextEdgeI,
                    pFacesHad
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
typename
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::surfaceTopo
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
surfaceType() const
{
    if (debug)
    {
        Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "surfaceType() : "
               "calculating patch topology"
            << endl;
    }

    const labelListList& edgeFcs = edgeFaces();

    surfaceTopo pType = MANIFOLD;

    forAll(edgeFcs, edgeI)
    {
        label nNbrs = edgeFcs[edgeI].size();

        if (nNbrs < 1 || nNbrs > 2)
        {
            pType = ILLEGAL;

            // Can exit now. Surface is illegal.
            return pType;
        }
        else if (nNbrs == 1)
        {
            // Surface might be open or illegal so keep looping.
            pType = OPEN;
        }
    }

    if (debug)
    {
        Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "surfaceType() : "
               "finished calculating patch topology"
            << endl;
    }

    return pType;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
bool
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
checkTopology
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "checkTopology(const bool, labelHashSet&) : "
               "checking patch topology"
            << endl;
    }

    // Check edgeFaces

    const labelListList& edgeFcs = edgeFaces();

    surfaceTopo surfaceType = MANIFOLD;

    forAll(edgeFcs, edgeI)
    {
        label nNbrs = edgeFcs[edgeI].size();

        if (nNbrs < 1 || nNbrs > 2)
        {
            surfaceType = ILLEGAL;

            if (report)
            {
                Info<< "Edge " << edgeI << " with vertices:" << edges()[edgeI]
                    << " has " << nNbrs << " face neighbours"
                    << endl;
            }

            if (setPtr)
            {
                const edge& e = edges()[edgeI];

                setPtr->insert(meshPoints()[e.start()]);
                setPtr->insert(meshPoints()[e.end()]);
            }
        }
        else if (nNbrs == 1)
        {
            surfaceType = OPEN;
        }
    }

    if (debug)
    {
        Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "checkTopology(const bool, labelHashSet&) : "
               "finished checking patch topology"
            << endl;
    }

    return surfaceType == ILLEGAL;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
bool
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
checkPointManifold
(
    const bool report,
    labelHashSet* setPtr
) const
{
    const labelListList& pf = pointFaces();
    const labelListList& pe = pointEdges();
    const labelListList& ef = edgeFaces();
    const labelList& mp = meshPoints();

    bool foundError = false;

    forAll(pf, pointI)
    {
        const labelList& pFaces = pf[pointI];

        // Visited faces (as indices into pFaces)
        boolList pFacesHad(pFaces.size(), false);

        // Starting edge
        const labelList& pEdges = pe[pointI];
        label startEdgeI = pEdges[0];

        const labelList& eFaces = ef[startEdgeI];

        forAll(eFaces, i)
        {
            // Visit all faces using pointI, starting from eFaces[i] and
            // startEdgeI. Mark off all faces visited in pFacesHad.
            this->visitPointRegion
            (
                pointI,
                pFaces,
                eFaces[i],  // starting face for walk
                startEdgeI, // starting edge for walk
                pFacesHad
            );
        }

        // After this all faces using pointI should have been visited and
        // marked off in pFacesHad.

        label unset = findIndex(pFacesHad, false);

        if (unset != -1)
        {
            foundError = true;

            label meshPointI = mp[pointI];

            if (setPtr)
            {
                setPtr->insert(meshPointI);
            }

            if (report)
            {
                Info<< "Point " << meshPointI
                    << " uses faces which are not connected through an edge"
                    << nl
                    << "This means that the surface formed by this patched"
                    << " is multiply connected at this point" << nl
                    << "Connected (patch) faces:" << nl;

                forAll(pFacesHad, i)
                {
                    if (pFacesHad[i])
                    {
                        Info<< "    " << pFaces[i] << endl;
                    }
                }

                Info<< nl << "Unconnected (patch) faces:" << nl;
                forAll(pFacesHad, i)
                {
                    if (!pFacesHad[i])
                    {
                        Info<< "    " << pFaces[i] << endl;
                    }
                }
            }
        }
    }

    return foundError;
}


// ************************************************************************* //
