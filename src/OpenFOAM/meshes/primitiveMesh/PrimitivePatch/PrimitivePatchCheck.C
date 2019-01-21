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
    Checks topology of the patch.

\*---------------------------------------------------------------------------*/

#include "PrimitivePatch.H"
#include "Map.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::visitPointRegion
(
    const label pointi,
    const labelList& pFaces,
    const label startFacei,
    const label startEdgeI,
    boolList& pFacesHad
) const
{
    label index = findIndex(pFaces, startFacei);

    if (!pFacesHad[index])
    {
        // Mark face as been visited.
        pFacesHad[index] = true;

        // Step to next edge on face which is still using pointi
        const labelList& fEdges = faceEdges()[startFacei];

        label nextEdgeI = -1;

        forAll(fEdges, i)
        {
            label edgeI = fEdges[i];

            const edge& e = edges()[edgeI];

            if (edgeI != startEdgeI && (e[0] == pointi || e[1] == pointi))
            {
                nextEdgeI = edgeI;

                break;
            }
        }

        if (nextEdgeI == -1)
        {
            FatalErrorInFunction
                << "Problem: cannot find edge out of " << fEdges
                << "on face " << startFacei << " that uses point " << pointi
                << " and is not edge " << startEdgeI << abort(FatalError);
        }

        // Walk to next face(s) across edge.
        const labelList& eFaces = edgeFaces()[nextEdgeI];

        forAll(eFaces, i)
        {
            if (eFaces[i] != startFacei)
            {
                visitPointRegion
                (
                    pointi,
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

template<class FaceList, class PointField>
typename
Foam::PrimitivePatch<FaceList, PointField>::surfaceTopo
Foam::PrimitivePatch<FaceList, PointField>::
surfaceType() const
{
    if (debug)
    {
        InfoInFunction << "Calculating patch topology" << endl;
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
        Info<< "    Finished." << endl;
    }

    return pType;
}


template<class FaceList, class PointField>
bool
Foam::PrimitivePatch<FaceList, PointField>::
checkTopology
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        InfoInFunction << "Checking patch topology" << endl;
    }

    // Check edgeFaces

    const labelListList& edgeFcs = edgeFaces();

    bool illegalTopo = false;

    forAll(edgeFcs, edgeI)
    {
        label nNbrs = edgeFcs[edgeI].size();

        if (nNbrs < 1 || nNbrs > 2)
        {
            illegalTopo = true;

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
    }

    if (debug)
    {
        Info<< "    Finished." << endl;
    }

    return illegalTopo;
}


template<class FaceList, class PointField>
bool
Foam::PrimitivePatch<FaceList, PointField>::
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

    forAll(pf, pointi)
    {
        const labelList& pFaces = pf[pointi];

        // Visited faces (as indices into pFaces)
        boolList pFacesHad(pFaces.size(), false);

        // Starting edge
        const labelList& pEdges = pe[pointi];
        label startEdgeI = pEdges[0];

        const labelList& eFaces = ef[startEdgeI];

        forAll(eFaces, i)
        {
            // Visit all faces using pointi, starting from eFaces[i] and
            // startEdgeI. Mark off all faces visited in pFacesHad.
            this->visitPointRegion
            (
                pointi,
                pFaces,
                eFaces[i],  // starting face for walk
                startEdgeI, // starting edge for walk
                pFacesHad
            );
        }

        // After this all faces using pointi should have been visited and
        // marked off in pFacesHad.

        label unset = findIndex(pFacesHad, false);

        if (unset != -1)
        {
            foundError = true;

            label meshPointi = mp[pointi];

            if (setPtr)
            {
                setPtr->insert(meshPointi);
            }

            if (report)
            {
                Info<< "Point " << meshPointi
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
