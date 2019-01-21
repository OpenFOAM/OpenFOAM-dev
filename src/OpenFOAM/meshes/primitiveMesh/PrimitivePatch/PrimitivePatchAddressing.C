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
    This function calculates the list of patch edges, defined on the list of
    points supporting the patch. The edges are ordered:
    - 0..nInternalEdges-1 : upper triangular order
    - nInternalEdges..    : boundary edges (no particular order)

    Other patch addressing information is also calculated:
    - faceFaces with neighbour faces in ascending order
    - edgeFaces with ascending face order
    - faceEdges sorted according to edges of a face

\*---------------------------------------------------------------------------*/

#include "PrimitivePatch.H"
#include "DynamicList.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcAddressing() const
{
    if (debug)
    {
        InfoInFunction << "calculating patch addressing" << endl;
    }

    if (edgesPtr_ || faceFacesPtr_ || edgeFacesPtr_ || faceEdgesPtr_)
    {
        // it is considered an error to attempt to recalculate
        // if already allocated
        FatalErrorInFunction
            << "addressing already calculated"
            << abort(FatalError);
    }

    // get reference to localFaces
    const List<FaceType>& locFcs = localFaces();

    // get reference to pointFaces
    const labelListList& pf = pointFaces();

    // Guess the max number of edges and neighbours for a face
    label maxEdges = 0;
    forAll(locFcs, facei)
    {
        maxEdges += locFcs[facei].size();
    }

    // create the lists for the various results. (resized on completion)
    edgesPtr_ = new edgeList(maxEdges);
    edgeList& edges = *edgesPtr_;

    edgeFacesPtr_ = new labelListList(maxEdges);
    labelListList& edgeFaces = *edgeFacesPtr_;

    // faceFaces created using a dynamic list.  Cannot guess size because
    // of multiple connections
    List<DynamicList<label>> ff(locFcs.size());

    faceEdgesPtr_ = new labelListList(locFcs.size());
    labelListList& faceEdges = *faceEdgesPtr_;

    // count the number of face neighbours
    labelList noFaceFaces(locFcs.size());

    // initialise the lists of subshapes for each face to avoid duplication
    edgeListList faceIntoEdges(locFcs.size());

    forAll(locFcs, facei)
    {
        faceIntoEdges[facei] = locFcs[facei].edges();

        labelList& curFaceEdges = faceEdges[facei];
        curFaceEdges.setSize(faceIntoEdges[facei].size());

        forAll(curFaceEdges, faceEdgeI)
        {
            curFaceEdges[faceEdgeI] = -1;
        }
    }

    // This algorithm will produce a separated list of edges, internal edges
    // starting from 0 and boundary edges starting from the top and
    // growing down.

    label nEdges = 0;

    bool found = false;

    // Note that faceIntoEdges is sorted acc. to local vertex numbering
    // in face (i.e. curEdges[0] is edge between f[0] and f[1])

    // For all local faces ...
    forAll(locFcs, facei)
    {
        // Get reference to vertices of current face and corresponding edges.
        const FaceType& curF = locFcs[facei];
        const edgeList& curEdges = faceIntoEdges[facei];

        // Record the neighbour face.  Multiple connectivity allowed
        List<DynamicList<label>> neiFaces(curF.size());
        List<DynamicList<label>> edgeOfNeiFace(curF.size());

        label nNeighbours = 0;

        // For all edges ...
        forAll(curEdges, edgeI)
        {
            // If the edge is already detected, skip
            if (faceEdges[facei][edgeI] >= 0) continue;

            found = false;

            // Set reference to the current edge
            const edge& e = curEdges[edgeI];

            // Collect neighbours for the current face vertex.

            const labelList& nbrFaces = pf[e.start()];

            forAll(nbrFaces, nbrFacei)
            {
                // set reference to the current neighbour
                label curNei = nbrFaces[nbrFacei];

                // Reject neighbours with the lower label
                if (curNei > facei)
                {
                    // get the reference to subshapes of the neighbour
                    const edgeList& searchEdges = faceIntoEdges[curNei];

                    forAll(searchEdges, neiEdgeI)
                    {
                        if (searchEdges[neiEdgeI] == e)
                        {
                            // Match
                            found = true;

                            neiFaces[edgeI].append(curNei);
                            edgeOfNeiFace[edgeI].append(neiEdgeI);

                            // Record faceFaces both ways
                            ff[facei].append(curNei);
                            ff[curNei].append(facei);

                            // Keep searching due to multiple connectivity
                        }
                    }
                }
            } // End of neighbouring faces

            if (found)
            {
                // Register another detected internal edge
                nNeighbours++;
            }
        } // End of current edges

        // Add the edges in increasing number of neighbours.
        // Note: for multiply connected surfaces, the lower index neighbour for
        // an edge will come first.

        // Add the faces in the increasing order of neighbours
        for (label neiSearch = 0; neiSearch < nNeighbours; neiSearch++)
        {
            // Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = locFcs.size();

            forAll(neiFaces, nfI)
            {
                if (neiFaces[nfI].size() && neiFaces[nfI][0] < minNei)
                {
                    nextNei = nfI;
                    minNei = neiFaces[nfI][0];
                }
            }

            if (nextNei > -1)
            {
                // Add the face to the list of faces
                edges[nEdges] = curEdges[nextNei];

                // Set face-edge and face-neighbour-edge to current face label
                faceEdges[facei][nextNei] = nEdges;

                DynamicList<label>& cnf = neiFaces[nextNei];
                DynamicList<label>& eonf = edgeOfNeiFace[nextNei];

                // Set edge-face addressing
                labelList& curEf = edgeFaces[nEdges];
                curEf.setSize(cnf.size() + 1);
                curEf[0] = facei;

                forAll(cnf, cnfI)
                {
                    faceEdges[cnf[cnfI]][eonf[cnfI]] = nEdges;

                    curEf[cnfI + 1] = cnf[cnfI];
                }

                // Stop the neighbour from being used again
                cnf.clear();
                eonf.clear();

                // Increment number of faces counter
                nEdges++;
            }
            else
            {
                FatalErrorInFunction
                    << "Error in internal edge insertion"
                    << abort(FatalError);
            }
        }
    }

    nInternalEdges_ = nEdges;

    // Do boundary faces

    forAll(faceEdges, facei)
    {
        labelList& curEdges = faceEdges[facei];

        forAll(curEdges, edgeI)
        {
            if (curEdges[edgeI] < 0)
            {
                // Grab edge and faceEdge
                edges[nEdges] = faceIntoEdges[facei][edgeI];
                curEdges[edgeI] = nEdges;

                // Add edgeFace
                labelList& curEf = edgeFaces[nEdges];
                curEf.setSize(1);
                curEf[0] = facei;

                nEdges++;
            }
        }
    }

    // edges
    edges.setSize(nEdges);

    // edgeFaces list
    edgeFaces.setSize(nEdges);

    // faceFaces list
    faceFacesPtr_ = new labelListList(locFcs.size());
    labelListList& faceFaces = *faceFacesPtr_;

    forAll(faceFaces, facei)
    {
        faceFaces[facei].transfer(ff[facei]);
    }


    if (debug)
    {
        InfoInFunction << "Finished calculating patch addressing" << endl;
    }
}


// ************************************************************************* //
