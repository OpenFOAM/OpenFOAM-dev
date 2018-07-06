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

#include "regionSide.H"
#include "meshTools.H"
#include "primitiveMesh.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(regionSide, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Step across edge onto other face on cell
Foam::label Foam::regionSide::otherFace
(
    const primitiveMesh& mesh,
    const label celli,
    const label facei,
    const label edgeI
)
{
    label f0I, f1I;

    meshTools::getEdgeFaces(mesh, celli, edgeI, f0I, f1I);

    if (f0I == facei)
    {
        return f1I;
    }
    else
    {
        return f0I;
    }
}


// Step across point to other edge on face
Foam::label Foam::regionSide::otherEdge
(
    const primitiveMesh& mesh,
    const label facei,
    const label edgeI,
    const label pointi
)
{
    const edge& e = mesh.edges()[edgeI];

    // Get other point on edge.
    label freePointi = e.otherVertex(pointi);

    const labelList& fEdges = mesh.faceEdges()[facei];

    forAll(fEdges, fEdgeI)
    {
        const label otherEdgeI = fEdges[fEdgeI];
        const edge& otherE = mesh.edges()[otherEdgeI];

        if
        (
            (
                otherE.start() == pointi
             && otherE.end() != freePointi
            )
         || (
                otherE.end() == pointi
             && otherE.start() != freePointi
            )
        )
        {
            // otherE shares one (but not two) points with e.
            return otherEdgeI;
        }
    }

    FatalErrorInFunction
        << "Cannot find other edge on face " << facei << " that uses point "
        << pointi << " but not point " << freePointi << endl
        << "Edges on face:" << fEdges
        << " verts:" << UIndirectList<edge>(mesh.edges(), fEdges)()
        << " Vertices on face:"
        << mesh.faces()[facei]
        << " Vertices on original edge:" << e << abort(FatalError);

    return -1;
}


// Step from facei (on side celli) to connected face & cell without crossing
// fenceEdges.
void Foam::regionSide::visitConnectedFaces
(
    const primitiveMesh& mesh,
    const labelHashSet& region,
    const labelHashSet& fenceEdges,
    const label celli,
    const label facei,
    labelHashSet& visitedFace
)
{
    if (!visitedFace.found(facei))
    {
        if (debug)
        {
            Info<< "visitConnectedFaces : celli:" << celli << " facei:"
                << facei << "  isOwner:" << (celli == mesh.faceOwner()[facei])
                << endl;
        }

        // Mark as visited
        visitedFace.insert(facei);

        // Mark which side of face was visited.
        if (celli == mesh.faceOwner()[facei])
        {
            sideOwner_.insert(facei);
        }


        // Visit all neighbouring faces on faceSet. Stay on this 'side' of
        // face by doing edge-face-cell walk.
        const labelList& fEdges = mesh.faceEdges()[facei];

        forAll(fEdges, fEdgeI)
        {
            label edgeI = fEdges[fEdgeI];

            if (!fenceEdges.found(edgeI))
            {
                // Step along faces on edge from cell to cell until
                // we hit face on faceSet.

                // Find face reachable from edge
                label otherFacei = otherFace(mesh, celli, facei, edgeI);

                if (mesh.isInternalFace(otherFacei))
                {
                    label otherCelli = celli;

                    // Keep on crossing faces/cells until back on face on
                    // surface
                    while (!region.found(otherFacei))
                    {
                        visitedFace.insert(otherFacei);

                        if (debug)
                        {
                            Info<< "visitConnectedFaces : celli:" << celli
                                << " found insideEdgeFace:" << otherFacei
                                << endl;
                        }


                        // Cross otherFacei into neighbouring cell
                        otherCelli =
                            meshTools::otherCell
                            (
                                mesh,
                                otherCelli,
                                otherFacei
                            );

                        otherFacei =
                                otherFace
                                (
                                    mesh,
                                    otherCelli,
                                    otherFacei,
                                    edgeI
                                );
                    }

                    visitConnectedFaces
                    (
                        mesh,
                        region,
                        fenceEdges,
                        otherCelli,
                        otherFacei,
                        visitedFace
                    );
                }
            }
        }
    }
}


// From edge on face connected to point on region (regionPointi) cross
// to all other edges using this point by walking across faces
// Does not cross regionEdges so stays on one side
// of region
void Foam::regionSide::walkPointConnectedFaces
(
    const primitiveMesh& mesh,
    const labelHashSet& regionEdges,
    const label regionPointi,
    const label startFacei,
    const label startEdgeI,
    labelHashSet& visitedEdges
)
{
    // Mark as visited
    insidePointFaces_.insert(startFacei);

    if (debug)
    {
        Info<< "walkPointConnectedFaces : regionPointi:" << regionPointi
            << " facei:" << startFacei
            << " edgeI:" << startEdgeI << " verts:"
            << mesh.edges()[startEdgeI]
            << endl;
    }

    // Cross facei i.e. get edge not startEdgeI which uses regionPointi
    label edgeI = otherEdge(mesh, startFacei, startEdgeI, regionPointi);

    if (!regionEdges.found(edgeI))
    {
        if (!visitedEdges.found(edgeI))
        {
            visitedEdges.insert(edgeI);

            if (debug)
            {
                Info<< "Crossed face from "
                    << " edgeI:" << startEdgeI << " verts:"
                    << mesh.edges()[startEdgeI]
                    << " to edge:" << edgeI << " verts:"
                    << mesh.edges()[edgeI]
                    << endl;
            }

            // Cross edge to all faces connected to it.

            const labelList& eFaces = mesh.edgeFaces()[edgeI];

            forAll(eFaces, eFacei)
            {
                label facei = eFaces[eFacei];

                walkPointConnectedFaces
                (
                    mesh,
                    regionEdges,
                    regionPointi,
                    facei,
                    edgeI,
                    visitedEdges
                );
            }
        }
    }
}


// Find all faces reachable from all non-fence points and staying on
// regionFaces side.
void Foam::regionSide::walkAllPointConnectedFaces
(
    const primitiveMesh& mesh,
    const labelHashSet& regionFaces,
    const labelHashSet& fencePoints
)
{
    //
    // Get all (internal and external) edges on region.
    //
    labelHashSet regionEdges(4*regionFaces.size());

    forAllConstIter(labelHashSet, regionFaces, iter)
    {
        const label facei = iter.key();
        const labelList& fEdges = mesh.faceEdges()[facei];

        forAll(fEdges, fEdgeI)
        {
            regionEdges.insert(fEdges[fEdgeI]);
        }
    }


    //
    // Visit all internal points on surface.
    //

    // Storage for visited points
    labelHashSet visitedPoint(4*regionFaces.size());

    // Insert fence points so we don't visit them
    forAllConstIter(labelHashSet, fencePoints, iter)
    {
        visitedPoint.insert(iter.key());
    }

    labelHashSet visitedEdges(2*fencePoints.size());


    if (debug)
    {
        Info<< "Excluding visit of points:" << visitedPoint << endl;
    }

    forAllConstIter(labelHashSet, regionFaces, iter)
    {
        const label facei = iter.key();

        // Get side of face.
        label celli;

        if (sideOwner_.found(facei))
        {
            celli = mesh.faceOwner()[facei];
        }
        else
        {
            celli = mesh.faceNeighbour()[facei];
        }

        // Find starting point and edge on face.
        const labelList& fEdges = mesh.faceEdges()[facei];

        forAll(fEdges, fEdgeI)
        {
            label edgeI = fEdges[fEdgeI];

            // Get the face 'perpendicular' to facei on region.
            label otherFacei = otherFace(mesh, celli, facei, edgeI);

            // Edge
            const edge& e = mesh.edges()[edgeI];

            if (!visitedPoint.found(e.start()))
            {
                Info<< "Determining visibility from point " << e.start()
                    << endl;

                visitedPoint.insert(e.start());

                // edgeI = otherEdge(mesh, otherFacei, edgeI, e.start());

                walkPointConnectedFaces
                (
                    mesh,
                    regionEdges,
                    e.start(),
                    otherFacei,
                    edgeI,
                    visitedEdges
                );
            }
            if (!visitedPoint.found(e.end()))
            {
                Info<< "Determining visibility from point " << e.end()
                    << endl;

                visitedPoint.insert(e.end());

                // edgeI = otherEdge(mesh, otherFacei, edgeI, e.end());

                walkPointConnectedFaces
                (
                    mesh,
                    regionEdges,
                    e.end(),
                    otherFacei,
                    edgeI,
                    visitedEdges
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::regionSide::regionSide
(
    const primitiveMesh& mesh,
    const labelHashSet& region,         // faces of region
    const labelHashSet& fenceEdges,     // outside edges
    const label startCelli,
    const label startFacei
)
:
    sideOwner_(region.size()),
    insidePointFaces_(region.size())
{
    // Storage for visited faces
    labelHashSet visitedFace(region.size());

    // Visit all faces on this side of region.
    // Mark for each face whether owner (or neighbour) has been visited.
    // Sets sideOwner_
    visitConnectedFaces
    (
        mesh,
        region,
        fenceEdges,
        startCelli,
        startFacei,
        visitedFace
    );

    //
    // Visit all internal points on region and mark faces visitable from these.
    // Sets insidePointFaces_.
    //

    labelHashSet fencePoints(fenceEdges.size());

    forAllConstIter(labelHashSet, fenceEdges, iter)
    {
        const edge& e = mesh.edges()[iter.key()];

        fencePoints.insert(e.start());
        fencePoints.insert(e.end());
    }

    walkAllPointConnectedFaces(mesh, region, fencePoints);
}


// ************************************************************************* //
