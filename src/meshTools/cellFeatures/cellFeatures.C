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

#include "cellFeatures.H"
#include "primitiveMesh.H"
#include "HashSet.H"
#include "Map.H"
#include "demandDrivenData.H"
#include "ListOps.H"
#include "meshTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// Return true if edge start and end are on increasing face vertices. (edge is
// guaranteed to be on face)
bool Foam::cellFeatures::faceAlignedEdge(const label facei, const label edgeI)
 const
{
    const edge& e = mesh_.edges()[edgeI];

    const face& f = mesh_.faces()[facei];

    forAll(f, fp)
    {
        if (f[fp] == e.start())
        {
            label fp1 = f.fcIndex(fp);

            return f[fp1] == e.end();
        }
    }

    FatalErrorInFunction
        << "Can not find edge " << mesh_.edges()[edgeI]
        << " on face " << facei << abort(FatalError);

    return false;
}


// Return edge in featureEdge that uses vertI and is on same superface
// but is not edgeI
Foam::label Foam::cellFeatures::nextEdge
(
    const Map<label>& toSuperFace,
    const label superFacei,
    const label thisEdgeI,
    const label thisVertI
) const
{
    const labelList& pEdges = mesh_.pointEdges()[thisVertI];

    forAll(pEdges, pEdgeI)
    {
        label edgeI = pEdges[pEdgeI];

        if ((edgeI != thisEdgeI) && featureEdge_.found(edgeI))
        {
            // Check that edge is used by a face on same superFace

            const labelList& eFaces = mesh_.edgeFaces()[edgeI];

            forAll(eFaces, eFacei)
            {
                label facei = eFaces[eFacei];

                if
                (
                    meshTools::faceOnCell(mesh_, celli_, facei)
                 && (toSuperFace[facei] == superFacei)
                )
                {
                    return edgeI;
                }
            }
        }
    }

    FatalErrorInFunction
        << "Can not find edge in " << featureEdge_ << " connected to edge "
        << thisEdgeI << " at vertex " << thisVertI << endl
        << "This might mean that the externalEdges do not form a closed loop"
        << abort(FatalError);

    return -1;
}


// Return true if angle between faces using it is larger than certain value.
bool Foam::cellFeatures::isCellFeatureEdge
(
    const scalar minCos,
    const label edgeI
) const
{
    // Get the two faces using this edge

    label face0;
    label face1;
    meshTools::getEdgeFaces(mesh_, celli_, edgeI, face0, face1);

    // Check the angle between them by comparing the face normals.

    vector n0 = mesh_.faceAreas()[face0];
    n0 /= mag(n0);

    vector n1 = mesh_.faceAreas()[face1];
    n1 /= mag(n1);

    scalar cosAngle = n0 & n1;


    const edge& e = mesh_.edges()[edgeI];

    const face& f0 = mesh_.faces()[face0];

    label face0Start = findIndex(f0, e.start());
    label face0End   = f0.fcIndex(face0Start);

    const face& f1 = mesh_.faces()[face1];

    label face1Start = findIndex(f1, e.start());
    label face1End   = f1.fcIndex(face1Start);

    if
    (
        (
            (f0[face0End] == e.end())
         && (f1[face1End] != e.end())
        )
     || (
            (f0[face0End] != e.end())
         && (f1[face1End] == e.end())
        )
    )
    {
    }
    else
    {
        cosAngle = -cosAngle;
    }

    if (cosAngle < minCos)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Recursively mark (on toSuperFace) all face reachable across non-feature
// edges.
void Foam::cellFeatures::walkSuperFace
(
    const label facei,
    const label superFacei,
    Map<label>& toSuperFace
) const
{
    if (!toSuperFace.found(facei))
    {
        toSuperFace.insert(facei, superFacei);

        const labelList& fEdges = mesh_.faceEdges()[facei];

        forAll(fEdges, fEdgeI)
        {
            label edgeI = fEdges[fEdgeI];

            if (!featureEdge_.found(edgeI))
            {
                label face0;
                label face1;
                meshTools::getEdgeFaces(mesh_, celli_, edgeI, face0, face1);

                if (face0 == facei)
                {
                    face0 = face1;
                }

                walkSuperFace
                (
                    face0,
                    superFacei,
                    toSuperFace
                );
            }
        }
    }
}


void Foam::cellFeatures::calcSuperFaces() const
{
    // Determine superfaces by edge walking across non-feature edges

    const labelList& cFaces = mesh_.cells()[celli_];

    // Mapping from old to super face:
    //    <not found> : not visited
    //    >=0 : superFace
    Map<label> toSuperFace(10*cFaces.size());

    label superFacei = 0;

    forAll(cFaces, cFacei)
    {
        label facei = cFaces[cFacei];

        if (!toSuperFace.found(facei))
        {
            walkSuperFace
            (
                facei,
                superFacei,
                toSuperFace
            );
            superFacei++;
        }
    }

    // Construct superFace-to-oldface mapping.

    faceMap_.setSize(superFacei);

    forAll(cFaces, cFacei)
    {
        label facei = cFaces[cFacei];

        faceMap_[toSuperFace[facei]].append(facei);
    }

    forAll(faceMap_, superI)
    {
        faceMap_[superI].shrink();
    }


    // Construct superFaces

    facesPtr_ = new faceList(superFacei);

    faceList& faces = *facesPtr_;

    forAll(cFaces, cFacei)
    {
        label facei = cFaces[cFacei];

        label superFacei = toSuperFace[facei];

        if (faces[superFacei].empty())
        {
            // Superface not yet constructed.

            // Find starting feature edge on face.
            label startEdgeI = -1;

            const labelList& fEdges = mesh_.faceEdges()[facei];

            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                if (featureEdge_.found(edgeI))
                {
                    startEdgeI = edgeI;

                    break;
                }
            }


            if (startEdgeI != -1)
            {
                // Walk point-edge-point along feature edges

                DynamicList<label> superFace(10*mesh_.faces()[facei].size());

                const edge& e = mesh_.edges()[startEdgeI];

                // Walk either start-end or end-start depending on orientation
                // of face. SuperFace will have celli as owner.
                bool flipOrientation =
                    (mesh_.faceOwner()[facei] == celli_)
                  ^ (faceAlignedEdge(facei, startEdgeI));

                label startVertI = -1;

                if (flipOrientation)
                {
                    startVertI = e.end();
                }
                else
                {
                    startVertI = e.start();
                }

                label edgeI = startEdgeI;

                label vertI = e.otherVertex(startVertI);

                do
                {
                    label newEdgeI = nextEdge
                    (
                        toSuperFace,
                        superFacei,
                        edgeI,
                        vertI
                    );

                    // Determine angle between edges.
                    if (isFeaturePoint(edgeI, newEdgeI))
                    {
                        superFace.append(vertI);
                    }

                    edgeI = newEdgeI;

                    if (vertI == startVertI)
                    {
                        break;
                    }

                    vertI = mesh_.edges()[edgeI].otherVertex(vertI);
                }
                while (true);

                if (superFace.size() <= 2)
                {
                    WarningInFunction
                        << " Can not collapse faces " << faceMap_[superFacei]
                        << " into one big face on cell " << celli_ << endl
                        << "Try decreasing minCos:" << minCos_ << endl;
                }
                else
                {
                    faces[superFacei].transfer(superFace);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::cellFeatures::cellFeatures
(
    const primitiveMesh& mesh,
    const scalar minCos,
    const label celli
)
:
    mesh_(mesh),
    minCos_(minCos),
    celli_(celli),
    featureEdge_(10*mesh.cellEdges()[celli].size()),
    facesPtr_(nullptr),
    faceMap_(0)
{
    const labelList& cEdges = mesh_.cellEdges()[celli_];

    forAll(cEdges, cEdgeI)
    {
        label edgeI = cEdges[cEdgeI];

        if (isCellFeatureEdge(minCos_, edgeI))
        {
            featureEdge_.insert(edgeI);
        }
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellFeatures::~cellFeatures()
{
    deleteDemandDrivenData(facesPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cellFeatures::isFeaturePoint(const label edge0, const label edge1)
 const
{
    if
    (
        (edge0 < 0)
     || (edge0 >= mesh_.nEdges())
     || (edge1 < 0)
     || (edge1 >= mesh_.nEdges())
    )
    {
        FatalErrorInFunction
            << "Illegal edge labels : edge0:" << edge0 << " edge1:" << edge1
            << abort(FatalError);
    }

    const edge& e0 = mesh_.edges()[edge0];

    vector e0Vec = e0.vec(mesh_.points());
    e0Vec /= mag(e0Vec);

    const edge& e1 = mesh_.edges()[edge1];

    vector e1Vec = e1.vec(mesh_.points());
    e1Vec /= mag(e1Vec);

    scalar cosAngle;

    if
    (
        (e0.start() == e1.end())
     || (e0.end() == e1.start())
    )
    {
        // Same direction
        cosAngle = e0Vec & e1Vec;
    }
    else if
    (
        (e0.start() == e1.start())
     || (e0.end() == e1.end())
    )
    {
        // back on back
        cosAngle = - e0Vec & e1Vec;
    }
    else
    {
        cosAngle = great;   // satisfy compiler

        FatalErrorInFunction
            << "Edges do not share common vertex. e0:" << e0
            << " e1:" << e1 << abort(FatalError);
    }

    if (cosAngle < minCos_)
    {
        // Angle larger than criterium
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::cellFeatures::isFeatureVertex(const label facei, const label vertI)
 const
{
    if
    (
        (facei < 0)
     || (facei >= mesh_.nFaces())
     || (vertI < 0)
     || (vertI >= mesh_.nPoints())
    )
    {
        FatalErrorInFunction
            << "Illegal face " << facei << " or vertex " << vertI
            << abort(FatalError);
    }

    const labelList& pEdges = mesh_.pointEdges()[vertI];

    label edge0 = -1;
    label edge1 = -1;

    forAll(pEdges, pEdgeI)
    {
        label edgeI = pEdges[pEdgeI];

        if (meshTools::edgeOnFace(mesh_, facei, edgeI))
        {
            if (edge0 == -1)
            {
                edge0 = edgeI;
            }
            else
            {
                edge1 = edgeI;

                // Found the two edges.
                break;
            }
        }
    }

    if (edge1 == -1)
    {
        FatalErrorInFunction
            << "Did not find two edges sharing vertex " << vertI
            << " on face " << facei << " vertices:" << mesh_.faces()[facei]
            << abort(FatalError);
    }

    return isFeaturePoint(edge0, edge1);
}


// ************************************************************************* //
