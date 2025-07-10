/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "meshBoundarySearch.H"
#include "meshSearchBoundBox.H"

/*---------------------------------------------------------------------------*\
           Class meshBoundarySearch::findUniqueIntersectOp Declaration
\*---------------------------------------------------------------------------*/

class Foam::meshBoundarySearch::findUniqueIntersectOp
:
    public treeDataFace::findIntersectOp
{
public:

    //- Reference to the boundary tree
    const indexedOctree<treeDataFace>& tree_;

    //- List of current hits
    const List<pointIndexHit>& hits_;


public:

    //- Construct from components
    findUniqueIntersectOp
    (
        const indexedOctree<treeDataFace>& tree,
        const List<pointIndexHit>& hits
    )
    :
        treeDataFace::findIntersectOp(tree),
        tree_(tree),
        hits_(hits)
    {}

    //- Calculate intersection of triangle with ray. Sets result
    //  accordingly
    bool operator()
    (
        const label index,
        const point& start,
        const point& end,
        point& intersectionPoint
    ) const
    {
        const primitiveMesh& mesh = tree_.shapes().mesh();

        // Check whether this hit has already happened. If the new face
        // index is the same as an existing hit then return no new hit. If
        // the new face shares a point with an existing hit face and the
        // line passes through both faces in the same direction, then this
        // is also assumed to be a duplicate hit.
        const label newFacei = tree_.shapes().faceLabels()[index];
        const face& newFace = mesh.faces()[newFacei];
        const scalar newDot = mesh.faceAreas()[newFacei] & (end - start);
        forAll(hits_, hiti)
        {
            const label oldFacei = hits_[hiti].index();
            const face& oldFace = mesh.faces()[oldFacei];
            const scalar oldDot =
                mesh.faceAreas()[oldFacei] & (end - start);

            if
            (
                hits_[hiti].index() == newFacei
             || (
                    newDot*oldDot > 0
                 && (labelHashSet(newFace) & labelHashSet(oldFace)).size()
                )
            )
            {
                return false;
            }
        }

        const bool hit =
            treeDataFace::findIntersectOp::operator()
            (
                index,
                start,
                end,
                intersectionPoint
            );

        return hit;
    }
};


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshBoundarySearch, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshBoundarySearch::meshBoundarySearch(const polyMesh& mesh)
:
    DemandDrivenMeshObject
    <
        polyMesh,
        DeletableMeshObject,
        meshBoundarySearch
    >(mesh),
    boundaryTree_
    (
        treeDataFace
        (
            false,
            mesh,
            identityMap
            (
                mesh.nInternalFaces(),
                mesh.nFaces() - mesh.nInternalFaces()
            )
        ),
        meshSearchBoundBox::New(mesh).bb(),
        8,
        10,
        scalar(3)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshBoundarySearch::~meshBoundarySearch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::meshBoundarySearch::findNearestBoundaryFace
(
    const point& p
) const
{
    pointIndexHit info =
        boundaryTree().findNearest
        (
            p,
            magSqr
            (
                boundaryTree().bb().max() - boundaryTree().bb().min()
            )
        );

    if (!info.hit())
    {
        info = boundaryTree().findNearest(p, sqr(great));
    }

    return boundaryTree().shapes().faceLabels()[info.index()];
}


Foam::pointIndexHit Foam::meshBoundarySearch::intersection
(
    const point& pStart,
    const point& pEnd
) const
{
    pointIndexHit curHit = boundaryTree().findLine(pStart, pEnd);

    if (curHit.hit())
    {
        // Change index into octreeData into face label
        curHit.setIndex(boundaryTree().shapes().faceLabels()[curHit.index()]);
    }

    return curHit;
}


Foam::List<Foam::pointIndexHit> Foam::meshBoundarySearch::intersections
(
    const point& pStart,
    const point& pEnd
) const
{
    DynamicList<pointIndexHit> hits;
    pointIndexHit curHit;

    findUniqueIntersectOp iop(boundaryTree(), hits);

    while (true)
    {
        // Get the next hit, or quit
        curHit = boundaryTree().findLine(pStart, pEnd, iop);
        if (!curHit.hit()) break;

        // Change index into octreeData into face label
        curHit.setIndex(boundaryTree().shapes().faceLabels()[curHit.index()]);

        hits.append(curHit);
    }

    hits.shrink();

    return hits;
}


bool Foam::meshBoundarySearch::isInside(const point& p) const
{
    return boundaryTree().getVolumeType(p) == volumeType::inside;
}


// ************************************************************************* //
