/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "treeDataCell.H"
#include "indexedOctree.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(treeDataCell, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::treeDataCell::update()
{
    if (cacheBb_)
    {
        bbs_.setSize(cellLabels_.size());

        forAll(cellLabels_, i)
        {
            bbs_[i] =
                treeBoundBox
                (
                    mesh_.cells()[cellLabels_[i]].bb
                    (
                        mesh_.points(),
                        mesh_.faces()
                    )
                );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::treeDataCell::treeDataCell
(
    const bool cacheBb,
    const polyMesh& mesh,
    const labelUList& cellLabels,
    const polyMesh::cellDecomposition decompMode
)
:
    mesh_(mesh),
    cellLabels_(cellLabels),
    cacheBb_(cacheBb),
    decompMode_(decompMode)
{
    update();
}


Foam::treeDataCell::treeDataCell
(
    const bool cacheBb,
    const polyMesh& mesh,
    labelList&& cellLabels,
    const polyMesh::cellDecomposition decompMode
)
:
    mesh_(mesh),
    cellLabels_(move(cellLabels)),
    cacheBb_(cacheBb),
    decompMode_(decompMode)
{
    update();
}


Foam::treeDataCell::treeDataCell
(
    const bool cacheBb,
    const polyMesh& mesh,
    const polyMesh::cellDecomposition decompMode
)
:
    mesh_(mesh),
    cellLabels_(identity(mesh_.nCells())),
    cacheBb_(cacheBb),
    decompMode_(decompMode)
{
    update();
}


Foam::treeDataCell::findNearestOp::findNearestOp
(
    const indexedOctree<treeDataCell>& tree
)
:
    tree_(tree)
{}


Foam::treeDataCell::findIntersectOp::findIntersectOp
(
    const indexedOctree<treeDataCell>& tree
)
:
    tree_(tree)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::treeDataCell::shapePoints() const
{
    pointField cc(cellLabels_.size());

    forAll(cellLabels_, i)
    {
        cc[i] = mesh_.cellCentres()[cellLabels_[i]];
    }

    return cc;
}


bool Foam::treeDataCell::overlaps
(
    const label index,
    const treeBoundBox& cubeBb
) const
{
    if (cacheBb_)
    {
        return cubeBb.overlaps(bbs_[index]);
    }
    else
    {
        return cubeBb.overlaps
        (
            mesh_.cells()[cellLabels_[index]].bb
            (
                mesh_.points(),
                mesh_.faces()
            )
        );
    }
}


bool Foam::treeDataCell::contains
(
    const label index,
    const point& sample
) const
{
    return mesh_.pointInCell(sample, cellLabels_[index], decompMode_);
}


void Foam::treeDataCell::findNearestOp::operator()
(
    const labelUList& indices,
    const point& sample,

    scalar& nearestDistSqr,
    label& minIndex,
    point& nearestPoint
) const
{
    const treeDataCell& shape = tree_.shapes();

    forAll(indices, i)
    {
        label index = indices[i];
        label celli = shape.cellLabels()[index];
        scalar distSqr = magSqr(sample - shape.mesh().cellCentres()[celli]);

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            minIndex = index;
            nearestPoint = shape.mesh().cellCentres()[celli];
        }
    }
}


void Foam::treeDataCell::findNearestOp::operator()
(
    const labelUList& indices,
    const linePointRef& ln,

    treeBoundBox& tightest,
    label& minIndex,
    point& linePoint,
    point& nearestPoint
) const
{
    NotImplemented;
}


bool Foam::treeDataCell::findIntersectOp::operator()
(
    const label index,
    const point& start,
    const point& end,
    point& intersectionPoint
) const
{
    const treeDataCell& shape = tree_.shapes();

    // Do quick rejection test
    if (shape.cacheBb_)
    {
        const treeBoundBox& cellBb = shape.bbs_[index];

        if ((cellBb.posBits(start) & cellBb.posBits(end)) != 0)
        {
            // Start and end in same block outside of cellBb.
            return false;
        }
    }
    else
    {
        const treeBoundBox cellBb
        (
            shape.mesh_.cells()[shape.cellLabels_[index]].bb
            (
                shape.mesh_.points(),
                shape.mesh_.faces()
            )
        );

        if ((cellBb.posBits(start) & cellBb.posBits(end)) != 0)
        {
            // Start and end in same block outside of cellBb.
            return false;
        }
    }


    // Do intersection with all faces of cell
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Disable picking up intersections behind us.
    scalar oldTol = intersection::setPlanarTol(0.0);

    const cell& cFaces = shape.mesh_.cells()[shape.cellLabels_[index]];

    const vector dir(end - start);
    scalar minDistSqr = magSqr(dir);
    bool hasMin = false;

    forAll(cFaces, i)
    {
        const face& f = shape.mesh_.faces()[cFaces[i]];

        pointHit inter = f.ray
        (
            start,
            dir,
            shape.mesh_.points(),
            intersection::algorithm::halfRay
        );

        if (inter.hit() && sqr(inter.distance()) <= minDistSqr)
        {
            // Note: no extra test on whether intersection is in front of us
            // since using half_ray AND zero tolerance. (note that tolerance
            // is used to look behind us)
            minDistSqr = sqr(inter.distance());
            intersectionPoint = inter.hitPoint();
            hasMin = true;
        }
    }

    // Restore picking tolerance
    intersection::setPlanarTol(oldTol);

    return hasMin;
}


// ************************************************************************* //
