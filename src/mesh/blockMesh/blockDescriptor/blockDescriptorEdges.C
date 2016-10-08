/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "blockDescriptor.H"
#include "lineEdge.H"
#include "lineDivide.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockDescriptor::makeBlockEdges()
{
    const label ni = density_.x();
    const label nj = density_.y();
    const label nk = density_.z();

    // These edges correspond to the "hex" cellModel

    // X-direction
    setEdge(0,  0, 1, ni);
    setEdge(1,  3, 2, ni);
    setEdge(2,  7, 6, ni);
    setEdge(3,  4, 5, ni);

    // Y-direction
    setEdge(4,  0, 3, nj);
    setEdge(5,  1, 2, nj);
    setEdge(6,  5, 6, nj);
    setEdge(7,  4, 7, nj);

    // Z-direction
    setEdge(8,  0, 4, nk);
    setEdge(9,  1, 5, nk);
    setEdge(10, 2, 6, nk);
    setEdge(11, 3, 7, nk);
}


void Foam::blockDescriptor::setEdge
(
    label edgei,
    label start,
    label end,
    label nDiv
)
{
    // Set reference to the list of labels defining the block
    const labelList& blockLabels = blockShape_;

    // Get list of points for this block
    const pointField blockPoints = blockShape_.points(vertices_);

    // Set the edge points/weights
    // The edge is a straight-line if it is not in the list of blockEdges

    forAll(edges_, cedgei)
    {
        const blockEdge& cedge = edges_[cedgei];

        int cmp = cedge.compare(blockLabels[start], blockLabels[end]);

        if (cmp)
        {
            if (cmp > 0)
            {
                // Curve has the same orientation

                // Divide the line
                lineDivide divEdge(cedge, nDiv, expand_[edgei]);

                edgePoints_[edgei]  = divEdge.points();
                edgeWeights_[edgei] = divEdge.lambdaDivisions();
            }
            else
            {
                // Curve has the opposite orientation

                // Divide the line
                lineDivide divEdge(cedge, nDiv, expand_[edgei].inv());

                const pointField& p = divEdge.points();
                const scalarList& d = divEdge.lambdaDivisions();

                edgePoints_[edgei].setSize(p.size());
                edgeWeights_[edgei].setSize(d.size());

                label pMax = p.size() - 1;
                forAll(p, pI)
                {
                    edgePoints_[edgei][pI]  = p[pMax - pI];
                    edgeWeights_[edgei][pI] = 1.0 - d[pMax - pI];
                }
            }

            // Found curved-edge: done
            return;
        }
    }


    // Not curved-edge: divide the edge as a straight line
    lineDivide divEdge
    (
        lineEdge(blockPoints, start, end),
        nDiv,
        expand_[edgei]
    );

    edgePoints_[edgei]  = divEdge.points();
    edgeWeights_[edgei] = divEdge.lambdaDivisions();
}


// ************************************************************************* //
