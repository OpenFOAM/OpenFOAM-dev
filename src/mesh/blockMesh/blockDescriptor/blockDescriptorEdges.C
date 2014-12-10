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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "blockDescriptor.H"

#include "lineEdge.H"
#include "lineDivide.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //! \cond fileScope
    //  Calculate the geometric expension factor from the expansion ratio
    inline scalar calcGexp(const scalar expRatio, const label dim)
    {
        return dim > 1 ? pow(expRatio, 1.0/(dim - 1)) : 0.0;
    }
    //! \endcond

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockDescriptor::makeBlockEdges()
{
    const label ni = meshDensity_.x();
    const label nj = meshDensity_.y();
    const label nk = meshDensity_.z();

    // these edges correspond to the "hex" cellModel

    // x-direction
    setEdge(0,  0, 1, ni);
    setEdge(1,  3, 2, ni);
    setEdge(2,  7, 6, ni);
    setEdge(3,  4, 5, ni);

    // y-direction
    setEdge(4,  0, 3, nj);
    setEdge(5,  1, 2, nj);
    setEdge(6,  5, 6, nj);
    setEdge(7,  4, 7, nj);

    // z-direction
    setEdge(8,  0, 4, nk);
    setEdge(9,  1, 5, nk);
    setEdge(10, 2, 6, nk);
    setEdge(11, 3, 7, nk);
}


void Foam::blockDescriptor::setEdge
(
    label edgeI,
    label start,
    label end,
    label dim
)
{
    // set reference to the list of labels defining the block
    const labelList& blockLabels = blockShape_;

    // set reference to global list of points
    const pointField blockPoints = blockShape_.points(blockPointField_);

    // Set the edge points/weights
    // The edge is a straight-line if it is not in the list of curvedEdges

    // calc geometric expension factor from the expansion ratio
    const scalar gExp = calcGexp(expand_[edgeI], dim);

    forAll(curvedEdges_, cedgeI)
    {
        const curvedEdge& cedge = curvedEdges_[cedgeI];

        int cmp = cedge.compare(blockLabels[start], blockLabels[end]);

        if (cmp)
        {
            if (cmp > 0)
            {
                // curve has the same orientation

                // divide the line
                lineDivide divEdge(cedge, dim, gExp);

                edgePoints_[edgeI]  = divEdge.points();
                edgeWeights_[edgeI] = divEdge.lambdaDivisions();
            }
            else
            {
                // curve has the opposite orientation

                // divide the line
                lineDivide divEdge(cedge, dim, 1.0/(gExp+SMALL));

                const pointField& p = divEdge.points();
                const scalarList& d = divEdge.lambdaDivisions();

                edgePoints_[edgeI].setSize(p.size());
                edgeWeights_[edgeI].setSize(d.size());

                label pMax = p.size() - 1;
                forAll(p, pI)
                {
                    edgePoints_[edgeI][pI]  = p[pMax - pI];
                    edgeWeights_[edgeI][pI] = 1.0 - d[pMax - pI];
                }

            }

            // found curved-edge: done
            return;
        }
    }


    // not found: divide the edge as a straight line

    lineDivide divEdge
    (
        lineEdge(blockPoints, start, end),
        dim,
        gExp
    );

    edgePoints_[edgeI]  = divEdge.points();
    edgeWeights_[edgeI] = divEdge.lambdaDivisions();
}


// ************************************************************************* //
