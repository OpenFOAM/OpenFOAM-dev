/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "sweepFace.H"
#include "blockDescriptor.H"
#include "transform.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockFaces
{
    defineTypeNameAndDebug(sweepFace, 0);
    addToRunTimeSelectionTable(blockFace, sweepFace, Istream);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockFaces::sweepFace::sweepFace
(
    const dictionary& dict,
    const label index,
    const searchableSurfaceList& geometry,
    Istream& is
)
:
    blockFace(dict, index, is),
    rotationFraction_(readScalar(is))
{
    if (rotationFraction_ < 0 || rotationFraction_ > 1)
    {
        FatalIOErrorInFunction(is)
            << "Rotation fraction for face " << vertices()
            << " must be between 0 and 1"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockFaces::sweepFace::project
(
    const blockDescriptor& desc,
    const label blockFacei,
    pointField& points
) const
{
    // Determine the extents of the two-dimensional grid of points
    static const List<direction> blockFaceDirection0({2, 2, 2, 2, 1, 1});
    static const List<direction> blockFaceDirection1({1, 1, 0, 0, 0, 0});
    label nRows = desc.density()[blockFaceDirection0[blockFacei]] + 1;
    label nCols = desc.density()[blockFaceDirection1[blockFacei]] + 1;

    // Based on the ordering of the user-supplied vertices, determine whether
    // or not the grid is to be transposed, and therefore in which direction
    // across the grid to sweep
    bool transpose = false;
    {
        static const labelList blockFaceAnchorEdgeis({4, 5, 0, 1, 0, 3});
        const edge anchorEdge =
            desc.blockShape().edges()[blockFaceAnchorEdgeis[blockFacei]];
        forAll(vertices(), i)
        {
            if (edge::compare(anchorEdge, vertices().faceEdge(i)) != 0)
            {
                transpose = i % 2;
                break;
            }
        }
        if (transpose) Swap(nRows, nCols);
    }

    // Access a point by index. The index may be negative, in which case it
    // will index in a reverse direction from the end, a-la python.
    auto pi = [&](label i, label j)
    {
        if (i < 0) i = nRows + i;
        if (j < 0) j = nCols + j;
        return transpose ? j*nRows + i : i*nCols + j;
    };
    auto p = [&](label i, label j) -> point&
    {
        return points[pi(i, j)];
    };

    // Determine the rotations *around* the edge for every row
    scalarList dThetaAround(nRows);
    for (label i = 1; i < nRows - 1; ++ i)
    {
        const vector axisi = normalised(p(i,-1) - p(i,0));
        const tensor A = tensor::I - sqr(axisi);
        const face fPrev({pi(i,0), pi(i,-1), pi(i-1,-1), pi(i-1,0)});
        const face fNext({pi(i+1,0), pi(i+1,-1), pi(i,-1), pi(i,0)});
        const vector aPrev = face::area(UIndirectList<point>(points, fPrev));
        const vector aNext = face::area(UIndirectList<point>(points, fNext));
        dThetaAround[i] =
            asin((normalised(A & aPrev) ^ normalised(A & aNext)) & axisi);
    }
    scalarList thetaAround0(nRows), thetaAround1(nRows);
    thetaAround0[0] = thetaAround1[nRows-1] = 0;
    for (label i = 1; i < nRows - 1; ++ i)
    {
        thetaAround0[i] = thetaAround0[i-1] + dThetaAround[i];
    }
    for (label i = nRows - 2; i > 0; -- i)
    {
        thetaAround1[i] = thetaAround1[i+1] - dThetaAround[i];
    }

    // Calculate the positions of the points in each row in turn
    const vector axis0 = normalised(p(0,-1) - p(0,0));
    const vector axis1 = normalised(p(-1,-1) - p(-1,0));
    for (label i = 1; i < nRows - 1; ++ i)
    {
        // Rotation *of* the edge
        const tensor Rof0 =
            rotationTensor
            (
                normalised(p(0,-1) - p(0,0)),
                normalised(p(i,-1) - p(i,0))
            );
        const tensor Rof1 =
            rotationTensor
            (
                normalised(p(-1,-1) - p(-1,0)),
                normalised(p(i,-1) - p(i,0))
            );

        // Rotation *around* the edge
        const tensor Raround0 = Ra(axis0, rotationFraction_*thetaAround0[i]);
        const tensor Raround1 = Ra(axis1, rotationFraction_*thetaAround1[i]);

        // Scale factor
        const scalar scale0 = mag(p(i,-1) - p(i,0))/mag(p(0,-1) - p(0,0));
        const scalar scale1 = mag(p(i,-1) - p(i,0))/mag(p(-1,-1) - p(-1,0));

        // Combined transforms
        const tensor T0 = Rof0 & Raround0 * scale0;
        const tensor T1 = Rof1 & Raround1 * scale1;

        // Set the positions in the row
        for (label j = 1; j < nCols - 1; ++ j)
        {
            const scalar fi = scalar(i)/(nRows - 1);
            const scalar fj = scalar(j)/(nCols - 1);

            p(i,j) =
                (1 - fi)*(1 - fj)*(p(i,0) + (T0 & (p(0,j) - p(0,0))))
              + (1 - fi)*fj*(p(i,-1) + (T0 & (p(0,j) - p(0,-1))))
              + fi*(1 - fj)*(p(i,0) + (T1 & (p(-1,j) - p(-1,0))))
              + fi*fj*(p(i,-1) + (T1 & (p(-1,j) - p(-1,-1))));
        }
    }
}


// ************************************************************************* //
