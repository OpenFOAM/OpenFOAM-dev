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

#include "block.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::block::createPoints() const
{
    // Set local variables for mesh specification
    const label ni = density().x();
    const label nj = density().y();
    const label nk = density().z();

    const point& p000 = blockPoint(0);
    const point& p100 = blockPoint(1);
    const point& p110 = blockPoint(2);
    const point& p010 = blockPoint(3);

    const point& p001 = blockPoint(4);
    const point& p101 = blockPoint(5);
    const point& p111 = blockPoint(6);
    const point& p011 = blockPoint(7);

    // List of edge point and weighting factors
    List<point> p[12];
    scalarList w[12];
    edgesPointsWeights(p, w);

    //
    // Generate vertices
    //
    points_.clear();
    points_.setSize(nPoints());

    for (label k = 0; k <= nk; k++)
    {
        for (label j = 0; j <= nj; j++)
        {
            for (label i = 0; i <= ni; i++)
            {
                const label vertexNo = vtxLabel(i, j, k);

                // Calculate the importance factors for all edges

                // x-direction
                scalar impx1 =
                (
                    (1 - w[0][i])*(1 - w[4][j])*(1 - w[8][k])
                  + w[0][i]*(1 - w[5][j])*(1 - w[9][k])
                );

                scalar impx2 =
                (
                    (1 - w[1][i])*w[4][j]*(1 - w[11][k])
                  + w[1][i]*w[5][j]*(1 - w[10][k])
                );

                scalar impx3 =
                (
                     (1 - w[2][i])*w[7][j]*w[11][k]
                   + w[2][i]*w[6][j]*w[10][k]
                );

                scalar impx4 =
                (
                    (1 - w[3][i])*(1 - w[7][j])*w[8][k]
                  + w[3][i]*(1 - w[6][j])*w[9][k]
                );

                const scalar magImpx = impx1 + impx2 + impx3 + impx4;
                impx1 /= magImpx;
                impx2 /= magImpx;
                impx3 /= magImpx;
                impx4 /= magImpx;


                // y-direction
                scalar impy1 =
                (
                    (1 - w[4][j])*(1 - w[0][i])*(1 - w[8][k])
                  + w[4][j]*(1 - w[1][i])*(1 - w[11][k])
                );

                scalar impy2 =
                (
                    (1 - w[5][j])*w[0][i]*(1 - w[9][k])
                  + w[5][j]*w[1][i]*(1 - w[10][k])
                );

                scalar impy3 =
                (
                    (1 - w[6][j])*w[3][i]*w[9][k]
                  + w[6][j]*w[2][i]*w[10][k]
                );

                scalar impy4 =
                (
                    (1 - w[7][j])*(1 - w[3][i])*w[8][k]
                  + w[7][j]*(1 - w[2][i])*w[11][k]
                );

                const scalar magImpy = impy1 + impy2 + impy3 + impy4;
                impy1 /= magImpy;
                impy2 /= magImpy;
                impy3 /= magImpy;
                impy4 /= magImpy;


                // z-direction
                scalar impz1 =
                (
                    (1 - w[8][k])*(1 - w[0][i])*(1 - w[4][j])
                  + w[8][k]*(1 - w[3][i])*(1 - w[7][j])
                );

                scalar impz2 =
                (
                    (1 - w[9][k])*w[0][i]*(1 - w[5][j])
                  + w[9][k]*w[3][i]*(1 - w[6][j])
                );

                scalar impz3 =
                (
                    (1 - w[10][k])*w[1][i]*w[5][j]
                  + w[10][k]*w[2][i]*w[6][j]
                );

                scalar impz4 =
                (
                    (1 - w[11][k])*(1 - w[1][i])*w[4][j]
                  + w[11][k]*(1 - w[2][i])*w[7][j]
                );

                const scalar magImpz = impz1 + impz2 + impz3 + impz4;
                impz1 /= magImpz;
                impz2 /= magImpz;
                impz3 /= magImpz;
                impz4 /= magImpz;


                // Points on straight edges
                const vector edgex1 = p000 + (p100 - p000)*w[0][i];
                const vector edgex2 = p010 + (p110 - p010)*w[1][i];
                const vector edgex3 = p011 + (p111 - p011)*w[2][i];
                const vector edgex4 = p001 + (p101 - p001)*w[3][i];

                const vector edgey1 = p000 + (p010 - p000)*w[4][j];
                const vector edgey2 = p100 + (p110 - p100)*w[5][j];
                const vector edgey3 = p101 + (p111 - p101)*w[6][j];
                const vector edgey4 = p001 + (p011 - p001)*w[7][j];

                const vector edgez1 = p000 + (p001 - p000)*w[8][k];
                const vector edgez2 = p100 + (p101 - p100)*w[9][k];
                const vector edgez3 = p110 + (p111 - p110)*w[10][k];
                const vector edgez4 = p010 + (p011 - p010)*w[11][k];

                // Add the contributions
                points_[vertexNo] =
                (
                    impx1*edgex1 + impx2*edgex2 + impx3*edgex3 + impx4*edgex4
                  + impy1*edgey1 + impy2*edgey2 + impy3*edgey3 + impy4*edgey4
                  + impz1*edgez1 + impz2*edgez2 + impz3*edgez3 + impz4*edgez4
                )/3.0;


                // Calculate the correction vectors
                const vector corx1 = impx1*(p[0][i] - edgex1);
                const vector corx2 = impx2*(p[1][i] - edgex2);
                const vector corx3 = impx3*(p[2][i] - edgex3);
                const vector corx4 = impx4*(p[3][i] - edgex4);

                const vector cory1 = impy1*(p[4][j] - edgey1);
                const vector cory2 = impy2*(p[5][j] - edgey2);
                const vector cory3 = impy3*(p[6][j] - edgey3);
                const vector cory4 = impy4*(p[7][j] - edgey4);

                const vector corz1 = impz1*(p[8][k] - edgez1);
                const vector corz2 = impz2*(p[9][k] - edgez2);
                const vector corz3 = impz3*(p[10][k] - edgez3);
                const vector corz4 = impz4*(p[11][k] - edgez4);

                points_[vertexNo] +=
                (
                    corx1 + corx2 + corx3 + corx4
                  + cory1 + cory2 + cory3 + cory4
                  + corz1 + corz2 + corz3 + corz4
                );
            }
        }
    }
}


void Foam::block::createCells() const
{
    const label ni = density().x();
    const label nj = density().y();
    const label nk = density().z();

    //
    // Generate cells
    //
    cells_.clear();
    cells_.setSize(nCells());

    label cellNo = 0;

    for (label k = 0; k < nk; k++)
    {
        for (label j = 0; j < nj; j++)
        {
            for (label i = 0; i < ni; i++)
            {
                cells_[cellNo].setSize(8);

                cells_[cellNo][0] =  vtxLabel(i, j, k);
                cells_[cellNo][1] =  vtxLabel(i+1, j, k);
                cells_[cellNo][2] =  vtxLabel(i+1, j+1, k);
                cells_[cellNo][3] =  vtxLabel(i, j+1, k);
                cells_[cellNo][4] =  vtxLabel(i, j, k+1);
                cells_[cellNo][5] =  vtxLabel(i+1, j, k+1);
                cells_[cellNo][6] =  vtxLabel(i+1, j+1, k+1);
                cells_[cellNo][7] =  vtxLabel(i, j+1, k+1);
                cellNo++;
            }
        }
    }
}


void Foam::block::createBoundary() const
{
    const label ni = density().x();
    const label nj = density().y();
    const label nk = density().z();

    //
    // Generate boundaries on each side of the hex
    //
    boundaryPatches_.clear();
    boundaryPatches_.setSize(6);


    // x-direction

    label wallLabel = 0;
    label wallCellLabel = 0;

    // x-min
    boundaryPatches_[wallLabel].setSize(nj*nk);
    for (label k = 0; k < nk; k++)
    {
        for (label j = 0; j < nj; j++)
        {
            boundaryPatches_[wallLabel][wallCellLabel].setSize(4);

            // Set the points
            boundaryPatches_[wallLabel][wallCellLabel][0] =
                vtxLabel(0, j, k);
            boundaryPatches_[wallLabel][wallCellLabel][1] =
                vtxLabel(0, j, k + 1);
            boundaryPatches_[wallLabel][wallCellLabel][2] =
                vtxLabel(0, j + 1, k + 1);
            boundaryPatches_[wallLabel][wallCellLabel][3] =
                vtxLabel(0, j + 1, k);

            // Update the counter
            wallCellLabel++;
        }
    }

    // x-max
    wallLabel++;
    wallCellLabel = 0;

    boundaryPatches_[wallLabel].setSize(nj*nk);

    for (label k = 0; k < nk; k++)
    {
        for (label j = 0; j < nj; j++)
        {
            boundaryPatches_[wallLabel][wallCellLabel].setSize(4);

            // Set the points
            boundaryPatches_[wallLabel][wallCellLabel][0] =
                vtxLabel(ni, j, k);
            boundaryPatches_[wallLabel][wallCellLabel][1] =
                vtxLabel(ni, j+1, k);
            boundaryPatches_[wallLabel][wallCellLabel][2] =
                vtxLabel(ni, j+1, k+1);
            boundaryPatches_[wallLabel][wallCellLabel][3] =
                vtxLabel(ni, j, k+1);

            // Update the counter
            wallCellLabel++;
        }
    }

    // y-direction

    // y-min
    wallLabel++;
    wallCellLabel = 0;

    boundaryPatches_[wallLabel].setSize(ni*nk);
    for (label i = 0; i < ni; i++)
    {
        for (label k = 0; k < nk; k++)
        {
            boundaryPatches_[wallLabel][wallCellLabel].setSize(4);

            // Set the points
            boundaryPatches_[wallLabel][wallCellLabel][0] =
                vtxLabel(i, 0, k);
            boundaryPatches_[wallLabel][wallCellLabel][1] =
                vtxLabel(i + 1, 0, k);
            boundaryPatches_[wallLabel][wallCellLabel][2] =
                vtxLabel(i + 1, 0, k + 1);
            boundaryPatches_[wallLabel][wallCellLabel][3] =
                vtxLabel(i, 0, k + 1);

            // Update the counter
            wallCellLabel++;
        }
    }

    // y-max
    wallLabel++;
    wallCellLabel = 0;

    boundaryPatches_[wallLabel].setSize(ni*nk);

    for (label i = 0; i < ni; i++)
    {
        for (label k = 0; k < nk; k++)
        {
            boundaryPatches_[wallLabel][wallCellLabel].setSize(4);

            // Set the points
            boundaryPatches_[wallLabel][wallCellLabel][0] =
                vtxLabel(i, nj, k);
            boundaryPatches_[wallLabel][wallCellLabel][1] =
                vtxLabel(i, nj, k + 1);
            boundaryPatches_[wallLabel][wallCellLabel][2] =
                vtxLabel(i + 1, nj, k + 1);
            boundaryPatches_[wallLabel][wallCellLabel][3] =
                vtxLabel(i + 1, nj, k);

            // Update the counter
            wallCellLabel++;
        }
    }

    // z-direction

    // z-min
    wallLabel++;
    wallCellLabel = 0;

    boundaryPatches_[wallLabel].setSize(ni*nj);

    for (label i = 0; i < ni; i++)
    {
        for (label j = 0; j < nj; j++)
        {
            boundaryPatches_[wallLabel][wallCellLabel].setSize(4);

            // Set the points
            boundaryPatches_[wallLabel][wallCellLabel][0] =
                vtxLabel(i, j, 0);
            boundaryPatches_[wallLabel][wallCellLabel][1] =
                vtxLabel(i, j + 1, 0);
            boundaryPatches_[wallLabel][wallCellLabel][2] =
                vtxLabel(i + 1, j + 1, 0);
            boundaryPatches_[wallLabel][wallCellLabel][3] =
                vtxLabel(i + 1, j, 0);

            // Update the counter
            wallCellLabel++;
        }
    }

    // z-max
    wallLabel++;
    wallCellLabel = 0;

    boundaryPatches_[wallLabel].setSize(ni*nj);

    for (label i = 0; i < ni; i++)
    {
        for (label j = 0; j < nj; j++)
        {
            boundaryPatches_[wallLabel][wallCellLabel].setSize(4);

            // Set the points
            boundaryPatches_[wallLabel][wallCellLabel][0] =
                vtxLabel(i, j, nk);
            boundaryPatches_[wallLabel][wallCellLabel][1] =
                vtxLabel(i + 1, j, nk);
            boundaryPatches_[wallLabel][wallCellLabel][2] =
                vtxLabel(i + 1, j + 1, nk);
            boundaryPatches_[wallLabel][wallCellLabel][3] =
                vtxLabel(i, j + 1, nk);

            // Update the counter
            wallCellLabel++;
        }
    }
}


void Foam::block::clearGeom()
{
    points_.clear();
    cells_.clear();
    boundaryPatches_.clear();
}


// ************************************************************************* //
