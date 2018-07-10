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

#include "block.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#define w0 w[0][i]
#define w1 w[1][i]
#define w2 w[2][i]
#define w3 w[3][i]

#define w4 w[4][j]
#define w5 w[5][j]
#define w6 w[6][j]
#define w7 w[7][j]

#define w8 w[8][k]
#define w9 w[9][k]
#define w10 w[10][k]
#define w11 w[11][k]

void Foam::block::createPoints()
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
    pointField p[12];
    scalarList w[12];
    label nCurvedEdges = edgesPointsWeights(p, w);

    points_.setSize(nPoints());

    points_[pointLabel(0,  0,  0)] = p000;
    points_[pointLabel(ni, 0,  0)] = p100;
    points_[pointLabel(ni, nj, 0)] = p110;
    points_[pointLabel(0,  nj, 0)] = p010;
    points_[pointLabel(0,  0,  nk)] = p001;
    points_[pointLabel(ni, 0,  nk)] = p101;
    points_[pointLabel(ni, nj, nk)] = p111;
    points_[pointLabel(0,  nj, nk)] = p011;

    for (label k=0; k<=nk; k++)
    {
        for (label j=0; j<=nj; j++)
        {
            for (label i=0; i<=ni; i++)
            {
                // Skip block vertices
                if (vertex(i, j, k)) continue;

                const label vijk = pointLabel(i, j, k);

                // Calculate the weighting factors for all edges

                // x-direction
                scalar wx1 = (1 - w0)*(1 - w4)*(1 - w8) + w0*(1 - w5)*(1 - w9);
                scalar wx2 = (1 - w1)*w4*(1 - w11)      + w1*w5*(1 - w10);
                scalar wx3 = (1 - w2)*w7*w11            + w2*w6*w10;
                scalar wx4 = (1 - w3)*(1 - w7)*w8       + w3*(1 - w6)*w9;

                const scalar sumWx = wx1 + wx2 + wx3 + wx4;
                wx1 /= sumWx;
                wx2 /= sumWx;
                wx3 /= sumWx;
                wx4 /= sumWx;


                // y-direction
                scalar wy1 = (1 - w4)*(1 - w0)*(1 - w8) + w4*(1 - w1)*(1 - w11);
                scalar wy2 = (1 - w5)*w0*(1 - w9)       + w5*w1*(1 - w10);
                scalar wy3 = (1 - w6)*w3*w9             + w6*w2*w10;
                scalar wy4 = (1 - w7)*(1 - w3)*w8       + w7*(1 - w2)*w11;

                const scalar sumWy = wy1 + wy2 + wy3 + wy4;
                wy1 /= sumWy;
                wy2 /= sumWy;
                wy3 /= sumWy;
                wy4 /= sumWy;


                // z-direction
                scalar wz1 = (1 - w8)*(1 - w0)*(1 - w4) + w8*(1 - w3)*(1 - w7);
                scalar wz2 = (1 - w9)*w0*(1 - w5)       + w9*w3*(1 - w6);
                scalar wz3 = (1 - w10)*w1*w5            + w10*w2*w6;
                scalar wz4 = (1 - w11)*(1 - w1)*w4      + w11*(1 - w2)*w7;

                const scalar sumWz = wz1 + wz2 + wz3 + wz4;
                wz1 /= sumWz;
                wz2 /= sumWz;
                wz3 /= sumWz;
                wz4 /= sumWz;


                // Points on straight edges
                const vector edgex1 = p000 + (p100 - p000)*w0;
                const vector edgex2 = p010 + (p110 - p010)*w1;
                const vector edgex3 = p011 + (p111 - p011)*w2;
                const vector edgex4 = p001 + (p101 - p001)*w3;

                const vector edgey1 = p000 + (p010 - p000)*w4;
                const vector edgey2 = p100 + (p110 - p100)*w5;
                const vector edgey3 = p101 + (p111 - p101)*w6;
                const vector edgey4 = p001 + (p011 - p001)*w7;

                const vector edgez1 = p000 + (p001 - p000)*w8;
                const vector edgez2 = p100 + (p101 - p100)*w9;
                const vector edgez3 = p110 + (p111 - p110)*w10;
                const vector edgez4 = p010 + (p011 - p010)*w11;

                // Add the contributions
                points_[vijk] =
                (
                    wx1*edgex1 + wx2*edgex2 + wx3*edgex3 + wx4*edgex4
                  + wy1*edgey1 + wy2*edgey2 + wy3*edgey3 + wy4*edgey4
                  + wz1*edgez1 + wz2*edgez2 + wz3*edgez3 + wz4*edgez4
                )/3;


                // Apply curved-edge correction if block has curved edges
                if (nCurvedEdges)
                {
                    // Calculate the correction vectors
                    const vector corx1 = wx1*(p[0][i] - edgex1);
                    const vector corx2 = wx2*(p[1][i] - edgex2);
                    const vector corx3 = wx3*(p[2][i] - edgex3);
                    const vector corx4 = wx4*(p[3][i] - edgex4);

                    const vector cory1 = wy1*(p[4][j] - edgey1);
                    const vector cory2 = wy2*(p[5][j] - edgey2);
                    const vector cory3 = wy3*(p[6][j] - edgey3);
                    const vector cory4 = wy4*(p[7][j] - edgey4);

                    const vector corz1 = wz1*(p[8][k] - edgez1);
                    const vector corz2 = wz2*(p[9][k] - edgez2);
                    const vector corz3 = wz3*(p[10][k] - edgez3);
                    const vector corz4 = wz4*(p[11][k] - edgez4);

                    points_[vijk] +=
                    (
                        corx1 + corx2 + corx3 + corx4
                      + cory1 + cory2 + cory3 + cory4
                      + corz1 + corz2 + corz3 + corz4
                    );
                }
            }
        }
    }

    if (!nCurvedFaces()) return;

    // Apply curvature correction to face points
    FixedList<pointField, 6> facePoints(this->facePoints(points_));
    correctFacePoints(facePoints);

    // Loop over points and apply the correction from the from the i-faces
    for (label ii=0; ii<=ni; ii++)
    {
        // Process the points on the curved faces last
        label i = (ii + 1)%(ni + 1);

        for (label j=0; j<=nj; j++)
        {
            for (label k=0; k<=nk; k++)
            {
                // Skip non-curved faces and edges
                if (flatFaceOrEdge(i, j, k)) continue;

                const label vijk = pointLabel(i, j, k);

                scalar wf0 =
                (
                    (1 - w0)*(1 - w4)*(1 - w8)
                  + (1 - w1)*w4*(1 - w11)
                  + (1 - w2)*w7*w11
                  + (1 - w3)*(1 - w7)*w8
                );

                scalar wf1 =
                (
                    w0*(1 - w5)*(1 - w9)
                  + w1*w5*(1 - w10)
                  + w2*w5*w10
                  + w3*(1 - w6)*w9
                );

                const scalar sumWf = wf0 + wf1;
                wf0 /= sumWf;
                wf1 /= sumWf;

                points_[vijk] +=
                (
                    wf0
                   *(
                       facePoints[0][facePointLabel(0, j, k)]
                     - points_[pointLabel(0, j, k)]
                    )
                  + wf1
                   *(
                       facePoints[1][facePointLabel(1, j, k)]
                     - points_[pointLabel(ni, j, k)]
                    )
                );
            }
        }
    }

    // Loop over points and apply the correction from the from the j-faces
    for (label jj=0; jj<=nj; jj++)
    {
        // Process the points on the curved faces last
        label j = (jj + 1)%(nj + 1);

        for (label i=0; i<=ni; i++)
        {
            for (label k=0; k<=nk; k++)
            {
                // Skip flat faces and edges
                if (flatFaceOrEdge(i, j, k)) continue;

                const label vijk = pointLabel(i, j, k);

                scalar wf2 =
                (
                    (1 - w4)*(1 - w1)*(1 - w8)
                  + (1 - w5)*w0*(1 - w9)
                  + (1 - w6)*w3*w9
                  + (1 - w7)*(1 - w3)*w8
                );

                scalar wf3 =
                (
                    w4*(1 - w1)*(1 - w11)
                  + w5*w1*(1 - w10)
                  + w6*w2*w10
                  + w7*(1 - w2)*w11
                );

                const scalar sumWf = wf2 + wf3;
                wf2 /= sumWf;
                wf3 /= sumWf;

                points_[vijk] +=
                (
                    wf2
                   *(
                       facePoints[2][facePointLabel(2, i, k)]
                     - points_[pointLabel(i, 0, k)]
                    )
                  + wf3
                   *(
                       facePoints[3][facePointLabel(3, i, k)]
                     - points_[pointLabel(i, nj, k)]
                    )
                );
            }
        }
    }

    // Loop over points and apply the correction from the from the k-faces
    for (label kk=0; kk<=nk; kk++)
    {
        // Process the points on the curved faces last
        label k = (kk + 1)%(nk + 1);

        for (label i=0; i<=ni; i++)
        {
            for (label j=0; j<=nj; j++)
            {
                // Skip flat faces and edges
                if (flatFaceOrEdge(i, j, k)) continue;

                const label vijk = pointLabel(i, j, k);

                scalar wf4 =
                (
                    (1 - w8)*(1 - w0)*(1 - w4)
                  + (1 - w9)*w0*(1 - w5)
                  + (1 - w10)*w1*w5
                  + (1 - w11)*(1 - w1)*w4
                );

                scalar wf5 =
                (
                    w8*(1 - w3)*(1 - w7)
                  + w9*w3*(1 - w6)
                  + w10*w2*w6
                  + w11*(1 - w2)*w7
                );

                const scalar sumWf = wf4 + wf5;
                wf4 /= sumWf;
                wf5 /= sumWf;

                points_[vijk] +=
                (
                    wf4
                   *(
                       facePoints[4][facePointLabel(4, i, j)]
                     - points_[pointLabel(i, j, 0)]
                    )
                  + wf5
                   *(
                       facePoints[5][facePointLabel(5, i, j)]
                     - points_[pointLabel(i, j, nk)]
                    )
                );
            }
        }
    }
}


Foam::List<Foam::FixedList<Foam::label, 8>> Foam::block::cells() const
{
    const label ni = density().x();
    const label nj = density().y();
    const label nk = density().z();

    List<FixedList<label, 8>> cells(nCells());

    label celli = 0;

    for (label k=0; k<nk; k++)
    {
        for (label j=0; j<nj; j++)
        {
            for (label i=0; i<ni; i++)
            {
                cells[celli][0] = pointLabel(i,   j,   k);
                cells[celli][1] = pointLabel(i+1, j,   k);
                cells[celli][2] = pointLabel(i+1, j+1, k);
                cells[celli][3] = pointLabel(i,   j+1, k);
                cells[celli][4] = pointLabel(i,   j,   k+1);
                cells[celli][5] = pointLabel(i+1, j,   k+1);
                cells[celli][6] = pointLabel(i+1, j+1, k+1);
                cells[celli][7] = pointLabel(i,   j+1, k+1);

                celli++;
            }
        }
    }

    return cells;
}


void Foam::block::createBoundary()
{
    const label ni = density().x();
    const label nj = density().y();
    const label nk = density().z();

    label patchi = 0;
    label facei = 0;

    // x-direction

    // x-min
    boundaryPatches_[patchi].setSize(nj*nk);
    for (label k=0; k<nk; k++)
    {
        for (label j=0; j<nj; j++)
        {
            boundaryPatches_[patchi][facei][0] = pointLabel(0, j,   k);
            boundaryPatches_[patchi][facei][1] = pointLabel(0, j,   k+1);
            boundaryPatches_[patchi][facei][2] = pointLabel(0, j+1, k+1);
            boundaryPatches_[patchi][facei][3] = pointLabel(0, j+1, k);

            facei++;
        }
    }

    // x-max
    patchi++;
    facei = 0;

    boundaryPatches_[patchi].setSize(nj*nk);

    for (label k=0; k<nk; k++)
    {
        for (label j=0; j<nj; j++)
        {
            boundaryPatches_[patchi][facei][0] = pointLabel(ni, j,   k);
            boundaryPatches_[patchi][facei][1] = pointLabel(ni, j+1, k);
            boundaryPatches_[patchi][facei][2] = pointLabel(ni, j+1, k+1);
            boundaryPatches_[patchi][facei][3] = pointLabel(ni, j,   k+1);

            facei++;
        }
    }

    // y-direction

    // y-min
    patchi++;
    facei = 0;

    boundaryPatches_[patchi].setSize(ni*nk);
    for (label i=0; i<ni; i++)
    {
        for (label k=0; k<nk; k++)
        {
            boundaryPatches_[patchi][facei][0] = pointLabel(i,   0, k);
            boundaryPatches_[patchi][facei][1] = pointLabel(i+1, 0, k);
            boundaryPatches_[patchi][facei][2] = pointLabel(i+1, 0, k+1);
            boundaryPatches_[patchi][facei][3] = pointLabel(i,   0, k+1);

            facei++;
        }
    }

    // y-max
    patchi++;
    facei = 0;

    boundaryPatches_[patchi].setSize(ni*nk);

    for (label i=0; i<ni; i++)
    {
        for (label k=0; k<nk; k++)
        {
            boundaryPatches_[patchi][facei][0] = pointLabel(i,   nj, k);
            boundaryPatches_[patchi][facei][1] = pointLabel(i,   nj, k+1);
            boundaryPatches_[patchi][facei][2] = pointLabel(i+1, nj, k+1);
            boundaryPatches_[patchi][facei][3] = pointLabel(i+1, nj, k);

            facei++;
        }
    }

    // z-direction

    // z-min
    patchi++;
    facei = 0;

    boundaryPatches_[patchi].setSize(ni*nj);

    for (label i=0; i<ni; i++)
    {
        for (label j=0; j<nj; j++)
        {
            boundaryPatches_[patchi][facei][0] = pointLabel(i,   j,   0);
            boundaryPatches_[patchi][facei][1] = pointLabel(i,   j+1, 0);
            boundaryPatches_[patchi][facei][2] = pointLabel(i+1, j+1, 0);
            boundaryPatches_[patchi][facei][3] = pointLabel(i+1, j,   0);

            facei++;
        }
    }

    // z-max
    patchi++;
    facei = 0;

    boundaryPatches_[patchi].setSize(ni*nj);

    for (label i=0; i<ni; i++)
    {
        for (label j=0; j<nj; j++)
        {
            boundaryPatches_[patchi][facei][0] = pointLabel(i,   j,   nk);
            boundaryPatches_[patchi][facei][1] = pointLabel(i+1, j,   nk);
            boundaryPatches_[patchi][facei][2] = pointLabel(i+1, j+1, nk);
            boundaryPatches_[patchi][facei][3] = pointLabel(i,   j+1, nk);

            facei++;
        }
    }
}


// ************************************************************************* //
