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

#include "hexBlock.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

label hexBlock::vtxLabel(label a, label b, label c) const
{
    return (a + b*(xDim_ + 1) + c*(xDim_ + 1)*(yDim_ + 1));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
hexBlock::hexBlock(const label nx, const label ny, const label nz)
:
    xDim_(nx),
    yDim_(ny),
    zDim_(nz),
    blockHandedness_(noPoints),
    points_((xDim_ + 1)*(yDim_ + 1)*(zDim_ + 1))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void hexBlock::readPoints(Istream& is)
{
    forAll(points_, i)
    {
        is  >> points_[i].x() >> points_[i].y() >> points_[i].z();
    }

    // Calculate the handedness of the block
    vector i = points_[xDim_] - points_[0];
    vector j = points_[(xDim_ + 1)*yDim_] - points_[0];
    vector k = points_[(xDim_ + 1)*(yDim_ + 1)*zDim_] - points_[0];

    if (((i ^ j) & k) > 0)
    {
        Info<< "right-handed block" << endl;
        blockHandedness_ = right;
    }
    else
    {
        Info<< "left-handed block" << endl;
        blockHandedness_ = left;
    }
}


labelListList hexBlock::blockCells() const
{
    labelListList result(xDim_*yDim_*zDim_);

    label cellNo = 0;

    if (blockHandedness_ == right)
    {
        for (label k = 0; k <= zDim_ - 1; k++)
        {
            for (label j = 0; j <= yDim_ - 1; j++)
            {
                for (label i = 0; i <= xDim_ - 1; i++)
                {
                    labelList& hexLabels = result[cellNo];
                    hexLabels.setSize(8);

                    hexLabels[0] = vtxLabel(i, j, k);
                    hexLabels[1] = vtxLabel(i+1, j, k);
                    hexLabels[2] = vtxLabel(i+1, j+1, k);
                    hexLabels[3] = vtxLabel(i, j+1, k);
                    hexLabels[4] = vtxLabel(i, j, k+1);
                    hexLabels[5] = vtxLabel(i+1, j, k+1);
                    hexLabels[6] = vtxLabel(i+1, j+1, k+1);
                    hexLabels[7] = vtxLabel(i, j+1, k+1);

                    cellNo++;
                }
            }
        }
    }
    else if (blockHandedness_ == left)
    {
        for (label k = 0; k <= zDim_ - 1; k++)
        {
            for (label j = 0; j <= yDim_ - 1; j++)
            {
                for (label i = 0; i <= xDim_ - 1; i++)
                {
                    labelList& hexLabels = result[cellNo];
                    hexLabels.setSize(8);

                    hexLabels[0] = vtxLabel(i, j, k+1);
                    hexLabels[1] = vtxLabel(i+1, j, k+1);
                    hexLabels[2] = vtxLabel(i+1, j+1, k+1);
                    hexLabels[3] = vtxLabel(i, j+1, k+1);
                    hexLabels[4] = vtxLabel(i, j, k);
                    hexLabels[5] = vtxLabel(i+1, j, k);
                    hexLabels[6] = vtxLabel(i+1, j+1, k);
                    hexLabels[7] = vtxLabel(i, j+1, k);

                    cellNo++;
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unable to determine block handedness as points "
            << "have not been read in yet"
            << abort(FatalError);
    }

    return result;
}


// Return block patch faces given direction and range limits
// From the cfx manual: direction
// 0 = solid (3-D patch),
// 1 = high i, 2 = high j, 3 = high k
// 4 = low i, 5 = low j, 6 = low k
faceList hexBlock::patchFaces(const label direc, const labelList& range) const
{
    if (range.size() != 6)
    {
        FatalErrorInFunction
            << "Invalid size of the range array: " << range.size()
            << ". Should be 6 (xMin, xMax, yMin, yMax, zMin, zMax"
            << abort(FatalError);
    }

    label xMinRange = range[0];
    label xMaxRange = range[1];
    label yMinRange = range[2];
    label yMaxRange = range[3];
    label zMinRange = range[4];
    label zMaxRange = range[5];

    faceList result(0);

    switch (direc)
    {
        case 1:
        {
            // high i = xmax

            result.setSize
            (
                (yMaxRange - yMinRange + 1)*(zMaxRange - zMinRange + 1)
            );

            label p = 0;
            for (label k = zMinRange - 1; k <= zMaxRange - 1; k++)
            {
                for (label j = yMinRange - 1; j <= yMaxRange - 1; j++)
                {
                    result[p].setSize(4);

                    // set the points
                    result[p][0] = vtxLabel(xDim_, j, k);
                    result[p][1] = vtxLabel(xDim_, j+1, k);
                    result[p][2] = vtxLabel(xDim_, j+1, k+1);
                    result[p][3] = vtxLabel(xDim_, j, k+1);

                    p++;
                }
            }

            result.setSize(p);
            break;
        }

        case 2:
        {
            // high j = ymax
            result.setSize
            (
                (xMaxRange - xMinRange + 1)*(zMaxRange - zMinRange + 1)
            );

            label p = 0;
            for (label i = xMinRange - 1; i <= xMaxRange - 1; i++)
            {
                for (label k = zMinRange - 1; k <= zMaxRange - 1; k++)
                {
                    result[p].setSize(4);

                    // set the points
                    result[p][0] = vtxLabel(i, yDim_, k);
                    result[p][1] = vtxLabel(i, yDim_, k + 1);
                    result[p][2] = vtxLabel(i + 1, yDim_, k + 1);
                    result[p][3] = vtxLabel(i + 1, yDim_, k);

                    p++;
                }
            }

            result.setSize(p);
            break;
        }

        case 3:
        {
            // high k = zmax
            result.setSize
            (
                (xMaxRange - xMinRange + 1)*(yMaxRange - yMinRange + 1)
            );

            label p = 0;
            for (label i = xMinRange - 1; i <= xMaxRange - 1; i++)
            {
                for (label j = yMinRange - 1; j <= yMaxRange - 1; j++)
                {
                    result[p].setSize(4);

                    // set the points
                    result[p][0] = vtxLabel(i, j, zDim_);
                    result[p][1] = vtxLabel(i + 1, j, zDim_);
                    result[p][2] = vtxLabel(i + 1, j + 1, zDim_);
                    result[p][3] = vtxLabel(i, j + 1, zDim_);

                    p++;
                }
            }

            result.setSize(p);
            break;
        }

        case 4:
        {
            // low i = xmin
            result.setSize
            (
                (yMaxRange - yMinRange + 1)*(zMaxRange - zMinRange + 1)
            );

            label p = 0;
            for (label k = zMinRange - 1; k <= zMaxRange - 1; k++)
            {
                for (label j = yMinRange - 1; j <= yMaxRange - 1; j++)
                {
                    result[p].setSize(4);

                    // set the points
                    result[p][0] = vtxLabel(0, j, k);
                    result[p][1] = vtxLabel(0, j, k + 1);
                    result[p][2] = vtxLabel(0, j + 1, k + 1);
                    result[p][3] = vtxLabel(0, j + 1, k);

                    p++;
                }
            }

            result.setSize(p);
            break;
        }

        case 5:
        {
            // low j = ymin
            result.setSize
            (
                (xMaxRange - xMinRange + 1)*(zMaxRange - zMinRange + 1)
            );

            label p = 0;
            for (label i = xMinRange - 1; i <= xMaxRange - 1; i++)
            {
                for (label k = zMinRange - 1; k <= zMaxRange - 1; k++)
                {
                    result[p].setSize(4);

                    // set the points
                    result[p][0] = vtxLabel(i, 0, k);
                    result[p][1] = vtxLabel(i + 1, 0, k);
                    result[p][2] = vtxLabel(i + 1, 0, k + 1);
                    result[p][3] = vtxLabel(i, 0, k + 1);

                    p++;
                }
            }

            result.setSize(p);
            break;
        }

        case 6:
        {
            // low k = zmin
            result.setSize
            (
                (xMaxRange - xMinRange + 1)*(yMaxRange - yMinRange + 1)
            );

            label p = 0;
            for (label i = xMinRange - 1; i <= xMaxRange - 1; i++)
            {
                for (label j = yMinRange - 1; j <= yMaxRange - 1; j++)
                {
                    result[p].setSize(4);

                    // set the points
                    result[p][0] = vtxLabel(i, j, 0);
                    result[p][1] = vtxLabel(i, j + 1, 0);
                    result[p][2] = vtxLabel(i + 1, j + 1, 0);
                    result[p][3] = vtxLabel(i + 1, j, 0);

                    p++;
                }
            }

            result.setSize(p);
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "direction out of range (1 to 6): " << direc
                << abort(FatalError);
        }
    }

    // Correct the face orientation based on the handedness of the block.
    // Do nothing for the right-handed block
    if (blockHandedness_ == noPoints)
    {
        FatalErrorInFunction
            << "Unable to determine block handedness as points "
            << "have not been read in yet"
            << abort(FatalError);
    }
    else if (blockHandedness_ == left)
    {
        // turn all faces inside out
        forAll(result, facei)
        {
            result[facei].flip();
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
