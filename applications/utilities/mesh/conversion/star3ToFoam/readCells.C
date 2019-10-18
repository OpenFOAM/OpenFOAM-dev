/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Description
    Create intermediate mesh from Prostar files

\*---------------------------------------------------------------------------*/

#include "starMesh.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::starMesh::addRegularCell
(
    const labelList& labels,
    const label nCreatedCells
)
{
    // Memory management
    static labelList labelsHex(8);
    static labelList labelsPrism(6);
    static labelList labelsPyramid(5);
    static labelList labelsTet(4);
    static labelList labelsTetWedge(5);

    label regularTypeFlag = -1;

    // grab the shape from the table
    const cellModel* curModelPtr = reinterpret_cast<cellModel*>(0);

    if      // Tetrahedron
    (
        labels[2] == labels[3]
     && labels[4] == labels[5]
     && labels[5] == labels[6]
     && labels[6] == labels[7]
    )
    {
        regularTypeFlag = 0;
        curModelPtr = tetPtr_;
    }
    else if // Square-based pyramid
    (
        labels[4] == labels[5]
     && labels[5] == labels[6]
     && labels[6] == labels[7]
    )
    {
        regularTypeFlag = 1;
        curModelPtr = pyrPtr_;
    }
    else if // Tet Wedge
    (
        labels[2] == labels[3]
     && labels[4] == labels[5]
     && labels[6] == labels[7]
    )
    {
        regularTypeFlag = 2;
        curModelPtr = tetWedgePtr_;
    }
    else if // Triangular prism
    (
        labels[2] == labels[3]
     && labels[6] == labels[7]
    )
    {
        regularTypeFlag = 3;
        curModelPtr = prismPtr_;
    }
    else if // Wedge
    (
        labels[4] == labels[7]
    )
    {
        regularTypeFlag = 4;
        curModelPtr = wedgePtr_;
    }
    else    // Hex
    {
        regularTypeFlag = 5;
        curModelPtr = hexPtr_;
    }

    labelList regularCellLabels(curModelPtr->nPoints(), -1);
    // get reference to the addressing list
    const label* addressing = regularAddressingTable[regularTypeFlag];

    forAll(regularCellLabels, labelI)
    {
        regularCellLabels[labelI] = labels[addressing[labelI]];
    }

    cellShapes_[nCreatedCells] = cellShape(*curModelPtr, regularCellLabels);
}


void Foam::starMesh::addSAMMcell
(
    const labelList& labels,
    const label nCreatedCells
)
{
    // get type, reg and permutation flag
    label typeFlag = labels[21];
//     label regularityFlag = labels[22];  // Not used.
    label permutationFlag = labels[23];

    // grab the shape from the table
    label sammTypeFlag = -1;
    const cellModel* curModelPtr = reinterpret_cast<cellModel*>(0);

    switch (typeFlag)
    {
        case 1:
        {
            sammTypeFlag = 1;
            curModelPtr = sammTrim1Ptr_;
            break;
        }

        case 2:
        {
            sammTypeFlag = 2;
            curModelPtr = sammTrim2Ptr_;
            break;
        }

        case 7:
        {
            if (labels[0] != -1)
            {
                sammTypeFlag = 3;
                curModelPtr = sammTrim3Ptr_;
            }
            else
            {
                sammTypeFlag = 5;
                curModelPtr = sammTrim5Ptr_;
            }

            break;
        }

        case 8:
        {
            sammTypeFlag = 4;
            curModelPtr = sammTrim4Ptr_;
            break;
        }

        case 85:
        {
            sammTypeFlag = 8;
            curModelPtr = sammTrim8Ptr_;
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "SAMM type " << sammTypeFlag << " is invalid"
                << abort(FatalError);
        }
    }

    // make a list of labels
    labelList sammCellLabels(curModelPtr->nPoints(), -1);
    // get reference to the addressing list
    const label* addressing = sammAddressingTable[sammTypeFlag];

    forAll(sammCellLabels, labelI)
    {
        sammCellLabels[labelI] = labels[addressing[labelI]];
    }

    cellShapes_[nCreatedCells] = cellShape(*curModelPtr, sammCellLabels);

    // set permutation flag for cell
    starCellPermutation_[nCreatedCells] = permutationFlag;
}


void Foam::starMesh::readCells()
{
    label nCells = 0;
    label maxLabel = -1;

    fileName cellsFileName(casePrefix_ + ".cel");

    {
        IFstream cellsFile(cellsFileName);

        if (cellsFile.good())
        {
            label lineLabel, pointLabel, regionLabel, typeFlag;

            maxLabel = -1;
            while (!(cellsFile >> lineLabel).eof())
            {
                maxLabel = max(maxLabel, lineLabel);
                for (int i=0; i<8; i++)
                {
                    cellsFile >> pointLabel;
                }

                cellsFile >> regionLabel;
                cellsFile >> typeFlag;

                // lines with typeFlag of zero are continuation lines.
                if (typeFlag != 0)
                {
                    nCells++;
                }

                // backward compatibility: number of trailing rubbish in
                // STAR is unknown.
                // Fixed to cope with missing \n on last line.
                readToNl(cellsFile);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Cannot read file " << cellsFileName
                << abort(FatalError);
        }
    }

    Info<< "Number of cells = " << nCells << endl << endl;

    cellShapes_.setSize(nCells);
    starCellID_.setSize(nCells);
    starCellPermutation_.setSize(nCells);

    // reset permutation to invalid value
    forAll(starCellPermutation_, i)
    {
        starCellPermutation_[i] = -1;
    }

    starCellLabelLookup_.setSize(maxLabel+1);

    // reset point labels to invalid value
    forAll(starCellLabelLookup_, i)
    {
        starCellLabelLookup_[i] = -1;
    }

    if (nCells > 0)
    {
        IFstream cellsFile(cellsFileName);

        labelList labels(24, label(-1));
        label lineLabel, starLabel, regionLabel, typeFlag;

        for (label celli = 0; celli < nCells; celli++)
        {
            label nLabels = 0;

            label addOnToCell = 0;

            // reset the labels to -1. Debugging.
            forAll(labels, i)
            {
                labels[i] = -1;
            }

            do
            {
                if ((cellsFile >> lineLabel).eof())
                {
                    FatalErrorInFunction
                        << "Reached end of cells file before "
                        << "all cells are read in."
                        << abort(FatalError);
                }

                nLabels += 8;

                for (int i=nLabels-8; i<nLabels; i++)
                {
                    cellsFile >> starLabel;

                    if (i < 21)
                    {
                        if (starLabel != 0)
                        {
                            // Convert Star vertex number to point label
                            labels[i] = starPointLabelLookup_[starLabel];

                            if (labels[i] < 0)
                            {
                                Info<< "Cells not consistent with vertex file. "
                                    << "Star vertex number " << starLabel
                                    << " does not exist\n";
                            }
                        }
                        else
                        {
                            labels[i] = -1;
                        }
                    }
                    else
                    {
                        labels[i] = starLabel;
                    }
                }

                cellsFile >> regionLabel;
                cellsFile >> typeFlag;

                // check for continuation line
                if (typeFlag == -1)
                {
                    addOnToCell = 2;
                }

                // backward compatibility: number of trailing rubbish in
                // STAR is unknown.
                readToNl(cellsFile);

                addOnToCell--;

            } while (addOnToCell >= 0);

            // Record STAR cell number (used for debugging)
            starCellID_[celli] = lineLabel;

            // insert STAR lookup addressing
            starCellLabelLookup_[lineLabel] = celli;

            if (nLabels == 8)
            {
                addRegularCell(labels, celli);
            }
            else
            {
                addSAMMcell(labels, celli);
            }

            // check cell labels
            const labelList& curShapeLabels = cellShapes_[celli];

            forAll(curShapeLabels, i)
            {
                if (curShapeLabels[i] < 0)
                {
                    FatalErrorInFunction
                        << "Invalid vertex found in cell " << celli
                        << ". STAR cell no: " << lineLabel
                        << " labels: " << curShapeLabels
                        << abort(FatalError);
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "No cells in file " << cellsFileName
            << abort(FatalError);
    }
}


// ************************************************************************* //
