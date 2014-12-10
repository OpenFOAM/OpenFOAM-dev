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

Description
    Create intermediate mesh from PROSTAR files

\*---------------------------------------------------------------------------*/

#include "starMesh.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label starMesh::readVtxLabel(IFstream& is)
{
    char lcs[16];

    for (int i=0; i<15; i++)
    {
        is.get(lcs[i]);
    }

    lcs[15] = '\0';

    return atoi(lcs);
}


scalar starMesh::readVtxCmpt(IFstream& is)
{
    char lcs[17];

    for (int i=0; i<16; i++)
    {
        is.get(lcs[i]);
    }

    lcs[16] = '\0';

    return scalar(atof(lcs));
}


void starMesh::readToNl(IFstream& is)
{
    char c;
    do
    {
        is.get(c);
    } while (is && c != '\n');
}


void starMesh::readPoints(const scalar scaleFactor)
{
    label nPoints = 0;
    label maxLabel = -1;

    fileName pointsFileName(casePrefix_ + ".vrt");

    {
        IFstream pointsFile(pointsFileName);

        // Pass 1: get # points and maximum vertex label

        if (pointsFile.good())
        {
            label pointLabel;

            maxLabel = -1;
            while (pointsFile)
            {
                pointLabel = readVtxLabel(pointsFile);

                if (!pointsFile) break;

                maxLabel = max(maxLabel, pointLabel);

                readVtxCmpt(pointsFile);
                readVtxCmpt(pointsFile);
                readVtxCmpt(pointsFile);

                readToNl(pointsFile);

                nPoints++;
            }
        }
        else
        {
            FatalErrorIn("starMesh::readPoints()")
                << "Cannot read file " << pointsFileName
                << abort(FatalError);
        }
    }

    Info<< "Number of points = " << nPoints << endl << endl;

    points_.setSize(nPoints);

#   ifdef starMesh_H
    starPointID_.setSize(nPoints);

    // Reset STAR point ID, just in case
    starPointID_ = -1;
#   endif

    starPointLabelLookup_.setSize(maxLabel+1);

    // reset point labels to invalid value
    starPointLabelLookup_ = -1;

    if (nPoints > 0)
    {
        // Pass 2: construct pointlist and conversion table
        // from Star vertex numbers to Foam pointLabels

        IFstream pointsFile(pointsFileName);
        label pointLabel;

        forAll(points_, p)
        {
            pointLabel = readVtxLabel(pointsFile);
            points_[p].x() = readVtxCmpt(pointsFile);
            points_[p].y() = readVtxCmpt(pointsFile);
            points_[p].z() = readVtxCmpt(pointsFile);

            readToNl(pointsFile);

#           ifdef starMesh_H
            starPointID_[p] = pointLabel;
#           endif

            starPointLabelLookup_[pointLabel] = p;
        }

        if (scaleFactor > 1.0 + SMALL || scaleFactor < 1.0 - SMALL)
        {
            points_ *= scaleFactor;
        }
    }
    else
    {
        FatalError
            << "void starMesh::readPoints() : "
            << "no points in file "
            << pointsFileName
            << abort(FatalError);
    }
}


// ************************************************************************* //
