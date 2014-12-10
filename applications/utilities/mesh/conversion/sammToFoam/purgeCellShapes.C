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
    Purge cell shapes which have been rendered invalid by cell face collapse

\*---------------------------------------------------------------------------*/

#include "sammMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void sammMesh::purgeCellShapes()
{
    forAll(cellFaces_, cellI)
    {
        const faceList& curFaces = cellFaces_[cellI];

        // Get model faces
        faceList shapeFaces = cellShapes_[cellI].faces();

        forAll(shapeFaces, faceI)
        {
            bool found = false;

            forAll(curFaces, i)
            {
                if (shapeFaces[faceI] == curFaces[i])
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                Info<< "Purging cell shape " << cellI << endl;
                cellShapes_[cellI] = cellShape(*unknownPtr_, labelList(0));
                break;
            }
        }
    }
}


// ************************************************************************* //
