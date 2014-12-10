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

#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMesh::initMesh()
{
    if (debug)
    {
        Info<< "void polyMesh::initMesh() : "
            << "initialising primitiveMesh" << endl;
    }

    // For backward compatibility check if the neighbour array is the same
    // length as the owner and shrink to remove the -1s padding
    if (neighbour_.size() == owner_.size())
    {
        label nInternalFaces = 0;

        forAll(neighbour_, faceI)
        {
            if (neighbour_[faceI] == -1)
            {
                break;
            }
            else
            {
                nInternalFaces++;
            }
        }

        neighbour_.setSize(nInternalFaces);
    }

    label nCells = -1;

    forAll(owner_, facei)
    {
        if (owner_[facei] < 0)
        {
            FatalErrorIn("polyMesh::initMesh()")
                << "Illegal cell label " << owner_[facei]
                << " in neighbour addressing for face " << facei
                << exit(FatalError);
        }
        nCells = max(nCells, owner_[facei]);
    }

    // The neighbour array may or may not be the same length as the owner
    forAll(neighbour_, facei)
    {
        if (neighbour_[facei] < 0)
        {
            FatalErrorIn("polyMesh::initMesh()")
                << "Illegal cell label " << neighbour_[facei]
                << " in neighbour addressing for face " << facei
                << exit(FatalError);
        }
        nCells = max(nCells, neighbour_[facei]);
    }

    nCells++;

    // Reset the primitiveMesh with the sizes of the primitive arrays
    primitiveMesh::reset
    (
        points_.size(),
        neighbour_.size(),
        owner_.size(),
        nCells
    );

    string meshInfo =
        "nPoints:" + Foam::name(nPoints())
      + "  nCells:" + Foam::name(this->nCells())
      + "  nFaces:" + Foam::name(nFaces())
      + "  nInternalFaces:" + Foam::name(nInternalFaces());

    owner_.note() = meshInfo;
    neighbour_.note() = meshInfo;
}


void Foam::polyMesh::initMesh(cellList& c)
{
    if (debug)
    {
        Info<< "void polyMesh::initMesh(cellList& c) : "
            << "calculating owner-neighbour arrays" << endl;
    }

    owner_.setSize(faces_.size(), -1);
    neighbour_.setSize(faces_.size(), -1);

    boolList markedFaces(faces_.size(), false);

    label nInternalFaces = 0;

    forAll(c, cellI)
    {
        // get reference to face labels for current cell
        const labelList& cellfaces = c[cellI];

        forAll(cellfaces, faceI)
        {
            if (cellfaces[faceI] < 0)
            {
                FatalErrorIn("polyMesh::initMesh(cellList&)")
                    << "Illegal face label " << cellfaces[faceI]
                    << " in cell " << cellI
                    << exit(FatalError);
            }

            if (!markedFaces[cellfaces[faceI]])
            {
                // First visit: owner
                owner_[cellfaces[faceI]] = cellI;
                markedFaces[cellfaces[faceI]] = true;
            }
            else
            {
                // Second visit: neighbour
                neighbour_[cellfaces[faceI]] = cellI;
                nInternalFaces++;
            }
        }
    }

    // The neighbour array is initialised with the same length as the owner
    // padded with -1s and here it is truncated to the correct size of the
    // number of internal faces.
    neighbour_.setSize(nInternalFaces);

    // Reset the primitiveMesh
    primitiveMesh::reset
    (
        points_.size(),
        neighbour_.size(),
        owner_.size(),
        c.size(),
        c
    );

    string meshInfo =
        "nPoints: " + Foam::name(nPoints())
      + " nCells: " + Foam::name(nCells())
      + " nFaces: " + Foam::name(nFaces())
      + " nInternalFaces: " + Foam::name(this->nInternalFaces());

    owner_.note() = meshInfo;
    neighbour_.note() = meshInfo;
}

// ************************************************************************* //
