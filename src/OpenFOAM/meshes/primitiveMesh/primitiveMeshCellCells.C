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

#include "primitiveMesh.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcCellCells() const
{
    // Loop through faceCells and mark up neighbours

    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCells() : calculating cellCells"
            << endl;

        if (debug == -1)
        {
            // For checking calls:abort so we can quickly hunt down
            // origin of call
            FatalErrorIn("primitiveMesh::calcCellCells()")
                << abort(FatalError);
        }
    }

    // It is an error to attempt to recalculate cellCells
    // if the pointer is already set
    if (ccPtr_)
    {
        FatalErrorIn("primitiveMesh::calcCellCells() const")
            << "cellCells already calculated"
            << abort(FatalError);
    }
    else
    {
        // 1. Count number of internal faces per cell

        labelList ncc(nCells(), 0);

        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();

        forAll(nei, faceI)
        {
            ncc[own[faceI]]++;
            ncc[nei[faceI]]++;
        }

        // Create the storage
        ccPtr_ = new labelListList(ncc.size());
        labelListList& cellCellAddr = *ccPtr_;



        // 2. Size and fill cellFaceAddr

        forAll(cellCellAddr, cellI)
        {
            cellCellAddr[cellI].setSize(ncc[cellI]);
        }
        ncc = 0;

        forAll(nei, faceI)
        {
            label ownCellI = own[faceI];
            label neiCellI = nei[faceI];

            cellCellAddr[ownCellI][ncc[ownCellI]++] = neiCellI;
            cellCellAddr[neiCellI][ncc[neiCellI]++] = ownCellI;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList& Foam::primitiveMesh::cellCells() const
{
    if (!ccPtr_)
    {
        calcCellCells();
    }

    return *ccPtr_;
}


const Foam::labelList& Foam::primitiveMesh::cellCells
(
    const label cellI,
    DynamicList<label>& storage
) const
{
    if (hasCellCells())
    {
        return cellCells()[cellI];
    }
    else
    {
        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();
        const cell& cFaces = cells()[cellI];

        storage.clear();

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            if (faceI < nInternalFaces())
            {
                if (own[faceI] == cellI)
                {
                    storage.append(nei[faceI]);
                }
                else
                {
                    storage.append(own[faceI]);
                }
            }
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::cellCells(const label cellI) const
{
    return cellCells(cellI, labels_);
}


// ************************************************************************* //
