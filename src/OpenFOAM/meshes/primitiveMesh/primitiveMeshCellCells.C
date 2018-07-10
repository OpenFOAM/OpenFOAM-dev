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
            FatalErrorInFunction
                << abort(FatalError);
        }
    }

    // It is an error to attempt to recalculate cellCells
    // if the pointer is already set
    if (ccPtr_)
    {
        FatalErrorInFunction
            << "cellCells already calculated"
            << abort(FatalError);
    }
    else
    {
        // 1. Count number of internal faces per cell

        labelList ncc(nCells(), 0);

        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();

        forAll(nei, facei)
        {
            ncc[own[facei]]++;
            ncc[nei[facei]]++;
        }

        // Create the storage
        ccPtr_ = new labelListList(ncc.size());
        labelListList& cellCellAddr = *ccPtr_;



        // 2. Size and fill cellFaceAddr

        forAll(cellCellAddr, celli)
        {
            cellCellAddr[celli].setSize(ncc[celli]);
        }
        ncc = 0;

        forAll(nei, facei)
        {
            label ownCelli = own[facei];
            label neiCelli = nei[facei];

            cellCellAddr[ownCelli][ncc[ownCelli]++] = neiCelli;
            cellCellAddr[neiCelli][ncc[neiCelli]++] = ownCelli;
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
    const label celli,
    DynamicList<label>& storage
) const
{
    if (hasCellCells())
    {
        return cellCells()[celli];
    }
    else
    {
        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();
        const cell& cFaces = cells()[celli];

        storage.clear();

        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            if (facei < nInternalFaces())
            {
                if (own[facei] == celli)
                {
                    storage.append(nei[facei]);
                }
                else
                {
                    storage.append(own[facei]);
                }
            }
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::cellCells(const label celli) const
{
    return cellCells(celli, labels_);
}


// ************************************************************************* //
