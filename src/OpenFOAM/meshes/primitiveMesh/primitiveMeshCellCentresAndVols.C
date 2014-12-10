/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Efficient cell-centre calculation using face-addressing, face-centres and
    face-areas.

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcCellCentresAndVols() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Calculating cell centres and cell volumes"
            << endl;
    }

    // It is an error to attempt to recalculate cellCentres
    // if the pointer is already set
    if (cellCentresPtr_ || cellVolumesPtr_)
    {
        FatalErrorIn("primitiveMesh::calcCellCentresAndVols() const")
            << "Cell centres or cell volumes already calculated"
            << abort(FatalError);
    }

    // set the accumulated cell centre to zero vector
    cellCentresPtr_ = new vectorField(nCells());
    vectorField& cellCtrs = *cellCentresPtr_;

    // Initialise cell volumes to 0
    cellVolumesPtr_ = new scalarField(nCells());
    scalarField& cellVols = *cellVolumesPtr_;

    // Make centres and volumes
    makeCellCentresAndVols(faceCentres(), faceAreas(), cellCtrs, cellVols);

    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Finished calculating cell centres and cell volumes"
            << endl;
    }
}


void Foam::primitiveMesh::makeCellCentresAndVols
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cellCtrs,
    scalarField& cellVols
) const
{
    // Clear the fields for accumulation
    cellCtrs = vector::zero;
    cellVols = 0.0;

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    // first estimate the approximate cell centre as the average of
    // face centres

    vectorField cEst(nCells(), vector::zero);
    labelField nCellFaces(nCells(), 0);

    forAll(own, facei)
    {
        cEst[own[facei]] += fCtrs[facei];
        nCellFaces[own[facei]] += 1;
    }

    forAll(nei, facei)
    {
        cEst[nei[facei]] += fCtrs[facei];
        nCellFaces[nei[facei]] += 1;
    }

    forAll(cEst, celli)
    {
        cEst[celli] /= nCellFaces[celli];
    }

    forAll(own, facei)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol =
            fAreas[facei] & (fCtrs[facei] - cEst[own[facei]]);

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*fCtrs[facei] + (1.0/4.0)*cEst[own[facei]];

        // Accumulate volume-weighted face-pyramid centre
        cellCtrs[own[facei]] += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        cellVols[own[facei]] += pyr3Vol;
    }

    forAll(nei, facei)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol =
            fAreas[facei] & (cEst[nei[facei]] - fCtrs[facei]);

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*fCtrs[facei] + (1.0/4.0)*cEst[nei[facei]];

        // Accumulate volume-weighted face-pyramid centre
        cellCtrs[nei[facei]] += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        cellVols[nei[facei]] += pyr3Vol;
    }

    forAll(cellCtrs, celli)
    {
        if (mag(cellVols[celli]) > VSMALL)
        {
            cellCtrs[celli] /= cellVols[celli];
        }
        else
        {
            cellCtrs[celli] = cEst[celli];
        }
    }

    cellVols *= (1.0/3.0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField& Foam::primitiveMesh::cellCentres() const
{
    if (!cellCentresPtr_)
    {
        calcCellCentresAndVols();
    }

    return *cellCentresPtr_;
}


const Foam::scalarField& Foam::primitiveMesh::cellVolumes() const
{
    if (!cellVolumesPtr_)
    {
        calcCellCentresAndVols();
    }

    return *cellVolumesPtr_;
}


// ************************************************************************* //
