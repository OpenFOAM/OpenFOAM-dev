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

Description
    Agglomerate one level using the MGridGen algorithm.

\*---------------------------------------------------------------------------*/

#include "MGridGenGAMGAgglomeration.H"
#include "fvMesh.H"
#include "syncTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::MGridGenGAMGAgglomeration::
makeCompactCellFaceAddressingAndFaceWeights
(
    const lduAddressing& fineAddressing,
    List<idxtype>& cellCells,
    List<idxtype>& cellCellOffsets,
    const scalarField& magSi,
    List<scalar>& faceWeights
)
{
    const label nFineCells = fineAddressing.size();
    const label nFineFaces = fineAddressing.upperAddr().size();

    const labelUList& upperAddr = fineAddressing.upperAddr();
    const labelUList& lowerAddr = fineAddressing.lowerAddr();

    // Number of neighbours for each cell
    labelList nNbrs(nFineCells, 0);

    forAll(upperAddr, facei)
    {
        nNbrs[upperAddr[facei]]++;
    }

    forAll(lowerAddr, facei)
    {
        nNbrs[lowerAddr[facei]]++;
    }

    // Set the sizes of the addressing and faceWeights arrays
    cellCellOffsets.setSize(nFineCells + 1);
    cellCells.setSize(2*nFineFaces);
    faceWeights.setSize(2*nFineFaces);


    cellCellOffsets[0] = 0;
    forAll(nNbrs, celli)
    {
        cellCellOffsets[celli+1] = cellCellOffsets[celli] + nNbrs[celli];
    }

    // reset the whole list to use as counter
    nNbrs = 0;

    forAll(upperAddr, facei)
    {
        label own = upperAddr[facei];
        label nei = lowerAddr[facei];

        label l1 = cellCellOffsets[own] + nNbrs[own]++;
        label l2 = cellCellOffsets[nei] + nNbrs[nei]++;

        cellCells[l1] = nei;
        cellCells[l2] = own;

        faceWeights[l1] = magSi[facei];
        faceWeights[l2] = magSi[facei];
    }
}


Foam::tmp<Foam::labelField> Foam::MGridGenGAMGAgglomeration::agglomerate
(
    label& nCoarseCells,
    const label minSize,
    const label maxSize,
    const lduAddressing& fineAddressing,
    const scalarField& V,
    const scalarField& magSf,
    const scalarField& magSb
)
{
    const label nFineCells = fineAddressing.size();

    // Compact addressing for cellCells
    List<idxtype> cellCells;
    List<idxtype> cellCellOffsets;

    // Face weights = face areas of the internal faces
    List<scalar> faceWeights;

    // Create the compact addressing for cellCells and faceWeights
    makeCompactCellFaceAddressingAndFaceWeights
    (
        fineAddressing,
        cellCells,
        cellCellOffsets,
        magSf,
        faceWeights
    );

    // agglomeration options.
    List<int> options(4, 0);
    options[0] = 4;                     // globular agglom
    options[1] = 6;                     // objective F3 and F2
    options[2] = 128;                   // debugging output level
    options[3] = fvMesh_.nGeometricD(); // Dimensionality of the grid


    // output: cell -> processor addressing
    List<int> finalAgglom(nFineCells);
    int nMoves = -1;

    MGridGen
    (
        nFineCells,
        cellCellOffsets.begin(),
        const_cast<scalar*>(V.begin()),
        const_cast<scalar*>(magSb.begin()),
        cellCells.begin(),
        faceWeights.begin(),
        minSize,
        maxSize,
        options.begin(),
        &nMoves,
        &nCoarseCells,
        finalAgglom.begin()
    );

    {
        label nNewCoarseCells = 0;
        labelList newRestrictAddr;
        bool ok = checkRestriction
        (
            newRestrictAddr,
            nNewCoarseCells,
            fineAddressing,
            finalAgglom,
            nCoarseCells
        );

        if (!ok)
        {
            nCoarseCells = nNewCoarseCells;
            finalAgglom.transfer(newRestrictAddr);
        }
    }

    return tmp<labelField>(new labelField(finalAgglom));
}


// ************************************************************************* //
