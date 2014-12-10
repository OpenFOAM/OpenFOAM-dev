/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "CPCCellToCellStencil.H"
#include "syncTools.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculates per point the neighbour data (= pointCells)
void Foam::CPCCellToCellStencil::calcPointBoundaryData
(
    const boolList& isValidBFace,
    const labelList& boundaryPoints,
    Map<labelList>& neiGlobal
) const
{
    neiGlobal.resize(2*boundaryPoints.size());

    labelHashSet pointGlobals;

    forAll(boundaryPoints, i)
    {
        label pointI = boundaryPoints[i];

        neiGlobal.insert
        (
            pointI,
            calcFaceCells
            (
                isValidBFace,
                mesh().pointFaces()[pointI],
                pointGlobals
            )
        );
    }

    syncTools::syncPointMap
    (
        mesh(),
        neiGlobal,
        unionEqOp(),
        Foam::dummyTransform()      // dummy transformation
    );
}


// Calculates per cell the neighbour data (= cell or boundary in global
// numbering). First element is always cell itself!
void Foam::CPCCellToCellStencil::calcCellStencil
(
    labelListList& globalCellCells
) const
{
    // Calculate points on coupled patches
    labelList boundaryPoints(allCoupledFacesPatch()().meshPoints());


    // Mark boundary faces to be included in stencil (i.e. not coupled or empty)
    boolList isValidBFace;
    validBoundaryFaces(isValidBFace);


    // Swap pointCells for coupled points
    Map<labelList> neiGlobal;
    calcPointBoundaryData
    (
        isValidBFace,
        boundaryPoints,
        neiGlobal
    );

    globalCellCells.setSize(mesh().nCells());

    // Do coupled points first

    forAll(boundaryPoints, i)
    {
        label pointI = boundaryPoints[i];

        const labelList& pGlobals = neiGlobal[pointI];

        // Distribute to all pointCells
        const labelList& pCells = mesh().pointCells(pointI);

        forAll(pCells, j)
        {
            label cellI = pCells[j];

            // Insert pGlobals into globalCellCells
            merge
            (
                globalNumbering().toGlobal(cellI),
                pGlobals,
                globalCellCells[cellI]
            );
        }
    }
    neiGlobal.clear();

    // Do remaining points cells
    labelHashSet pointGlobals;

    for (label pointI = 0; pointI < mesh().nPoints(); pointI++)
    {
        labelList pGlobals
        (
            calcFaceCells
            (
                isValidBFace,
                mesh().pointFaces()[pointI],
                pointGlobals
            )
        );

        const labelList& pCells = mesh().pointCells(pointI);

        forAll(pCells, j)
        {
            label cellI = pCells[j];

            merge
            (
                globalNumbering().toGlobal(cellI),
                pGlobals,
                globalCellCells[cellI]
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CPCCellToCellStencil::CPCCellToCellStencil(const polyMesh& mesh)
:
    cellToCellStencil(mesh)
{
    // Calculate per cell the (point) connected cells (in global numbering)
    labelListList globalCellCells;
    calcCellStencil(*this);
}


// ************************************************************************* //
