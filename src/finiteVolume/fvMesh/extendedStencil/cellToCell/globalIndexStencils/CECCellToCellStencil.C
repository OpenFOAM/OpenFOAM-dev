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

#include "CECCellToCellStencil.H"
#include "syncTools.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculates per edge the neighbour data (= edgeCells)
void Foam::CECCellToCellStencil::calcEdgeBoundaryData
(
    const boolList& isValidBFace,
    const labelList& boundaryEdges,
    EdgeMap<labelList>& neiGlobal
) const
{
    neiGlobal.resize(2*boundaryEdges.size());

    labelHashSet edgeGlobals;

    forAll(boundaryEdges, i)
    {
        label edgeI = boundaryEdges[i];

        neiGlobal.insert
        (
            mesh().edges()[edgeI],
            calcFaceCells
            (
                isValidBFace,
                mesh().edgeFaces(edgeI),
                edgeGlobals
            )
        );
    }

    syncTools::syncEdgeMap(mesh(), neiGlobal, unionEqOp(), dummyTransform());
}


// Calculates per cell the neighbour data (= cell or boundary in global
// numbering). First element is always cell itself!
void Foam::CECCellToCellStencil::calcCellStencil
(
    labelListList& globalCellCells
) const
{
    // Calculate edges on coupled patches
    labelList boundaryEdges
    (
        allCoupledFacesPatch()().meshEdges
        (
            mesh().edges(),
            mesh().pointEdges()
        )
    );

    //{
    //    OFstream str(mesh().time().path()/"boundaryEdges.obj");
    //    Pout<< "DUmping boundary edges to " << str.name() << endl;
    //
    //    label vertI = 0;
    //    forAll(boundaryEdges, i)
    //    {
    //        label edgeI = boundaryEdges[i];
    //        const edge& e = mesh().edges()[edgeI];
    //        const point& p0 = mesh().points()[e[0]];
    //        const point& p1 = mesh().points()[e[1]];
    //
    //        Pout<< "boundary edge " << edgeI << " between " << p0 << p1
    //            << endl;
    //
    //        meshTools::writeOBJ(str, p0);
    //        vertI++;
    //        meshTools::writeOBJ(str, p1);
    //        vertI++;
    //        str << "l " << vertI-1 << ' ' << vertI << nl;
    //    }
    //}


    // Mark boundary faces to be included in stencil (i.e. not coupled or empty)
    boolList isValidBFace;
    validBoundaryFaces(isValidBFace);


    // Swap edgeCells for coupled edges. Note: use EdgeMap for now since we've
    // got syncTools::syncEdgeMap for those. Should be replaced with Map and
    // syncTools functionality to handle those.
    EdgeMap<labelList> neiGlobal;
    calcEdgeBoundaryData
    (
        isValidBFace,
        boundaryEdges,
        neiGlobal
    );

    globalCellCells.setSize(mesh().nCells());

    // Do coupled edges first

    forAll(boundaryEdges, i)
    {
        label edgeI = boundaryEdges[i];

        const labelList& eGlobals = neiGlobal[mesh().edges()[edgeI]];

        // Distribute to all edgeCells
        const labelList& eCells = mesh().edgeCells(edgeI);

        forAll(eCells, j)
        {
            label cellI = eCells[j];

            // Insert pGlobals into globalCellCells
            merge
            (
                globalNumbering().toGlobal(cellI),
                eGlobals,
                globalCellCells[cellI]
            );
        }
    }
    neiGlobal.clear();

    // Do remaining edges cells
    labelHashSet edgeGlobals;

    for (label edgeI = 0; edgeI < mesh().nEdges(); edgeI++)
    {
        labelList eGlobals
        (
            calcFaceCells
            (
                isValidBFace,
                mesh().edgeFaces(edgeI),
                edgeGlobals
            )
        );

        const labelList& eCells = mesh().edgeCells(edgeI);

        forAll(eCells, j)
        {
            label cellI = eCells[j];

            merge
            (
                globalNumbering().toGlobal(cellI),
                eGlobals,
                globalCellCells[cellI]
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CECCellToCellStencil::CECCellToCellStencil(const polyMesh& mesh)
:
    cellToCellStencil(mesh)
{
    // Calculate per cell the (edge) connected cells (in global numbering)
    calcCellStencil(*this);
}


// ************************************************************************* //
