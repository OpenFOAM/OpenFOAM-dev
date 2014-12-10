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

#include "FECCellToFaceStencil.H"
#include "syncTools.H"
#include "emptyPolyPatch.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculates per edge the neighbour data (= edgeCells)
void Foam::FECCellToFaceStencil::calcEdgeBoundaryData
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


// Calculates per face the edge connected data (= cell or boundary in global
// numbering).
void Foam::FECCellToFaceStencil::calcFaceStencil
(
    labelListList& faceStencil
) const
{
    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    const label nBnd = mesh().nFaces()-mesh().nInternalFaces();
    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();



    // Determine neighbouring global cell
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList neiGlobalCell(nBnd);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                neiGlobalCell[faceI-mesh().nInternalFaces()] =
                    globalNumbering().toGlobal(own[faceI]);
                faceI++;
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh(), neiGlobalCell);



    // Determine on coupled edges the edge cells on the other side
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Calculate edges on coupled patches
    labelList boundaryEdges
    (
        allCoupledFacesPatch()().meshEdges
        (
            mesh().edges(),
            mesh().pointEdges()
        )
    );

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



    // Construct stencil in global numbering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    faceStencil.setSize(mesh().nFaces());

    // Do coupled edges first

    forAll(boundaryEdges, i)
    {
        label edgeI = boundaryEdges[i];

        const labelList& eGlobals = neiGlobal[mesh().edges()[edgeI]];

        // Distribute to all edgeFaces
        const labelList& eFaces = mesh().edgeFaces(edgeI);

        forAll(eFaces, j)
        {
            label faceI = eFaces[j];

            // Insert eGlobals into faceStencil.
            merge(-1, -1, eGlobals, faceStencil[faceI]);
        }
    }
    neiGlobal.clear();


    // Do remaining edges by looping over all faces

    // Work arrays
    DynamicList<label> fEdgesSet;
    DynamicList<label> eFacesSet;
    labelHashSet faceStencilSet;

    for (label faceI = 0; faceI < mesh().nInternalFaces(); faceI++)
    {
        label globalOwn = globalNumbering().toGlobal(own[faceI]);
        label globalNei = globalNumbering().toGlobal(nei[faceI]);

        // Convert any existing faceStencil (from coupled edges) into
        // set and operate on this.

        faceStencilSet.clear();

        // Insert all but global owner and neighbour
        forAll(faceStencil[faceI], i)
        {
            label globalI = faceStencil[faceI][i];
            if (globalI != globalOwn && globalI != globalNei)
            {
                faceStencilSet.insert(globalI);
            }
        }
        faceStencil[faceI].clear();

        // Collect all edge connected (internal) cells
        const labelList& fEdges = mesh().faceEdges(faceI, fEdgesSet);

        forAll(fEdges, i)
        {
            label edgeI = fEdges[i];

            insertFaceCells
            (
                globalOwn,
                globalNei,
                isValidBFace,
                mesh().edgeFaces(edgeI, eFacesSet),
                faceStencilSet
            );
        }

        // Extract, guarantee owner first, neighbour second.
        faceStencil[faceI].setSize(faceStencilSet.size()+2);
        label n = 0;
        faceStencil[faceI][n++] = globalOwn;
        faceStencil[faceI][n++] = globalNei;
        forAllConstIter(labelHashSet, faceStencilSet, iter)
        {
            if (iter.key() == globalOwn || iter.key() == globalNei)
            {
                FatalErrorIn("FECCellToFaceStencil::calcFaceStencil(..)")
                    << "problem:" << faceStencilSet
                    << abort(FatalError);
            }
            faceStencil[faceI][n++] = iter.key();
        }
    }
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        label faceI = pp.start();

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label globalOwn = globalNumbering().toGlobal(own[faceI]);
                label globalNei = neiGlobalCell[faceI-mesh().nInternalFaces()];

                // Convert any existing faceStencil (from coupled edges) into
                // set and operate on this.

                faceStencilSet.clear();

                // Insert all but global owner and neighbour
                forAll(faceStencil[faceI], i)
                {
                    label globalI = faceStencil[faceI][i];
                    if (globalI != globalOwn && globalI != globalNei)
                    {
                        faceStencilSet.insert(globalI);
                    }
                }
                faceStencil[faceI].clear();

                // Collect all edge connected (internal) cells
                const labelList& fEdges = mesh().faceEdges(faceI, fEdgesSet);

                forAll(fEdges, i)
                {
                    label edgeI = fEdges[i];

                    insertFaceCells
                    (
                        globalOwn,
                        globalNei,
                        isValidBFace,
                        mesh().edgeFaces(edgeI, eFacesSet),
                        faceStencilSet
                    );
                }

                // Extract, guarantee owner first, neighbour second.
                faceStencil[faceI].setSize(faceStencilSet.size()+2);
                label n = 0;
                faceStencil[faceI][n++] = globalOwn;
                faceStencil[faceI][n++] = globalNei;
                forAllConstIter(labelHashSet, faceStencilSet, iter)
                {
                    if (iter.key() == globalOwn || iter.key() == globalNei)
                    {
                        FatalErrorIn
                        (
                            "FECCellToFaceStencil::calcFaceStencil(..)"
                        )   << "problem:" << faceStencilSet
                            << abort(FatalError);
                    }
                    faceStencil[faceI][n++] = iter.key();
                }

                if (n != faceStencil[faceI].size())
                {
                    FatalErrorIn("problem") << "n:" << n
                        << " size:" << faceStencil[faceI].size()
                        << abort(FatalError);
                }

                faceI++;
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                label globalOwn = globalNumbering().toGlobal(own[faceI]);

                // Convert any existing faceStencil (from coupled edges) into
                // set and operate on this.

                faceStencilSet.clear();

                // Insert all but global owner and neighbour
                forAll(faceStencil[faceI], i)
                {
                    label globalI = faceStencil[faceI][i];
                    if (globalI != globalOwn)
                    {
                        faceStencilSet.insert(globalI);
                    }
                }
                faceStencil[faceI].clear();

                // Collect all edge connected (internal) cells
                const labelList& fEdges = mesh().faceEdges(faceI, fEdgesSet);

                forAll(fEdges, i)
                {
                    label edgeI = fEdges[i];

                    insertFaceCells
                    (
                        globalOwn,
                        -1,
                        isValidBFace,
                        mesh().edgeFaces(edgeI, eFacesSet),
                        faceStencilSet
                    );
                }

                // Extract, guarantee owner first, neighbour second.
                faceStencil[faceI].setSize(faceStencilSet.size()+1);
                label n = 0;
                faceStencil[faceI][n++] = globalOwn;
                forAllConstIter(labelHashSet, faceStencilSet, iter)
                {
                    if (iter.key() == globalOwn)
                    {
                        FatalErrorIn
                        (
                            "FECCellToFaceStencil::calcFaceStencil(..)"
                        )   << "problem:" << faceStencilSet
                            << abort(FatalError);
                    }
                    faceStencil[faceI][n++] = iter.key();
                }

                faceI++;
            }
        }
    }


    for (label faceI = 0; faceI < mesh().nInternalFaces(); faceI++)
    {
        label globalOwn = globalNumbering().toGlobal(own[faceI]);
        if (faceStencil[faceI][0] != globalOwn)
        {
            FatalErrorIn("FECCellToFaceStencil::calcFaceStencil(..)")
                << "problem:" << faceStencil[faceI]
                << " globalOwn:" << globalOwn
                << abort(FatalError);
        }
        label globalNei = globalNumbering().toGlobal(nei[faceI]);
        if (faceStencil[faceI][1] != globalNei)
        {
            FatalErrorIn("FECCellToFaceStencil::calcFaceStencil(..)")
                << "problem:" << faceStencil[faceI]
                << " globalNei:" << globalNei
                << abort(FatalError);
        }
    }


    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;

                label globalOwn = globalNumbering().toGlobal(own[faceI]);
                if (faceStencil[faceI][0] != globalOwn)
                {
                    FatalErrorIn("FECCellToFaceStencil::calcFaceStencil(..)")
                        << "problem:" << faceStencil[faceI]
                        << " globalOwn:" << globalOwn
                        << abort(FatalError);
                }
                label globalNei = neiGlobalCell[faceI-mesh().nInternalFaces()];
                if (faceStencil[faceI][1] != globalNei)
                {
                    FatalErrorIn("FECCellToFaceStencil::calcFaceStencil(..)")
                        << "problem:" << faceStencil[faceI]
                        << " globalNei:" << globalNei
                        << abort(FatalError);
                }
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;

                label globalOwn = globalNumbering().toGlobal(own[faceI]);
                if (faceStencil[faceI][0] != globalOwn)
                {
                    FatalErrorIn("FECCellToFaceStencil::calcFaceStencil(..)")
                        << "problem:" << faceStencil[faceI]
                        << " globalOwn:" << globalOwn
                        << abort(FatalError);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FECCellToFaceStencil::FECCellToFaceStencil(const polyMesh& mesh)
:
    cellToFaceStencil(mesh)
{
    // Calculate per face the (edge) connected cells (in global numbering)
    labelListList faceStencil;
    calcFaceStencil(faceStencil);

    // Transfer to *this
    transfer(faceStencil);
}


// ************************************************************************* //
