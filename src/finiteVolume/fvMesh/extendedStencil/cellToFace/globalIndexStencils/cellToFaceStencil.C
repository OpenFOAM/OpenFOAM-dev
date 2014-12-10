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

#include "cellToFaceStencil.H"
#include "SortableList.H"
#include "emptyPolyPatch.H"
#include "syncTools.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Merge two list and guarantee global0,global1 are first.
void Foam::cellToFaceStencil::merge
(
    const label global0,
    const label global1,
    const labelList& listA,
    labelList& listB
)
{
    sort(listB);

    // See if global0, global1 already present in listB
    label nGlobalInsert = 0;

    if (global0 != -1)
    {
        label index0 = findSortedIndex(listB, global0);
        if (index0 == -1)
        {
            nGlobalInsert++;
        }
    }

    if (global1 != -1)
    {
        label index1 = findSortedIndex(listB, global1);
        if (index1 == -1)
        {
            nGlobalInsert++;
        }
    }


    // For all in listA see if they are present
    label nInsert = 0;

    forAll(listA, i)
    {
        label elem = listA[i];

        if (elem != global0 && elem != global1)
        {
            if (findSortedIndex(listB, elem) == -1)
            {
                nInsert++;
            }
        }
    }

    // Extend B with nInsert and whether global0,global1 need to be inserted.
    labelList result(listB.size() + nGlobalInsert + nInsert);

    label resultI = 0;

    // Insert global0,1 first
    if (global0 != -1)
    {
        result[resultI++] = global0;
    }
    if (global1 != -1)
    {
        result[resultI++] = global1;
    }


    // Insert listB
    forAll(listB, i)
    {
        label elem = listB[i];

        if (elem != global0 && elem != global1)
        {
            result[resultI++] = elem;
        }
    }


    // Insert listA
    forAll(listA, i)
    {
        label elem = listA[i];

        if (elem != global0 && elem != global1)
        {
            if (findSortedIndex(listB, elem) == -1)
            {
                result[resultI++] = elem;
            }
        }
    }

    if (resultI != result.size())
    {
        FatalErrorIn("cellToFaceStencil::merge(..)")
            << "problem" << abort(FatalError);
    }

    listB.transfer(result);
}


// Merge two list and guarantee globalI is first.
void Foam::cellToFaceStencil::merge
(
    const label globalI,
    const labelList& pGlobals,
    labelList& cCells
)
{
    labelHashSet set;
    forAll(cCells, i)
    {
        if (cCells[i] != globalI)
        {
            set.insert(cCells[i]);
        }
    }

    forAll(pGlobals, i)
    {
        if (pGlobals[i] != globalI)
        {
            set.insert(pGlobals[i]);
        }
    }

    cCells.setSize(set.size()+1);
    label n = 0;
    cCells[n++] = globalI;

    forAllConstIter(labelHashSet, set, iter)
    {
        cCells[n++] = iter.key();
    }
}


void Foam::cellToFaceStencil::validBoundaryFaces(boolList& isValidBFace) const
{
    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    isValidBFace.setSize(mesh().nFaces()-mesh().nInternalFaces(), true);

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled() || isA<emptyPolyPatch>(pp))
        {
            label bFaceI = pp.start()-mesh().nInternalFaces();
            forAll(pp, i)
            {
                isValidBFace[bFaceI++] = false;
            }
        }
    }
}


Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::cellToFaceStencil::allCoupledFacesPatch() const
{
    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    label nCoupled = 0;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            nCoupled += pp.size();
        }
    }
    labelList coupledFaces(nCoupled);
    nCoupled = 0;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                coupledFaces[nCoupled++] = faceI++;
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>
            (
                mesh().faces(),
                coupledFaces
            ),
            mesh().points()
        )
    );
}


void Foam::cellToFaceStencil::unionEqOp::operator()
(
    labelList& x,
    const labelList& y
) const
{
    if (y.size())
    {
        if (x.empty())
        {
            x = y;
        }
        else
        {
            labelHashSet set(x);
            forAll(y, i)
            {
                set.insert(y[i]);
            }
            x = set.toc();
        }
    }
}


void Foam::cellToFaceStencil::insertFaceCells
(
    const label exclude0,
    const label exclude1,
    const boolList& isValidBFace,
    const labelList& faceLabels,
    labelHashSet& globals
) const
{
    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();

    forAll(faceLabels, i)
    {
        label faceI = faceLabels[i];

        label globalOwn = globalNumbering().toGlobal(own[faceI]);
        if (globalOwn != exclude0 && globalOwn != exclude1)
        {
            globals.insert(globalOwn);
        }

        if (mesh().isInternalFace(faceI))
        {
            label globalNei = globalNumbering().toGlobal(nei[faceI]);
            if (globalNei != exclude0 && globalNei != exclude1)
            {
                globals.insert(globalNei);
            }
        }
        else
        {
            label bFaceI = faceI-mesh().nInternalFaces();

            if (isValidBFace[bFaceI])
            {
                label globalI = globalNumbering().toGlobal
                (
                    mesh().nCells()
                  + bFaceI
                );

                if (globalI != exclude0 && globalI != exclude1)
                {
                    globals.insert(globalI);
                }
            }
        }
    }
}


Foam::labelList Foam::cellToFaceStencil::calcFaceCells
(
    const boolList& isValidBFace,
    const labelList& faceLabels,
    labelHashSet& globals
) const
{
    globals.clear();

    insertFaceCells
    (
        -1,
        -1,
        isValidBFace,
        faceLabels,
        globals
    );

    return globals.toc();
}


// Calculates per face a list of global cell/face indices.
void Foam::cellToFaceStencil::calcFaceStencil
(
    const labelListList& globalCellCells,
    labelListList& faceStencil
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const label nBnd = mesh_.nFaces()-mesh_.nInternalFaces();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();


    // Determine neighbouring global cell Cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList neiGlobalCellCells(nBnd);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                neiGlobalCellCells[faceI-mesh_.nInternalFaces()] =
                    globalCellCells[own[faceI]];
                faceI++;
            }
        }
    }
    //syncTools::swapBoundaryFaceList(mesh_, neiGlobalCellCells);
    syncTools::syncBoundaryFaceList
    (
        mesh_,
        neiGlobalCellCells,
        eqOp<labelList>(),
        dummyTransform()
    );



    // Construct stencil in global numbering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    faceStencil.setSize(mesh_.nFaces());

    labelHashSet faceStencilSet;

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        faceStencilSet.clear();

        const labelList& ownCCells = globalCellCells[own[faceI]];
        label globalOwn = ownCCells[0];
        // Insert cellCells
        forAll(ownCCells, i)
        {
            faceStencilSet.insert(ownCCells[i]);
        }

        const labelList& neiCCells = globalCellCells[nei[faceI]];
        label globalNei = neiCCells[0];
        // Insert cellCells
        forAll(neiCCells, i)
        {
            faceStencilSet.insert(neiCCells[i]);
        }

        // Guarantee owner first, neighbour second.
        faceStencil[faceI].setSize(faceStencilSet.size());
        label n = 0;
        faceStencil[faceI][n++] = globalOwn;
        faceStencil[faceI][n++] = globalNei;
        forAllConstIter(labelHashSet, faceStencilSet, iter)
        {
            if (iter.key() != globalOwn && iter.key() != globalNei)
            {
                faceStencil[faceI][n++] = iter.key();
            }
        }
        //Pout<< "internalface:" << faceI << " toc:" << faceStencilSet.toc()
        //    << " faceStencil:" << faceStencil[faceI] << endl;
    }
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        label faceI = pp.start();

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                faceStencilSet.clear();

                const labelList& ownCCells = globalCellCells[own[faceI]];
                label globalOwn = ownCCells[0];
                forAll(ownCCells, i)
                {
                    faceStencilSet.insert(ownCCells[i]);
                }

                // And the neighbours of the coupled cell
                const labelList& neiCCells =
                    neiGlobalCellCells[faceI-mesh_.nInternalFaces()];
                label globalNei = neiCCells[0];
                forAll(neiCCells, i)
                {
                    faceStencilSet.insert(neiCCells[i]);
                }

                // Guarantee owner first, neighbour second.
                faceStencil[faceI].setSize(faceStencilSet.size());
                label n = 0;
                faceStencil[faceI][n++] = globalOwn;
                faceStencil[faceI][n++] = globalNei;
                forAllConstIter(labelHashSet, faceStencilSet, iter)
                {
                    if (iter.key() != globalOwn && iter.key() != globalNei)
                    {
                        faceStencil[faceI][n++] = iter.key();
                    }
                }

                //Pout<< "coupledface:" << faceI
                //    << " toc:" << faceStencilSet.toc()
                //    << " faceStencil:" << faceStencil[faceI] << endl;

                faceI++;
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                faceStencilSet.clear();

                const labelList& ownCCells = globalCellCells[own[faceI]];
                label globalOwn = ownCCells[0];
                forAll(ownCCells, i)
                {
                    faceStencilSet.insert(ownCCells[i]);
                }

                // Guarantee owner first
                faceStencil[faceI].setSize(faceStencilSet.size());
                label n = 0;
                faceStencil[faceI][n++] = globalOwn;
                forAllConstIter(labelHashSet, faceStencilSet, iter)
                {
                    if (iter.key() != globalOwn)
                    {
                        faceStencil[faceI][n++] = iter.key();
                    }
                }

                //Pout<< "boundaryface:" << faceI
                //    << " toc:" << faceStencilSet.toc()
                //    << " faceStencil:" << faceStencil[faceI] << endl;

                faceI++;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellToFaceStencil::cellToFaceStencil(const polyMesh& mesh)
:
    mesh_(mesh),
    globalNumbering_(mesh_.nCells()+mesh_.nFaces()-mesh_.nInternalFaces())
{}


// ************************************************************************* //
