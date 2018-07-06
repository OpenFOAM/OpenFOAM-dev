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

#include "undoableMeshCutter.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "DynamicList.H"
#include "meshCutter.H"
#include "cellCuts.H"
#include "splitCell.H"
#include "mapPolyMesh.H"
#include "unitConversion.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(undoableMeshCutter, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// For debugging
void Foam::undoableMeshCutter::printCellRefTree
(
    Ostream& os,
    const word& indent,
    const splitCell* splitCellPtr
) const
{
    if (splitCellPtr)
    {
        os << indent << splitCellPtr->cellLabel() << endl;

        word subIndent = indent + "--";

        printCellRefTree(os, subIndent, splitCellPtr->master());

        printCellRefTree(os, subIndent, splitCellPtr->slave());
    }
}


// For debugging
void Foam::undoableMeshCutter::printRefTree(Ostream& os) const
{
    forAllConstIter(Map<splitCell*>, liveSplitCells_, iter)
    {
        const splitCell* splitPtr = iter();

        // Walk to top (master path only)
        while (splitPtr->parent())
        {
            if (!splitPtr->isMaster())
            {
                splitPtr = nullptr;

                break;
            }
            else
            {
                splitPtr = splitPtr->parent();
            }
        }

        // If we have reached top along master path start printing.
        if (splitPtr)
        {
            // Print from top down
            printCellRefTree(os, word(""), splitPtr);
        }
    }
}


// Update all (cell) labels on splitCell structure.
// Do in two passes to prevent allocation if nothing changed.
void Foam::undoableMeshCutter::updateLabels
(
    const labelList& map,
    Map<splitCell*>& liveSplitCells
)
{
    // Pass1 : check if changed

    bool changed = false;

    forAllConstIter(Map<splitCell*>, liveSplitCells, iter)
    {
        const splitCell* splitPtr = iter();

        if (!splitPtr)
        {
            FatalErrorInFunction
                << "Problem: null pointer on liveSplitCells list"
                << abort(FatalError);
        }

        label celli = splitPtr->cellLabel();

        if (celli != map[celli])
        {
            changed = true;

            break;
        }
    }


    // Pass2: relabel

    if (changed)
    {
        // Build new liveSplitCells
        // since new labels (= keys in Map) might clash with existing ones.
        Map<splitCell*> newLiveSplitCells(2*liveSplitCells.size());

        forAllIter(Map<splitCell*>, liveSplitCells, iter)
        {
            splitCell* splitPtr = iter();

            label celli = splitPtr->cellLabel();

            label newCelli = map[celli];

            if (debug && (celli != newCelli))
            {
                Pout<< "undoableMeshCutter::updateLabels :"
                    << " Updating live (split)cell from " << celli
                    << " to " << newCelli << endl;
            }

            if (newCelli >= 0)
            {
                // Update splitCell. Can do inplace since only one celli will
                // refer to this structure.
                splitPtr->cellLabel() = newCelli;

                // Update liveSplitCells
                newLiveSplitCells.insert(newCelli, splitPtr);
            }
        }
        liveSplitCells = newLiveSplitCells;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::undoableMeshCutter::undoableMeshCutter
(
    const polyMesh& mesh,
    const bool undoable
)
:
    meshCutter(mesh),
    undoable_(undoable),
    liveSplitCells_(mesh.nCells()/100 + 100),
    faceRemover_
    (
        mesh,
        Foam::cos(degToRad(30.0))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::undoableMeshCutter::~undoableMeshCutter()
{
    // Clean split cell tree.

    forAllIter(Map<splitCell*>, liveSplitCells_, iter)
    {
        splitCell* splitPtr = iter();

        while (splitPtr)
        {
            splitCell* parentPtr = splitPtr->parent();

            // Sever ties with parent. Also of other side of refinement since
            // we are handling rest of tree so other side will not have to.
            if (parentPtr)
            {
                splitCell* otherSidePtr = splitPtr->getOther();

                otherSidePtr->parent() = nullptr;

                splitPtr->parent() = nullptr;
            }

            // Delete splitCell (updates pointer on parent to itself)
            delete splitPtr;

            splitPtr = parentPtr;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::undoableMeshCutter::setRefinement
(
    const cellCuts& cuts,
    polyTopoChange& meshMod
)
{
    // Insert commands to actually cut cells
    meshCutter::setRefinement(cuts, meshMod);

    if (undoable_)
    {
        // Use cells cut in this iteration to update splitCell tree.
        forAllConstIter(Map<label>, addedCells(), iter)
        {
            label celli = iter.key();

            label addedCelli = iter();


            // Newly created split cell. (celli ->  celli + addedCelli)

            // Check if celli already part of split.
            Map<splitCell*>::iterator findCell =
                liveSplitCells_.find(celli);

            if (findCell == liveSplitCells_.end())
            {
                // Celli not yet split. It cannot be unlive split cell
                // since that would be illegal to split in the first
                // place.

                // Create 0th level. Null parent to denote this.
                splitCell* parentPtr = new splitCell(celli, nullptr);

                splitCell* masterPtr = new splitCell(celli, parentPtr);

                splitCell* slavePtr = new splitCell(addedCelli, parentPtr);

                // Store newly created cells on parent together with face
                // that splits them
                parentPtr->master() = masterPtr;
                parentPtr->slave() = slavePtr;

                // Insert master and slave into live splitcell list

                if (liveSplitCells_.found(addedCelli))
                {
                    FatalErrorInFunction
                        << "problem addedCell:" << addedCelli
                        << abort(FatalError);
                }

                liveSplitCells_.insert(celli, masterPtr);
                liveSplitCells_.insert(addedCelli, slavePtr);
            }
            else
            {
                // Cell that was split has been split again.
                splitCell* parentPtr = findCell();

                // It is no longer live
                liveSplitCells_.erase(findCell);

                splitCell* masterPtr = new splitCell(celli, parentPtr);

                splitCell* slavePtr = new splitCell(addedCelli, parentPtr);

                // Store newly created cells on parent together with face
                // that splits them
                parentPtr->master() = masterPtr;
                parentPtr->slave() = slavePtr;

                // Insert master and slave into live splitcell list

                if (liveSplitCells_.found(addedCelli))
                {
                    FatalErrorInFunction
                        << "problem addedCell:" << addedCelli
                        << abort(FatalError);
                }

                liveSplitCells_.insert(celli, masterPtr);
                liveSplitCells_.insert(addedCelli, slavePtr);
            }
        }

        if (debug & 2)
        {
            Pout<< "** After refinement: liveSplitCells_:" << endl;

            printRefTree(Pout);
        }
    }
}


void Foam::undoableMeshCutter::updateMesh(const mapPolyMesh& morphMap)
{
    // Update mesh cutter for new labels.
    meshCutter::updateMesh(morphMap);

    // No need to update cell walker for new labels since does not store any.

    // Update faceRemover for new labels
    faceRemover_.updateMesh(morphMap);

    if (undoable_)
    {
        // Update all live split cells for mesh mapper.
        updateLabels(morphMap.reverseCellMap(), liveSplitCells_);
    }
}


Foam::labelList Foam::undoableMeshCutter::getSplitFaces() const
{
    if (!undoable_)
    {
        FatalErrorInFunction
            << "Only call if constructed with unrefinement capability"
            << abort(FatalError);
    }

    DynamicList<label> liveSplitFaces(liveSplitCells_.size());

    forAllConstIter(Map<splitCell*>, liveSplitCells_, iter)
    {
        const splitCell* splitPtr = iter();

        if (!splitPtr->parent())
        {
            FatalErrorInFunction
                << "Live split cell without parent" << endl
                << "splitCell:" << splitPtr->cellLabel()
                << abort(FatalError);
        }

        // Check if not top of refinement and whether it is the master side
        if (splitPtr->isMaster())
        {
            splitCell* slavePtr = splitPtr->getOther();

            if
            (
                liveSplitCells_.found(slavePtr->cellLabel())
             && splitPtr->isUnrefined()
             && slavePtr->isUnrefined()
            )
            {
                // Both master and slave are live and are not refined.
                // Find common face.

                label celli = splitPtr->cellLabel();

                label slaveCelli = slavePtr->cellLabel();

                label commonFacei =
                    meshTools::getSharedFace
                    (
                        mesh(),
                        celli,
                        slaveCelli
                    );

                liveSplitFaces.append(commonFacei);
            }
        }
    }

    return liveSplitFaces.shrink();
}


Foam::Map<Foam::label> Foam::undoableMeshCutter::getAddedCells() const
{
    // (code copied from getSplitFaces)

    if (!undoable_)
    {
        FatalErrorInFunction
            << "Only call if constructed with unrefinement capability"
            << abort(FatalError);
    }

    Map<label> addedCells(liveSplitCells_.size());

    forAllConstIter(Map<splitCell*>, liveSplitCells_, iter)
    {
        const splitCell* splitPtr = iter();

        if (!splitPtr->parent())
        {
            FatalErrorInFunction
                << "Live split cell without parent" << endl
                << "splitCell:" << splitPtr->cellLabel()
                << abort(FatalError);
        }

        // Check if not top of refinement and whether it is the master side
        if (splitPtr->isMaster())
        {
            splitCell* slavePtr = splitPtr->getOther();

            if
            (
                liveSplitCells_.found(slavePtr->cellLabel())
             && splitPtr->isUnrefined()
             && slavePtr->isUnrefined()
            )
            {
                // Both master and slave are live and are not refined.
                addedCells.insert(splitPtr->cellLabel(), slavePtr->cellLabel());
            }
        }
    }
    return addedCells;
}


Foam::labelList Foam::undoableMeshCutter::removeSplitFaces
(
    const labelList& splitFaces,
    polyTopoChange& meshMod
)
{
    if (!undoable_)
    {
        FatalErrorInFunction
            << "Only call if constructed with unrefinement capability"
            << abort(FatalError);
    }

    // Check with faceRemover what faces will get removed. Note that this can
    // be more (but never less) than splitFaces provided.
    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;

    faceRemover().compatibleRemoves
    (
        splitFaces,         // pierced faces
        cellRegion,         // per cell -1 or region it is merged into
        cellRegionMaster,   // per region the master cell
        facesToRemove       // new faces to be removed.
    );

    if (facesToRemove.size() != splitFaces.size())
    {
        Pout<< "cellRegion:" << cellRegion << endl;
        Pout<< "cellRegionMaster:" << cellRegionMaster << endl;

        FatalErrorInFunction
            << "Faces to remove:" << splitFaces << endl
            << "to be removed:" << facesToRemove
            << abort(FatalError);
    }


    // Every face removed will result in neighbour and owner being merged
    // into owner.
    forAll(facesToRemove, facesToRemoveI)
    {
        label facei = facesToRemove[facesToRemoveI];

        if (!mesh().isInternalFace(facei))
        {
            FatalErrorInFunction
                << "Trying to remove face that is not internal"
                << abort(FatalError);
        }

        label own = mesh().faceOwner()[facei];

        label nbr = mesh().faceNeighbour()[facei];

        Map<splitCell*>::iterator ownFind = liveSplitCells_.find(own);

        Map<splitCell*>::iterator nbrFind = liveSplitCells_.find(nbr);

        if
        (
            (ownFind == liveSplitCells_.end())
         || (nbrFind == liveSplitCells_.end())
        )
        {
            // Can happen because of removeFaces adding extra faces to
            // original splitFaces
        }
        else
        {
            // Face is original splitFace.

            splitCell* ownPtr = ownFind();

            splitCell* nbrPtr = nbrFind();

            splitCell* parentPtr = ownPtr->parent();

            // Various checks on sanity.

            if (debug)
            {
                Pout<< "Updating for removed splitFace " << facei
                    << " own:" << own <<  " nbr:" << nbr
                    << " ownPtr:" << ownPtr->cellLabel()
                    << " nbrPtr:" << nbrPtr->cellLabel()
                    << endl;
            }
            if (!parentPtr)
            {
                FatalErrorInFunction
                    << "No parent for owner " << ownPtr->cellLabel()
                    << abort(FatalError);
            }

            if (!nbrPtr->parent())
            {
                FatalErrorInFunction
                    << "No parent for neighbour " << nbrPtr->cellLabel()
                    << abort(FatalError);
            }

            if (parentPtr != nbrPtr->parent())
            {
                FatalErrorInFunction
                    << "Owner and neighbour liveSplitCell entries do not have"
                    << " same parent. facei:" << facei << "  owner:" << own
                    << "  ownparent:" << parentPtr->cellLabel()
                    << " neighbour:" << nbr
                    << "  nbrparent:" << nbrPtr->parent()->cellLabel()
                    << abort(FatalError);
            }

            if
            (
                !ownPtr->isUnrefined()
             || !nbrPtr->isUnrefined()
             || parentPtr->isUnrefined()
            )
            {
                // Live owner and neighbour are refined themselves.
                FatalErrorInFunction
                    << "Owner and neighbour liveSplitCell entries are"
                    << " refined themselves or the parent is not refined"
                    << endl
                    << "owner unrefined:" << ownPtr->isUnrefined()
                    << "  neighbour unrefined:" << nbrPtr->isUnrefined()
                    << "  master unrefined:" << parentPtr->isUnrefined()
                    << abort(FatalError);
            }

            // Delete from liveSplitCell
            liveSplitCells_.erase(ownFind);

            //!important: Redo search since ownFind entry deleted.
            liveSplitCells_.erase(liveSplitCells_.find(nbr));

            // Delete entries themselves
            delete ownPtr;
            delete nbrPtr;

            //
            // Update parent:
            //   - has parent itself: is part of split cell. Update cellLabel
            //     with merged cell one.
            //   - has no parent: is start of tree. Completely remove.

            if (parentPtr->parent())
            {
                // Update parent cell label to be new merged cell label
                // (will be owner)
                parentPtr->cellLabel() = own;

                // And insert into live cells (is ok since old entry with
                // own as key has been removed above)
                liveSplitCells_.insert(own, parentPtr);
            }
            else
            {
                // No parent so is start of tree. No need to keep splitCell
                // tree.
                delete parentPtr;
            }
        }
    }

    // Insert all commands to combine cells. Never fails so don't have to
    // test for success.
    faceRemover().setRefinement
    (
        facesToRemove,
        cellRegion,
        cellRegionMaster,
        meshMod
    );

    return facesToRemove;
}


// ************************************************************************* //
