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

\*---------------------------------------------------------------------------*/

#include "DynamicList.H"
#include "refinementHistory.H"
#include "ListOps.H"
#include "mapPolyMesh.H"
#include "mapDistributePolyMesh.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(refinementHistory, 0);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::refinementHistory::writeEntry
(
    const List<splitCell8>& splitCells,
    const splitCell8& split
)
{
    // Write me:
    if (split.addedCellsPtr_.valid())
    {
        Pout<< "parent:" << split.parent_
            << " subCells:" << split.addedCellsPtr_()
            << endl;
    }
    else
    {
        Pout<< "parent:" << split.parent_
            << " no subcells"
            << endl;
    }

    if (split.parent_ >= 0)
    {
        Pout<< "parent data:" << endl;
        // Write my parent
        string oldPrefix = Pout.prefix();
        Pout.prefix() = "  " + oldPrefix;
        writeEntry(splitCells, splitCells[split.parent_]);
        Pout.prefix() = oldPrefix;
    }
}


void Foam::refinementHistory::writeDebug
(
    const labelList& visibleCells,
    const List<splitCell8>& splitCells
)
{
    string oldPrefix = Pout.prefix();
    Pout.prefix() = "";

    forAll(visibleCells, cellI)
    {
        label index = visibleCells[cellI];

        if (index >= 0)
        {
            Pout<< "Cell from refinement:" << cellI << " index:" << index
                << endl;

            string oldPrefix = Pout.prefix();
            Pout.prefix() = "  " + oldPrefix;
            writeEntry(splitCells, splitCells[index]);
            Pout.prefix() = oldPrefix;
        }
        else
        {
            Pout<< "Unrefined cell:" << cellI << " index:" << index << endl;
        }
    }
    Pout.prefix() = oldPrefix;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
Foam::refinementHistory::splitCell8::splitCell8()
:
    parent_(-1),
    addedCellsPtr_(NULL)
{}


//- Construct as child element of parent
Foam::refinementHistory::splitCell8::splitCell8(const label parent)
:
    parent_(parent),
    addedCellsPtr_(NULL)
{}


//- Construct from Istream
Foam::refinementHistory::splitCell8::splitCell8(Istream& is)
{
    is >> *this;
}


//- Construct as (deep) copy.
Foam::refinementHistory::splitCell8::splitCell8(const splitCell8& sc)
:
    parent_(sc.parent_),
    addedCellsPtr_
    (
        sc.addedCellsPtr_.valid()
      ? new FixedList<label, 8>(sc.addedCellsPtr_())
      : NULL
    )
{}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, refinementHistory::splitCell8& sc)
{
    labelList addedCells;

    is >> sc.parent_ >> addedCells;

    if (addedCells.size())
    {
        sc.addedCellsPtr_.reset(new FixedList<label, 8>(addedCells));
    }
    else
    {
        sc.addedCellsPtr_.reset(NULL);
    }

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const refinementHistory::splitCell8& sc
)
{
    // Output as labelList so we can have 0 sized lists. Alternative is to
    // output as fixedlist with e.g. -1 elements and check for this upon
    // reading. However would cause much more data to be transferred.

    if (sc.addedCellsPtr_.valid())
    {
        return os
            << sc.parent_
            << token::SPACE
            << labelList(sc.addedCellsPtr_());
    }
    else
    {
        return os << sc.parent_ << token::SPACE << labelList(0);
    }
}


void Foam::refinementHistory::checkIndices() const
{
    // Check indices.
    forAll(visibleCells_, i)
    {
        if (visibleCells_[i] < 0 && visibleCells_[i] >= splitCells_.size())
        {
            FatalErrorIn("refinementHistory::checkIndices() const")
                << "Illegal entry " << visibleCells_[i]
                << " in visibleCells at location" << i << nl
                << "It points outside the range of splitCells : 0.."
                << splitCells_.size()-1
                << abort(FatalError);
        }
    }
}


Foam::label Foam::refinementHistory::allocateSplitCell
(
    const label parent,
    const label i
)
{
    label index = -1;

    if (freeSplitCells_.size())
    {
        index = freeSplitCells_.remove();

        splitCells_[index] = splitCell8(parent);
    }
    else
    {
        index = splitCells_.size();

        splitCells_.append(splitCell8(parent));
    }


    // Update the parent field
    if (parent >= 0)
    {
        splitCell8& parentSplit = splitCells_[parent];

        if (parentSplit.addedCellsPtr_.empty())
        {
            // Allocate storage on parent for the 8 subcells.
            parentSplit.addedCellsPtr_.reset(new FixedList<label, 8>(-1));
        }


        // Store me on my parent
        FixedList<label, 8>& parentSplits = parentSplit.addedCellsPtr_();

        parentSplits[i] = index;
    }

    return index;
}


void Foam::refinementHistory::freeSplitCell(const label index)
{
    splitCell8& split = splitCells_[index];

    // Make sure parent does not point to me anymore.
    if (split.parent_ >= 0)
    {
        autoPtr<FixedList<label, 8> >& subCellsPtr =
            splitCells_[split.parent_].addedCellsPtr_;

        if (subCellsPtr.valid())
        {
            FixedList<label, 8>& subCells = subCellsPtr();

            label myPos = findIndex(subCells, index);

            if (myPos == -1)
            {
                FatalErrorIn("refinementHistory::freeSplitCell")
                    << "Problem: cannot find myself in"
                    << " parents' children" << abort(FatalError);
            }
            else
            {
                subCells[myPos] = -1;
            }
        }
    }

    // Mark splitCell as free
    split.parent_ = -2;

    // Add to cache of free splitCells
    freeSplitCells_.append(index);
}


// Mark entry in splitCells. Recursively mark its parent and subs.
void Foam::refinementHistory::markSplit
(
    const label index,
    labelList& oldToNew,
    DynamicList<splitCell8>& newSplitCells
) const
{
    if (oldToNew[index] == -1)
    {
        // Not yet compacted.

        const splitCell8& split = splitCells_[index];

        oldToNew[index] = newSplitCells.size();
        newSplitCells.append(split);

        if (split.parent_ >= 0)
        {
            markSplit(split.parent_, oldToNew, newSplitCells);
        }
        if (split.addedCellsPtr_.valid())
        {
            const FixedList<label, 8>& splits = split.addedCellsPtr_();

            forAll(splits, i)
            {
                if (splits[i] >= 0)
                {
                    markSplit(splits[i], oldToNew, newSplitCells);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementHistory::refinementHistory(const IOobject& io)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
            "refinementHistory::refinementHistory(const IOobject&)"
        )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " does not support automatic rereading."
            << endl;
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }

    if (debug)
    {
        Pout<< "refinementHistory::refinementHistory :"
            << " constructed history from IOobject :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


//- Read or construct
Foam::refinementHistory::refinementHistory
(
    const IOobject& io,
    const List<splitCell8>& splitCells,
    const labelList& visibleCells
)
:
    regIOobject(io),
    splitCells_(splitCells),
    freeSplitCells_(0),
    visibleCells_(visibleCells)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
            "refinementHistory::refinementHistory"
            "(const IOobject&, const List<splitCell8>&, const labelList&)"
        )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " does not support automatic rereading."
            << endl;
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementHistory::refinementHistory :"
            << " constructed history from IOobject or components :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// Construct from initial number of cells (all visible)
Foam::refinementHistory::refinementHistory
(
    const IOobject& io,
    const label nCells
)
:
    regIOobject(io),
    freeSplitCells_(0)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
            "refinementHistory::refinementHistory"
            "(const IOobject&, const label)"
        )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " does not support automatic rereading."
            << endl;
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
    else
    {
        visibleCells_.setSize(nCells);
        splitCells_.setCapacity(nCells);

        for (label cellI = 0; cellI < nCells; cellI++)
        {
            visibleCells_[cellI] = cellI;
            splitCells_.append(splitCell8());
        }
    }

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementHistory::refinementHistory :"
            << " constructed history from IOobject or initial size :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// Construct as copy
Foam::refinementHistory::refinementHistory
(
    const IOobject& io,
    const refinementHistory& rh
)
:
    regIOobject(io),
    splitCells_(rh.splitCells()),
    freeSplitCells_(rh.freeSplitCells()),
    visibleCells_(rh.visibleCells())
{
    if (debug)
    {
        Pout<< "refinementHistory::refinementHistory : constructed initial"
            << " history." << endl;
    }
}


// Construct from Istream
Foam::refinementHistory::refinementHistory(const IOobject& io, Istream& is)
:
    regIOobject(io),
    splitCells_(is),
    freeSplitCells_(0),
    visibleCells_(is)
{
    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementHistory::refinementHistory :"
            << " constructed history from Istream"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::refinementHistory::resize(const label size)
{
    label oldSize = visibleCells_.size();

    if (debug)
    {
        Pout<< "refinementHistory::resize from " << oldSize << " to " << size
            << " cells" << endl;
    }

    visibleCells_.setSize(size);

    // Set additional elements to -1.
    for (label i = oldSize; i < visibleCells_.size(); i++)
    {
        visibleCells_[i] = -1;
    }
}


void Foam::refinementHistory::updateMesh(const mapPolyMesh& map)
{
    if (active())
    {
        const labelList& reverseCellMap = map.reverseCellMap();

        // Note that only the live cells need to be renumbered.

        labelList newVisibleCells(map.cellMap().size(), -1);

        forAll(visibleCells_, cellI)
        {
            if (visibleCells_[cellI] != -1)
            {
                label index = visibleCells_[cellI];

                // Check not already set
                if (splitCells_[index].addedCellsPtr_.valid())
                {
                    FatalErrorIn
                    (
                        "refinementHistory::updateMesh(const mapPolyMesh&)"
                    )   << "Problem" << abort(FatalError);
                }

                label newCellI = reverseCellMap[cellI];

                if (newCellI >= 0)
                {
                    newVisibleCells[newCellI] = index;
                }
            }
        }

        if (debug)
        {
            Pout<< "refinementHistory::updateMesh : from "
                << visibleCells_.size()
                << " to " << newVisibleCells.size()
                << " cells" << endl;
        }

        visibleCells_.transfer(newVisibleCells);
    }
}


// Update numbering for subsetting
void Foam::refinementHistory::subset
(
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap
)
{
    if (active())
    {
        labelList newVisibleCells(cellMap.size(), -1);

        forAll(newVisibleCells, cellI)
        {
            label oldCellI = cellMap[cellI];

            label index = visibleCells_[oldCellI];

            // Check that cell is live (so its parent has no refinement)
            if (index >= 0 && splitCells_[index].addedCellsPtr_.valid())
            {
                FatalErrorIn
                (
                    "refinementHistory::subset"
                    "(const labelList&, const labelList&, const labelList&)"
                )   << "Problem" << abort(FatalError);
            }

            newVisibleCells[cellI] = index;
        }

        if (debug)
        {
            Pout<< "refinementHistory::updateMesh : from "
                << visibleCells_.size()
                << " to " << newVisibleCells.size()
                << " cells" << endl;
        }

        visibleCells_.transfer(newVisibleCells);
    }
}


void Foam::refinementHistory::countProc
(
    const label index,
    const label newProcNo,
    labelList& splitCellProc,
    labelList& splitCellNum
) const
{
    if (splitCellProc[index] != newProcNo)
    {
        // Different destination processor from other cells using this
        // parent. Reset count.
        splitCellProc[index] = newProcNo;
        splitCellNum[index] = 1;
    }
    else
    {
        splitCellNum[index]++;

        // Increment parent if whole splitCell moves to same processor
        if (splitCellNum[index] == 8)
        {
            if (debug)
            {
                Pout<< "Moving " << splitCellNum[index]
                    << " cells originating from cell " << index
                    << " from processor " << Pstream::myProcNo()
                    << " to processor " << splitCellProc[index]
                    << endl;
            }

            label parent = splitCells_[index].parent_;

            if (parent >= 0)
            {
                countProc(parent, newProcNo, splitCellProc, splitCellNum);
            }
        }
    }
}


void Foam::refinementHistory::distribute(const mapDistributePolyMesh& map)
{
    if (!active())
    {
        FatalErrorIn
        (
            "refinementHistory::distribute(const mapDistributePolyMesh&)"
        )   << "Calling distribute on inactive history" << abort(FatalError);
    }


    if (!Pstream::parRun())
    {
        return;
    }

    // Remove unreferenced history.
    compact();

    //Pout<< nl << "--BEFORE:" << endl;
    //writeDebug();
    //Pout<< "---------" << nl << endl;


    // Distribution is only partially functional.
    // If all 8 cells resulting from a single parent are sent across in one
    // go it will also send across that part of the refinement history.
    // If however e.g. first 1 and then the other 7 are sent across the
    // history will not be reconstructed.

    // Determine clusters. This is per every entry in splitCells_ (that is
    // a parent of some refinement) a label giving the processor it goes to
    // if all its children are going to the same processor.

    // Per visible cell the processor it goes to.
    labelList destination(visibleCells_.size());

    const labelListList& subCellMap = map.cellMap().subMap();

    forAll(subCellMap, procI)
    {
        const labelList& newToOld = subCellMap[procI];

        forAll(newToOld, i)
        {
            label oldCellI = newToOld[i];

            destination[oldCellI] = procI;
        }
    }

//Pout<< "refinementHistory::distribute :"
//    << " destination:" << destination << endl;

    // Per splitCell entry the processor it moves to
    labelList splitCellProc(splitCells_.size(), -1);
    // Per splitCell entry the number of live cells that move to that processor
    labelList splitCellNum(splitCells_.size(), 0);

    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            countProc
            (
                splitCells_[index].parent_,
                destination[cellI],
                splitCellProc,
                splitCellNum
            );
        }
    }

    //Pout<< "refinementHistory::distribute :"
    //    << " splitCellProc:" << splitCellProc << endl;
    //
    //Pout<< "refinementHistory::distribute :"
    //    << " splitCellNum:" << splitCellNum << endl;


    // Create subsetted refinement tree consisting of all parents that
    // move in their whole to other processor.
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        //Pout<< "-- Subetting for processor " << procI << endl;

        // From uncompacted to compacted splitCells.
        labelList oldToNew(splitCells_.size(), -1);

        // Compacted splitCells. Similar to subset routine below.
        DynamicList<splitCell8> newSplitCells(splitCells_.size());

        // Loop over all entries. Note: could recurse like countProc so only
        // visit used entries but is probably not worth it.

        forAll(splitCells_, index)
        {
//            Pout<< "oldCell:" << index
//                << " proc:" << splitCellProc[index]
//                << " nCells:" << splitCellNum[index]
//                << endl;

            if (splitCellProc[index] == procI && splitCellNum[index] == 8)
            {
                // Entry moves in its whole to procI
                oldToNew[index] = newSplitCells.size();
                newSplitCells.append(splitCells_[index]);

                //Pout<< "Added oldCell " << index
                //    << " info " << newSplitCells.last()
                //    << " at position " << newSplitCells.size()-1
                //    << endl;
            }
        }

        // Add live cells that are subsetted.
        forAll(visibleCells_, cellI)
        {
            label index = visibleCells_[cellI];

            if (index >= 0 && destination[cellI] == procI)
            {
                label parent = splitCells_[index].parent_;

                //Pout<< "Adding refined cell " << cellI
                //    << " since moves to "
                //    << procI << " old parent:" << parent << endl;

                // Create new splitCell with parent
                oldToNew[index] = newSplitCells.size();
                newSplitCells.append(splitCell8(parent));
            }
        }

        //forAll(oldToNew, index)
        //{
        //    Pout<< "old:" << index << " new:" << oldToNew[index]
        //        << endl;
        //}

        newSplitCells.shrink();

        // Renumber contents of newSplitCells
        forAll(newSplitCells, index)
        {
            splitCell8& split = newSplitCells[index];

            if (split.parent_ >= 0)
            {
                split.parent_ = oldToNew[split.parent_];
            }
            if (split.addedCellsPtr_.valid())
            {
                FixedList<label, 8>& splits = split.addedCellsPtr_();

                forAll(splits, i)
                {
                    if (splits[i] >= 0)
                    {
                        splits[i] = oldToNew[splits[i]];
                    }
                }
            }
        }


        const labelList& subMap = subCellMap[procI];

        // New visible cells.
        labelList newVisibleCells(subMap.size(), -1);

        forAll(subMap, newCellI)
        {
            label oldCellI = subMap[newCellI];

            label oldIndex = visibleCells_[oldCellI];

            if (oldIndex >= 0)
            {
                newVisibleCells[newCellI] = oldToNew[oldIndex];
            }
        }

        //Pout<< nl << "--Subset for domain:" << procI << endl;
        //writeDebug(newVisibleCells, newSplitCells);
        //Pout<< "---------" << nl << endl;


        // Send to neighbours
        OPstream toNbr(Pstream::blocking, procI);
        toNbr << newSplitCells << newVisibleCells;
    }


    // Receive from neighbours and merge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Remove all entries. Leave storage intact.
    splitCells_.clear();

    visibleCells_.setSize(map.mesh().nCells());
    visibleCells_ = -1;

    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        IPstream fromNbr(Pstream::blocking, procI);
        List<splitCell8> newSplitCells(fromNbr);
        labelList newVisibleCells(fromNbr);

        //Pout<< nl << "--Received from domain:" << procI << endl;
        //writeDebug(newVisibleCells, newSplitCells);
        //Pout<< "---------" << nl << endl;


        // newSplitCells contain indices only into newSplitCells so
        // renumbering can be done here.
        label offset = splitCells_.size();

        //Pout<< "**Renumbering data from proc " << procI << " with offset "
        //    << offset << endl;

        forAll(newSplitCells, index)
        {
            splitCell8& split = newSplitCells[index];

            if (split.parent_ >= 0)
            {
                split.parent_ += offset;
            }
            if (split.addedCellsPtr_.valid())
            {
                FixedList<label, 8>& splits = split.addedCellsPtr_();

                forAll(splits, i)
                {
                    if (splits[i] >= 0)
                    {
                        splits[i] += offset;
                    }
                }
            }

            splitCells_.append(split);
        }


        // Combine visibleCell.
        const labelList& constructMap = map.cellMap().constructMap()[procI];

        forAll(newVisibleCells, i)
        {
            if (newVisibleCells[i] >= 0)
            {
                visibleCells_[constructMap[i]] = newVisibleCells[i] + offset;
            }
        }
    }
    splitCells_.shrink();

    //Pout<< nl << "--AFTER:" << endl;
    //writeDebug();
    //Pout<< "---------" << nl << endl;
}


void Foam::refinementHistory::compact()
{
    if (debug)
    {
        Pout<< "refinementHistory::compact() Entering with:"
            << " freeSplitCells_:" << freeSplitCells_.size()
            << " splitCells_:" << splitCells_.size()
            << " visibleCells_:" << visibleCells_.size()
            << endl;

        // Check all free splitCells are marked as such
        forAll(freeSplitCells_, i)
        {
            label index = freeSplitCells_[i];

            if (splitCells_[index].parent_ != -2)
            {
                FatalErrorIn("refinementHistory::compact()")
                    << "Problem index:" << index
                    << abort(FatalError);
            }
        }

        // Check none of the visible cells are marked as free
        forAll(visibleCells_, cellI)
        {
            if
            (
                visibleCells_[cellI] >= 0
             && splitCells_[visibleCells_[cellI]].parent_ == -2
            )
            {
                FatalErrorIn("refinementHistory::compact()")
                    << "Problem : visible cell:" << cellI
                    << " is marked as being free." << abort(FatalError);
            }
        }
    }

    DynamicList<splitCell8> newSplitCells(splitCells_.size());

    // From uncompacted to compacted splitCells.
    labelList oldToNew(splitCells_.size(), -1);

    // Mark all used splitCell entries. These are either indexed by visibleCells
    // or indexed from other splitCell entries.

    // Mark from visibleCells
    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            // Make sure we only mark visible indices if they either have a
            // parent or subsplits.
            if
            (
                splitCells_[index].parent_ != -1
             || splitCells_[index].addedCellsPtr_.valid()
            )
            {
                markSplit(index, oldToNew, newSplitCells);
            }
        }
    }

    // Mark from splitCells
    forAll(splitCells_, index)
    {
        if (splitCells_[index].parent_ == -2)
        {
            // freed cell.
        }
        else if
        (
            splitCells_[index].parent_ == -1
         && splitCells_[index].addedCellsPtr_.empty()
        )
        {
            // recombined cell. No need to keep since no parent and no subsplits
            // Note that gets marked if reachable from other index!
        }
        else
        {
            // Is used element.
            markSplit(index, oldToNew, newSplitCells);
        }
    }


    // Now oldToNew is fully complete and compacted elements are in
    // newSplitCells.
    // Renumber contents of newSplitCells and visibleCells.
    forAll(newSplitCells, index)
    {
        splitCell8& split = newSplitCells[index];

        if (split.parent_ >= 0)
        {
            split.parent_ = oldToNew[split.parent_];
        }
        if (split.addedCellsPtr_.valid())
        {
            FixedList<label, 8>& splits = split.addedCellsPtr_();

            forAll(splits, i)
            {
                if (splits[i] >= 0)
                {
                    splits[i] = oldToNew[splits[i]];
                }
            }
        }
    }


    if (debug)
    {
        Pout<< "refinementHistory::compact : compacted splitCells from "
            << splitCells_.size() << " to " << newSplitCells.size() << endl;
    }

    splitCells_.transfer(newSplitCells);
    freeSplitCells_.clearStorage();


    if (debug)
    {
        Pout<< "refinementHistory::compact() NOW:"
            << " freeSplitCells_:" << freeSplitCells_.size()
            << " splitCells_:" << splitCells_.size()
            << " newSplitCells:" << newSplitCells.size()
            << " visibleCells_:" << visibleCells_.size()
            << endl;
    }


    // Adapt indices in visibleCells_
    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            // Note that oldToNew can be -1 so it resets newVisibleCells.
            visibleCells_[cellI] = oldToNew[index];
        }
        else
        {
            // Keep -1 value.
        }
    }
}


void Foam::refinementHistory::writeDebug() const
{
    writeDebug(visibleCells_, splitCells_);
}


void Foam::refinementHistory::storeSplit
(
    const label cellI,
    const labelList& addedCells
)
{
    label parentIndex = -1;

    if (visibleCells_[cellI] != -1)
    {
        // Was already live. The current live cell becomes the
        // parent of the cells split off from it.

        parentIndex = visibleCells_[cellI];

        // It is no longer live (note that actually cellI gets alive
        // again below since is addedCells[0])
        visibleCells_[cellI] = -1;
    }
    else
    {
        // Create 0th level. -1 parent to denote this.
        parentIndex = allocateSplitCell(-1, -1);
    }

    // Create live entries for added cells that point to the
    // cell they were created from (parentIndex)
    forAll(addedCells, i)
    {
        label addedCellI = addedCells[i];

        // Create entries for the split off cells. All of them
        // are visible.
        visibleCells_[addedCellI] = allocateSplitCell(parentIndex, i);
    }
}


void Foam::refinementHistory::combineCells
(
    const label masterCellI,
    const labelList& combinedCells
)
{
    // Save the parent structure
    label parentIndex = splitCells_[visibleCells_[masterCellI]].parent_;

    // Remove the information for the combined cells
    forAll(combinedCells, i)
    {
        label cellI = combinedCells[i];

        freeSplitCell(visibleCells_[cellI]);
        visibleCells_[cellI] = -1;
    }

    splitCell8& parentSplit = splitCells_[parentIndex];
    parentSplit.addedCellsPtr_.reset(NULL);
    visibleCells_[masterCellI] = parentIndex;
}


bool Foam::refinementHistory::readData(Istream& is)
{
    is >> *this;
    return !is.bad();
}


bool Foam::refinementHistory::writeData(Ostream& os) const
{
    os << *this;

    return os.good();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, refinementHistory& rh)
{
    rh.freeSplitCells_.clearStorage();

    is >> rh.splitCells_ >> rh.visibleCells_;

    // Check indices.
    rh.checkIndices();

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const refinementHistory& rh)
{
    const_cast<refinementHistory&>(rh).compact();

    return os   << "// splitCells" << nl
                << rh.splitCells_ << nl
                << "// visibleCells" << nl
                << rh.visibleCells_;
}


// ************************************************************************* //
