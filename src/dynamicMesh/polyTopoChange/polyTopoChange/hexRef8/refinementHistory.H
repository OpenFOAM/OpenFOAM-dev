/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::refinementHistory

Description
    All refinement history. Used in unrefinement.

    - visibleCells: valid for the current mesh and contains per cell -1
      (cell unrefined) or an index into splitCells_.
    - splitCells: for every split contains the parent (also index into
      splitCells) and optionally a subsplit as 8 indices into splitCells.
      Note that the numbers in splitCells are not cell labels, they are purely
      indices into splitCells.

    E.g. 2 cells, cell 1 gets refined so end up with 9 cells:
    \verbatim
        // splitCells
        9
        (
        -1 (1 2 3 4 5 6 7 8)
        0 0()
        0 0()
        0 0()
        0 0()
        0 0()
        0 0()
        0 0()
        0 0()
        )

        // visibleCells
        9(-1 1 2 3 4 5 6 7 8)
    \endverbatim


    So cell0 (visibleCells=-1) is unrefined.
    Cells 1-8 have all valid splitCells entries which are:
      - parent:0
      - subsplits:0()

    The parent 0 refers back to the splitcell entries.


SourceFiles
    refinementHistory.C

\*---------------------------------------------------------------------------*/

#ifndef refinementHistory_H
#define refinementHistory_H

#include "DynamicList.H"
#include "labelList.H"
#include "FixedList.H"
#include "autoPtr.H"
#include "regIOobject.H"
#include "boolList.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mapPolyMesh;
class mapDistributePolyMesh;


// Forward declaration of friend functions and operators

class refinementHistory;

Istream& operator>>(Istream&, refinementHistory&);
Ostream& operator<<(Ostream&, const refinementHistory&);


/*---------------------------------------------------------------------------*\
                      Class refinementHistory Declaration
\*---------------------------------------------------------------------------*/

class refinementHistory
:
    public regIOobject
{
public:

    class splitCell8
    {
    public:

        // Index to original splitCell this cell was refined off from
        // -1: top level cell
        // -2: free splitCell (so should also be in freeSplitCells_)
        label parent_;

        //- Cells this cell was refined into
        autoPtr<FixedList<label, 8>> addedCellsPtr_;

        //- Construct null (parent = -1)
        splitCell8();

        //- Construct from parent
        splitCell8(const label parent);

        //- Construct from Istream
        splitCell8(Istream& is);

        //- Construct as deep copy
        splitCell8(const splitCell8&);

        //- Copy operator since autoPtr otherwise 'steals' storage.
        void operator=(const splitCell8& s);

        bool operator==(const splitCell8& s) const;

        bool operator!=(const splitCell8& s) const;

        friend Istream& operator>>(Istream&, splitCell8&);
        friend Ostream& operator<<(Ostream&, const splitCell8&);
    };


private:

    // Private Data

        //- Is active?
        bool active_;

        //- Storage for splitCells
        DynamicList<splitCell8> splitCells_;

        //- Unused indices in splitCells
        DynamicList<label> freeSplitCells_;

        //- Currently visible cells. Indices into splitCells.
        labelList visibleCells_;


    // Private Member Functions

        //- Debug write
        static void writeEntry
        (
            const List<splitCell8>&,
            const splitCell8&
        );
        //- Debug write
        static void writeDebug
        (
            const labelList&,
            const List<splitCell8>&
        );

        //- Check consistency of structure, i.e. indices into splitCells_.
        void checkIndices() const;

        //- Allocate a splitCell. Return index in splitCells_.
        label allocateSplitCell(const label parent, const label i);

        //- Free a splitCell.
        void freeSplitCell(const label index);

        //- Mark entry in splitCells. Recursively mark its parent and subs.
        void markSplit
        (
            const label,
            labelList& oldToNew,
            DynamicList<splitCell8>&
        ) const;

        void countProc
        (
            const label index,
            const label newProcNo,
            labelList& splitCellProc,
            labelList& splitCellNum
        ) const;

        // For distribution:

            //- Mark index and all its descendants
            void mark(const label, const label, labelList&) const;

            //- Mark cells according to top parent. Return number of clusters
            //  (set of cells originating from same parent)
            label markCommonCells(labelList& cellToCluster) const;

public:

    // Declare name of the class and its debug switch
    TypeName("refinementHistory");


    // Constructors

        //- Construct (read) given an IOobject. If global number of visible
        //  cells > 0 becomes active
        refinementHistory(const IOobject&);

        //- Construct (read) or construct from components
        refinementHistory
        (
            const IOobject&,
            const List<splitCell8>& splitCells,
            const labelList& visibleCells,
            const bool active
        );

        //- Construct (read) or construct from initial number of cells
        //  (all visible). If global number of visible
        //  cells > 0 becomes active
        refinementHistory(const IOobject&, const label nCells);

        //- Construct (read) or construct from initial number of cells
        //  (all visible) and active flag
        refinementHistory
        (
            const IOobject&,
            const label nCells,
            const bool active
        );

        //- Copy constructor
        refinementHistory(const IOobject&, const refinementHistory&);

        //- Construct from multiple refinement histories. If global number of
        //  visible cells > 0 becomes active
        refinementHistory
        (
            const IOobject&,
            const UPtrList<const labelList>& cellMaps,
            const UPtrList<const refinementHistory>&
        );

        //- Construct from Istream. If global number of
        //  visible cells > 0 becomes active
        refinementHistory(const IOobject&, Istream&);


    // Member Functions


        //- Per cell in the current mesh (i.e. visible) either -1 (unrefined)
        //  or an index into splitCells.
        const labelList& visibleCells() const
        {
            return visibleCells_;
        }

        //- Storage for splitCell8s.
        const DynamicList<splitCell8>& splitCells() const
        {
            return splitCells_;
        }

        //- Cache of unused indices in splitCells
        const DynamicList<label>& freeSplitCells() const
        {
            return freeSplitCells_;
        }

        //- Is there unrefinement history?
        bool active() const
        {
            return active_;
        }

        //- Is there unrefinement history?
        bool& active()
        {
            return active_;
        }

        //- Get parent of cell
        label parentIndex(const label celli) const
        {
            label index = visibleCells_[celli];

            if (index < 0)
            {
                FatalErrorInFunction
                    << "Cell " << celli << " is not visible"
                    << abort(FatalError);
            }
            return splitCells_[index].parent_;
        }

        //- Store splitting of cell into 8
        void storeSplit
        (
            const label celli,
            const labelList& addedCells
        );

        //- Store combining 8 cells into master
        void combineCells
        (
            const label masterCelli,
            const labelList& combinedCells
        );

        //- Low level clone
        autoPtr<refinementHistory> clone
        (
            const IOobject& io,
            const labelList& decomposition,
            const labelList& splitCellProc,
            const labelList& splitCellNum,
            const label procI,
            labelList& oldToNewSplit
        ) const;

        //- Create clone from subset
        autoPtr<refinementHistory> clone
        (
            const IOobject& io,
            const labelList& cellMap
        ) const;

        //- Update numbering for mesh changes
        void updateMesh(const mapPolyMesh&);

        //- Update numbering for subsetting
        void subset
        (
            const labelList& pointMap,
            const labelList& faceMap,
            const labelList& cellMap
        );

        //- Update local numbering for mesh redistribution.
        //  Can only distribute clusters sent across in one go; cannot
        //  handle parts recombined in multiple passes.
        void distribute(const mapDistributePolyMesh&);


        //- Compact splitCells_. Removes all freeSplitCells_ elements.
        void compact();

        //- Extend/shrink storage. additional visibleCells_ elements get
        //  set to -1.
        void resize(const label nCells);

        //- Debug write
        void writeDebug() const;


        //- Read object. If global number of visible cells > 0 becomes active
        virtual bool read();

        //- ReadData function required for regIOobject read operation. Note:
        //  does not do a reduction - does not set active_ flag
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;

        // Helpers for decompositionConstraint

            //- Add my decomposition constraints
            void add
            (
                boolList& blockedFace,
                PtrList<labelList>& specifiedProcessorFaces,
                labelList& specifiedProcessor,
                List<labelPair>& explicitConnections
            ) const;

            //- Apply any additional post-decomposition constraints
            void apply
            (
                const boolList& blockedFace,
                const PtrList<labelList>& specifiedProcessorFaces,
                const labelList& specifiedProcessor,
                const List<labelPair>& explicitConnections,
                labelList& decomposition
            ) const;


    // IOstream Operators

        //- Istream operator. Note: does not do a reduction - does not set
        //  active_ flag
        friend Istream& operator>>(Istream&, refinementHistory&);

        friend Ostream& operator<<(Ostream&, const refinementHistory&);
};


Istream& operator>>(Istream&, refinementHistory::splitCell8&);
Ostream& operator<<(Ostream&, const refinementHistory::splitCell8&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
