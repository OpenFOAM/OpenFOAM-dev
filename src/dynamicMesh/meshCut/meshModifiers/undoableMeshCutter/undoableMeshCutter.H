/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::undoableMeshCutter

Description
    The main refinement handler. Gets cellCuts which is structure that
    describes which cells are to be cut and in what way. Maintains an undo
    list (if told so during construction). Apart from undo list is just
    wrapper around meshCutter.

    Undo list: contains a refinement tree (of type splitCell; cell labels are
    of no consequence) and a list of visible splitCells, i.e. the top of the
    tree (where the cell labels are valid). Now every cell added gets put on
    the tree and every updateMesh action updates the labels of visible
    splitcells.

    We can now ask this structure for a list of visible split cells or the list
    of faces between these. These can be passed to removeFaces for actual
    deletion and we delete the top splitCell and update the now newly visible
    underlying cells for the new cell number (passed back from removeFaces).

    NOTE: Undoing note properly tested. Expect it to fail if the faces to
    be removed cause other faces to be additionally removed (i.e. removeFaces
    adds additional faces to remove).

    splitCell:
    - original cell number.
    - pointer to parent (null for first level splitCell)
    - two pointers to splitCell children. Both null (unrefined=visible cell) or
      both non-null.

    - live are:
        (-all unrefined cells (original cell without any splitCells))
        -all splitCells with null children

    - liveSplitCells contains pointers to splitCells with null children.



SourceFiles
    undoableMeshCutter.C

\*---------------------------------------------------------------------------*/

#ifndef undoableMeshCutter_H
#define undoableMeshCutter_H

#include "edgeVertex.H"
#include "refineCell.H"
#include "boolList.H"
#include "cellLooper.H"
#include "meshCutter.H"
#include "Map.H"
#include "typeInfo.H"
#include "removeFaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class refineCell;
class splitCell;

/*---------------------------------------------------------------------------*\
                     Class undoableMeshCutter Declaration
\*---------------------------------------------------------------------------*/

class undoableMeshCutter
:
    public meshCutter
{
    // Private data

        //- Whether or not to store actions for unplaying.
        const bool undoable_;

        //- Current split cells which are 'visible'. Only set if undoable.
        Map<splitCell*> liveSplitCells_;

        //- Face remover engine
        removeFaces faceRemover_;


    // Private Member Functions

        //- Debug print
        void printCellRefTree(Ostream& os, const word&, const splitCell*)
         const;

        //- Debug print
        void printRefTree(Ostream& os) const;

        //- Find shared face between two cells
        label sharedFace
        (
            const label cell0I,
            const label cell1I
        ) const;


        //- Update labels on splitCell structure after morphing.
        static void updateLabels(const labelList& map, Map<splitCell*>&);


        //- Disallow default bitwise copy construct
        undoableMeshCutter(const undoableMeshCutter&);

        //- Disallow default bitwise assignment
        void operator=(const undoableMeshCutter&);


public:

    //- Runtime type information
    ClassName("undoableMeshCutter");



    // Constructors

        //- Construct from mesh and flag whether refinement pattern needs
        //  to be stored.
        undoableMeshCutter(const polyMesh& mesh, const bool undoable = true);


    //- Destructor
    ~undoableMeshCutter();


    // Member Functions

        // Access

            //- All current live split cells. Warning: cell labels will change
            //  during morphing. Only this map is guaranteed to hold uptodate
            //  info.
            const Map<splitCell*>& liveSplitCells() const
            {
                return liveSplitCells_;
            }

            const removeFaces& faceRemover() const
            {
                return faceRemover_;
            }


        // Edit

            //- Refine cells acc. to cellCuts. Plays topology changes
            //  into polyTopoChange.
            void setRefinement(const cellCuts& cuts, polyTopoChange&);

            //- Update stored refinement pattern for changes to mesh. Only
            //  call if undoable set.
            void updateMesh(const mapPolyMesh& morphMap);

            //- Calculate split faces from current liveCells. Only
            //  call if undoable set.
            labelList getSplitFaces() const;

            //- Like getSplitFaces but returns map from original to added cell.
            //  Only call if undoable set.
            Map<label> getAddedCells() const;

            //- Remove some refinement. Needs to be supplied subset of
            //  getSplitFaces() output. Returns list of faces removed
            //  (can be more or equal but never less than splitFaces - since
            //  removeFaces might decide to take down unnecessary faces)
            //  Only call if undoable set.
            labelList removeSplitFaces
            (
                const labelList& splitFaces,
                polyTopoChange&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
