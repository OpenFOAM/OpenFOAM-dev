/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::multiDirRefinement

Description
    Does multiple pass refinement to refine cells in multiple directions.

    Gets a list of cells to refine and vectorFields for the whole mesh.
    It then tries to refine in one direction after the other the wanted cells.
    After construction the mesh will have been refined in multiple directions.

    Holds the list of cells to refine and the map from original to added for
    every refinement level.

    Gets constructed from a dictionary or from components.
    Uses an undoableMeshCutter which does the actual cutting. Undo facility
    is switched of unless constructed from external one which allows this.

    The cut cells get stored in addedCells which is for every vectorField
    to cut with the map from uncut to added cell (i.e. from master to slave).
    Note: map is only valid for a given direction.

    Parallel: should be ok. Uses 'reduce' whenever it needs to make a
    local decision.

SourceFiles
    multiDirRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef multiDirRefinement_H
#define multiDirRefinement_H

#include "refinementIterator.H"
#include "vectorField.H"
#include "Map.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class undoableMeshCutter;
class cellLooper;
class topoSet;

/*---------------------------------------------------------------------------*\
                           Class multiDirRefinement Declaration
\*---------------------------------------------------------------------------*/

class multiDirRefinement
{
    // Private data

        //- Current set of cells to refine. Extended with added cells.
        labelList cellLabels_;

        //- from original to added cells.
        //  Gives for every cell in the original mesh an empty list or the
        //  list of cells this one has been split into (note: will include
        //  itself so e.g. for hex will be 8 if 2x2x2 refinement)
        labelListList addedCells_;


    // Private Static Functions

        //- Given map from original to added cell set the refineCell for
        //  the added cells to be equal to the one on the original cells.
        static void addCells(const Map<label>&, List<refineCell>&);

        //- Given map from original to added cell set the vectorField for
        //  the added cells to be equal to the one on the original cells.
        static void update(const Map<label>&, vectorField&);

        //- Given map from original to added cell add the added cell to the
        //  list of labels
        static void addCells(const Map<label>&, labelList& labels);


    // Private Member Functions

        //- Add new cells from map to overall list (addedCells_).
        void addCells(const primitiveMesh&, const Map<label>&);

        //- Remove hexes from cellLabels_ and return these in a list.
        labelList splitOffHex(const primitiveMesh& mesh);


        //- Refine cells (hex only) in all 3 directions.
        void refineHex8
        (
            polyMesh& mesh,
            const labelList& hexCells,
            const bool writeMesh
        );

        //- Refine cells in cellLabels_ in directions mentioned.
        void refineAllDirs
        (
            polyMesh& mesh,
            List<vectorField>& cellDirections,
            const cellLooper& cellWalker,
            undoableMeshCutter& cutter,
            const bool writeMesh
        );

        //- Refine based on dictionary. Calls refineAllDirs.
        void refineFromDict
        (
            polyMesh& mesh,
            List<vectorField>& cellDirections,
            const dictionary& dict,
            const bool writeMesh
        );


        //- Disallow default bitwise copy construct
        multiDirRefinement(const multiDirRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const multiDirRefinement&);


public:

    //- Runtime type information
    ClassName("multiDirRefinement");


    // Constructors

        //- Construct from dictionary. After construction all refinement will
        //  have been done (and runTime will have increased a few time steps if
        //  writeMesh = true)
        multiDirRefinement
        (
            polyMesh& mesh,
            const labelList& cellLabels,    // cells to refine
            const dictionary& dict
        );

        //- Explicitly provided directions to split in.
        multiDirRefinement
        (
            polyMesh& mesh,
            const labelList& cellLabels,    // cells to refine
            const List<vectorField>&,       // Explicitly provided directions
            const dictionary& dict
        );

        //- Construct from components. Only this one would allow undo actions.
        multiDirRefinement
        (
            polyMesh& mesh,
            undoableMeshCutter& cutter,     // actual mesh modifier
            const cellLooper& cellCutter,   // how to cut a single cell with
                                            // a plane
            const labelList& cellLabels,    // list of cells to refine
            const List<vectorField>& directions,
            const bool writeMesh = false    // write intermediate meshes
        );


    // Member Functions

        //- Access to addedCells (on the original mesh; see above)
        const labelListList& addedCells() const
        {
            return addedCells_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
