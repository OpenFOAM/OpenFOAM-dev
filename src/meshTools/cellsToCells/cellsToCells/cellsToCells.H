/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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
    Foam::cellsToCells

Description
    Class to calculate interpolative addressing and weights between the cells
    of two overlapping meshes

SourceFiles
    cellsToCells.C
    cellsToCellsParallelOps.C
    cellsToCellsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cellsToCells_H
#define cellsToCells_H

#include "remote.H"
#include "distributionMap.H"
#include "polyMesh.H"
#include "runTimeSelectionTables.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cellsToCells Declaration
\*---------------------------------------------------------------------------*/

class cellsToCells
{
protected:

    // Protected Data

        //- Index of processor that holds all of both sides, or -1
        label singleProcess_;

        //- For each source cell, the coupled local target cells
        labelListList srcLocalTgtCells_;

        //- For each target cell, the coupled local source cells
        labelListList tgtLocalSrcCells_;

        //- For each source cell, the coupled target weights
        scalarListList srcWeights_;

        //- For each target cell, the coupled source weights
        scalarListList tgtWeights_;

        //- Map from source patch cells to target-local source patch cells
        autoPtr<distributionMap> srcMapPtr_;

        //- Map from target patch cells to source-local target patch cells
        autoPtr<distributionMap> tgtMapPtr_;

        //- When running in parallel, a map from local source cell index to
        //  source processor and cell index
        autoPtr<List<remote>> localSrcProcCellsPtr_;

        //- When running in parallel, a map from local target cell index to
        //  target processor and cell index
        autoPtr<List<remote>> localTgtProcCellsPtr_;

        //- The target mesh, distributed locally to the source
        autoPtr<polyMesh> localTgtMeshPtr_;


    // Protected Member Functions

        // Intersection

            //- Initialise the addressing and weights
            void initialise
            (
                const polyMesh& srcMesh,
                const polyMesh& tgtMesh
            );

            //- Calculate the addressing and weights
            virtual scalar calculate
            (
                const polyMesh& srcMesh,
                const polyMesh& tgtMesh
            ) = 0;

            //- Normalise the weights for a given mesh
            virtual void normalise
            (
                const polyMesh& mesh,
                labelListList& localOtherCells,
                scalarListList& weights
            ) const = 0;


        // Helpers

            //- Return src cell IDs for the overlap region
            labelList maskCells
            (
                const polyMesh& srcMesh,
                const polyMesh& tgtMesh
            ) const;

            //- Append target cell neighbour cells to cellIDs list
            void appendNbrCells
            (
                const label tgtCelli,
                const polyMesh& mesh,
                const DynamicList<label>& visitedTgtCells,
                DynamicList<label>& nbrTgtCellIDs
            ) const;


        // Parallel operations

            //- Determine which target cells need to be sent to the source.
            //  This is done before intersection. Bound boxes are used to
            //  estimate what cells will intersect.
            labelListList tgtMeshSendCells
            (
                const polyMesh& srcMesh,
                const polyMesh& tgtMesh
            ) const;

            //- Distribute a mesh given its distribution map
            static List<remote> distributeMesh
            (
                const distributionMap& map,
                const polyMesh& mesh,
                autoPtr<polyMesh>& localMeshPtr
            );

            //- Trim the local target addressing and mesh so that communication
            //  from the target to the source is optimised
            void trimLocalTgt();


public:

    //- Run-time type information
    TypeName("cellsToCells");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        cellsToCells,
        word,
        (),
        ()
    );


    // Constructors

        //- Construct null
        cellsToCells();

        //- Disallow default bitwise copy construction
        cellsToCells(const cellsToCells&) = delete;


    //- Destructor
    virtual ~cellsToCells();


    // Selector

        //- Select from name
        static autoPtr<cellsToCells> New(const word& cellsToCellsType);


    // Member Functions

        // Access

            //- Index of the processor holding all cells of the cellsToCells,
            //  or -1 if spread across multiple processors
            inline label singleProcess() const;

            //- Is this intersection on a single process?
            inline bool isSingleProcess() const;

            //- Return a list indicating which source cells are coupled
            PackedBoolList srcCoupled() const;

            //- Return a list indicating which target cells are coupled
            PackedBoolList tgtCoupled() const;


        // Interpolation

            //- Interpolate a source cell field to the target with no left
            //  over values specified. If the interpolation weight sum is less
            //  than one for a face then they will be normalised. If the
            //  interpolation weight sum is zero for a face then that face's
            //  value will be NaN.
            template<class Type>
            tmp<Field<Type>> srcToTgt(const Field<Type>& srcFld) const;

            //- Interpolate a source cell field to the target with left over
            //  values specified. If the interpolation weight sum is less than
            //  one for a face then the average will include the left over
            //  value multiplied by one minus the weight sum.
            template<class Type>
            tmp<Field<Type>> srcToTgt
            (
                const Field<Type>& srcFld,
                const Field<Type>& leftOverTgtFld
            ) const;

            //- Interpolate a target cell field to the source with no left
            //  over values specified. As the corresponding srcToTgt.
            template<class Type>
            tmp<Field<Type>> tgtToSrc(const Field<Type>& tgtFld) const;

            //- Interpolate a target cell field to the source with left
            //  over values specified. As the corresponding srcToTgt.
            template<class Type>
            tmp<Field<Type>> tgtToSrc
            (
                const Field<Type>& tgtFld,
                const Field<Type>& leftOverSrcFld
            ) const;


        // Source-to-target point finding

            //- Find the target processor and cell associated with a point in a
            //  source cell. Note that this will only work with derivations
            //  that fill a cell's stencil with everything that overlaps that
            //  cell. qt present this is just cellsToCellss::intersection, but
            //  we might add a cheaper bound-box based method like
            //  patchToPatches::rays in future.
            remote srcToTgtPoint
            (
                const polyMesh& tgtMesh,
                const label srcCelli,
                const point& p
            ) const;


        // Manipulation

            //- Update addressing and weights for the given meshes. Returns the
            //  overlapping volume (if that is relevant to the method).
            scalar update
            (
                const polyMesh& srcMesh,
                const polyMesh& tgtMesh
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellsToCells&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellsToCellsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cellsToCellsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
