/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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
    Foam::meshToMesh

Description
    Class to calculate the cell-addressing between two overlapping meshes

    Mapping is performed using a run-time selectable interpolation mothod

See also
    meshToMeshMethod

SourceFiles
    meshToMesh.C
    meshToMeshParallelOps.C
    meshToMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef meshToMesh_H
#define meshToMesh_H

#include "polyMesh.H"
#include "boundBox.H"
#include "distributionMap.H"
#include "NamedEnum.H"
#include "patchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshToMesh Declaration
\*---------------------------------------------------------------------------*/

class meshToMesh
{
    // Private Data

        //- Reference to the source mesh
        const polyMesh& srcMesh_;

        //- Reference to the target mesh
        const polyMesh& tgtMesh_;

        //- List of corresponding source and target patches that are to be
        //  mapped to each other
        List<labelPair> srcToTgtPatchIDs_;

        //- List of patchToPatch interpolation engines between source and
        //  target patches
        PtrList<patchToPatch> srcToTgtPatchToPatches_;

        //- Cutting patches whose values are set using a zero-gradient condition
        labelList tgtCuttingPatchIDs_;

        //- For each source cell, the coupled local target cells
        labelListList srcLocalTgtCells_;

        //- For each target cell, the coupled local source cells
        labelListList tgtLocalSrcCells_;

        //- For each source cell, the coupled target weights
        scalarListList srcWeights_;

        //- For each target cell, the coupled source weights
        scalarListList tgtWeights_;

        //- Index of processor that holds all of both sides, or -1
        label singleProcess_;

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


    // Private Member Functions

        //- Helper function to add a constant offset to a list
        template<class Type>
        static void add(UList<Type>& fld, const label offset);

        //- Calculate the cell intersections. Returns the overlap volume.
        scalar calculateCellToCells(const word& methodName);

        //- Calculate patch intersections
        void calculatePatchToPatches(const word& patchToPatchType);

        //- Constructor helper
        void constructNoCuttingPatches
        (
            const word& methodName,
            const bool interpAllPatches
        );

        //- Constructor helper
        void constructFromCuttingPatches
        (
            const word& methodName,
            const HashTable<word>& patchMap,
            const wordList& tgtCuttingPatches
        );


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


protected:

    // Protected Member Functions

        // Access

            //- Return the source patch indices
            inline const List<labelPair>& srcToTgtPatchIDs() const;

            //- Return the list of patchToPatch interpolation engines between
            //  source and target patches
            inline const PtrList<patchToPatch>& srcToTgtPatchToPatches() const;

            //- Return the cutting patch indices
            inline const labelList& tgtCuttingPatchIDs() const;


public:

    //- Run-time type information
    TypeName("meshToMesh");


    // Constructors

        //- Construct from source and target meshes, generic mapping methods
        meshToMesh
        (
            const polyMesh& src,
            const polyMesh& tgt,
            const word& methodName,
            const bool interpAllPatches = true
        );

        //- Construct from source and target meshes, generic mapping methods
        meshToMesh
        (
            const polyMesh& src,
            const polyMesh& tgt,
            const word& methodName,
            const HashTable<word>& patchMap,
            const wordList& cuttingPatches
        );

        //- Disallow default bitwise copy construction
        meshToMesh(const meshToMesh&) = delete;


    //- Destructor
    virtual ~meshToMesh();


    // Member Functions

        // Access

            //- Return const access to the source mesh
            inline const polyMesh& srcMesh() const;

            //- Return const access to the target mesh
            inline const polyMesh& tgtMesh() const;

            //- Index of the processor holding all cells of the meshToMesh,
            //  or -1 if spread across multiple processors
            inline label singleProcess() const;

            //- Is this intersection on a single process?
            inline bool isSingleProcess() const;


        // Evaluation

            // Source-to-target field mapping

                //- Map field from src to tgt mesh. Values passed in via
                //  'result' are used to initialise the return value.
                template<class Type>
                void mapSrcToTgt
                (
                    const UList<Type>& srcFld,
                    List<Type>& result
                ) const;

                //- Return the src field mapped to the tgt mesh. Initial values
                //  of the result are set to zero.
                template<class Type>
                tmp<Field<Type>> mapSrcToTgt(const Field<Type>& srcFld) const;


            // Target-to-source field mapping

                //- Map field from tgt to src mesh. Values passed in via
                //  'result' are used to initialise the return value.
                template<class Type>
                void mapTgtToSrc
                (
                    const UList<Type>& tgtFld,
                    List<Type>& result
                ) const;

                //- Return the tgt field mapped to the src mesh. Initial values
                //  of the result are set to zero.
                template<class Type>
                tmp<Field<Type>> mapTgtToSrc(const Field<Type>& tgtFld) const;


            // Source-to-target point finding

                //- Find the target processor and cell associated with a point
                //  in a source cell. Note that this will only work with
                //  methods that fill a cell's stencil with everything that
                //  overlaps that cell. At present this is just
                //  cellVolumeWeightMethod, but we might add a cheaper
                //  bound-box based method like patchToPatches::rays in future.
                remote srcToTgtPoint
                (
                    const label srcCelli,
                    const point& p
                ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshToMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshToMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshToMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
