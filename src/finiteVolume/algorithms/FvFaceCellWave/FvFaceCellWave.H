/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::FvFaceCellWave

Description
    Wave propagation of information through grid. Every iteration
    information goes through one layer of cells. Templated on information
    that is transferred.

    Handles parallel and cyclics and non-parallel cyclics.

    Note: whether to propagate depends on the return value of Type::update
    which returns true (i.e. propagate) if the value changes by more than a
    certain tolerance.
    This tolerance can be very strict for normal face-cell and parallel
    cyclics (we use a value of 0.01 just to limit propagation of small changes)
    but for non-parallel cyclics this tolerance can be critical and if chosen
    too small can lead to non-convergence.

SourceFiles
    FvFaceCellWave.C

\*---------------------------------------------------------------------------*/

#ifndef FvFaceCellWave_H
#define FvFaceCellWave_H

#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class fvPatch;
class transformer;

/*---------------------------------------------------------------------------*\
                      Class FvFaceCellWaveName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(FvFaceCellWave);


/*---------------------------------------------------------------------------*\
                        Class FvFaceCellWave Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class TrackingData = int>
class FvFaceCellWave
:
    public FvFaceCellWaveName
{
protected:

    // Protected static data

        //- ...
        static const scalar geomTol_;

        //- ...
        static scalar propagationTol_;


    // Protected data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Information for internal faces
        List<Type>& internalFaceInfo_;

        //- Information for patch faces
        List<List<Type>>& patchFaceInfo_;

        //- Information for all cells
        List<Type>& cellInfo_;

        //- Additional data to be passed into container
        TrackingData& td_;

        //- Has internal face changed?
        PackedBoolList internalFaceChanged_;

        //- Has patch face changed?
        List<PackedBoolList> patchFaceChanged_;

        //- Has cell changed?
        PackedBoolList cellChanged_;

        //- List of changed patch and faces
        DynamicList<labelPair> changedPatchAndFaces_;

        //- List of changed cells
        DynamicList<label> changedCells_;

        //- Contains cyclics
        const bool hasCyclicPatches_;


    // Protected member functions

        //- ...
        const Type& faceInfo(const labelPair& patchAndFacei) const;

        //- ...
        Type& faceInfo(const labelPair& patchAndFacei);

        //- ...
        bool faceChanged(const labelPair& patchAndFacei) const;

        //- ...
        PackedBoolList::iteratorBase faceChanged(const labelPair&);

        //- Updates cellInfo with information from neighbour. Updates all
        //  statistics.
        bool updateCell
        (
            const label celli,
            const labelPair& neighbourPatchAndFacei,
            const Type& neighbourInfo,
            const scalar tol,
            Type& cellInfo
        );

        //- Updates faceInfo with information from neighbour. Updates all
        //  statistics.
        bool updateFace
        (
            const labelPair& patchAndFacei,
            const label neighbourCelli,
            const Type& neighbourInfo,
            const scalar tol,
            Type& faceInfo
        );

        //- Updates faceInfo with information from same face. Updates all
        //  statistics.
        bool updateFace
        (
            const labelPair& patchAndFacei,
            const Type& neighbourInfo,
            const scalar tol,
            Type& faceInfo
        );


        // Parallel, cyclic

            //- Debugging: check info on both sides of cyclic
            void checkCyclic(const fvPatch& patch) const;

            //- Has cyclic patch?
            template<class PatchType>
            bool hasPatch() const;

            //- Merge received patch data into global data
            void mergeFaceInfo
            (
                const fvPatch& patch,
                const label nFaces,
                const labelList&,
                const List<Type>&
            );

            //- Extract info for single patch only
            label getChangedPatchFaces
            (
                const fvPatch& patch,
                labelList& changedPatchFaces,
                List<Type>& changedPatchFacesInfo
            ) const;

            //- Transform across an interface
            void transform
            (
                const fvPatch& patch,
                const label nFaces,
                const labelList& patchFaces,
                const transformer& transform,
                List<Type>& faceInfo
            );

            //- Merge data from across processor boundaries
            void handleProcPatches();

            //- Merge data from across cyclics
            void handleCyclicPatches();


public:

    // Static Data

        //- Default tracking data to go with default template argument
        static int defaultTrackingData_;


    // Static Functions

        //- Access to tolerance
        static scalar propagationTol()
        {
            return propagationTol_;
        }

        //- Change tolerance
        static void setPropagationTol(const scalar tol)
        {
            propagationTol_ = tol;
        }

        //- ...
        template<class ListList>
        static labelList listListSizes(const ListList& ll);

        //- ...
        template<class ListList, class Value>
        static ListList sizesListList(const labelList& s, const Value& value);


    // Constructors

        // Construct from mesh. Use setFaceInfo and iterate() to do actual
        // calculation.
        FvFaceCellWave
        (
            const fvMesh& mesh,
            List<Type>& internalFaceInfo,
            List<List<Type>>& patchFaceInfo,
            List<Type>& cellInfo,
            TrackingData& td = defaultTrackingData_
        );

        //- Construct from mesh and list of changed faces with the Type
        //  for these faces. Iterates until nothing changes or maxIter reached.
        //  (maxIter can be 0)
        FvFaceCellWave
        (
            const fvMesh& mesh,
            const List<labelPair>& initialChangedPatchAndFaces,
            const List<Type>& initialChangedFacesInfo,
            List<Type>& internalFaceInfo,
            List<List<Type>>& patchFaceInfo,
            List<Type>& cellInfo,
            const label maxIter,
            TrackingData& td = defaultTrackingData_
        );

        //- Disallow default bitwise copy construction
        FvFaceCellWave(const FvFaceCellWave&) = delete;


    //- Destructor
    virtual ~FvFaceCellWave()
    {}


    // Member Functions

        // Access

            //- Access internalFaceInfo
            List<Type>& internalFaceInfo()
            {
                return internalFaceInfo_;
            }

            //- Access patchFaceInfo
            List<List<Type>>& patchFaceInfo()
            {
                return patchFaceInfo_;
            }

            //- Access cellInfo
            List<Type>& cellInfo()
            {
                return cellInfo_;
            }

            //- Additional data to be passed into container
            const TrackingData& data() const
            {
                return td_;
            }

            //- Additional data to be passed into container
            TrackingData& data()
            {
                return td_;
            }

            //- Access mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }


        // Edit

            //- Set initial changed faces
            void setFaceInfo
            (
                const List<labelPair>& changedPatchAndFaces,
                const List<Type>& changedFacesInfo
            );

            //- Propagate from face to cell. Returns total number of cells
            //  (over all processors) changed.
            virtual label faceToCell();

            //- Propagate from cell to face. Returns total number of faces
            //  (over all processors) changed. (Faces on processorpatches are
            //  counted double)
            virtual label cellToFace();

            //- Iterate until no changes or maxIter reached.  Returns actual
            //  number of iterations.
            virtual label iterate(const label maxIter);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const FvFaceCellWave&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FvFaceCellWave.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
