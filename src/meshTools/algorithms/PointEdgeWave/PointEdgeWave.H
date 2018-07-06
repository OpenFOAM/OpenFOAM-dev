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

Class
    Foam::PointEdgeWave

Description
    Wave propagation of information through grid. Every iteration
    information goes through one layer of edges.

    Templated on information that is transferred.

    Handles parallel and cyclics. Only parallel reasonably tested. Cyclics
    hardly tested.

    Note: whether to propagate depends on the return value of Type::update
    which returns true (i.e. propagate) if the value changes by more than a
    certain tolerance.

    Note: parallel is done in two steps:
      -# transfer patch points in offset notation, i.e. every patch
         point is denoted by a patchface label and an index in this face.
         Receiving end uses that fact that f[0] is shared and order is
         reversed.
      -# do all non-local shared points by means of reduce of data on them.

    Note: cyclics is with offset in patchface as well. Patch is divided into
    two sub patches and the point-point addressing is never explicitly
    calculated but instead use is made of the face-face correspondence.
    (it probably is more efficient to calculate a point-point
    correspondence at the start and then reuse this; task to be done)

SourceFiles
    PointEdgeWave.C

\*---------------------------------------------------------------------------*/

#ifndef PointEdgeWave_H
#define PointEdgeWave_H

#include "boolList.H"
#include "scalarField.H"
#include "tensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyPatch;

/*---------------------------------------------------------------------------*\
                        Class PointEdgeWaveName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(PointEdgeWave);


/*---------------------------------------------------------------------------*\
                           Class PointEdgeWave Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class TrackingData = int>
class PointEdgeWave
:
    public PointEdgeWaveName
{
  // Private static data

        //- Relative tolerance. Stop propagation if relative changes
        //  less than this tolerance (responsibility for checking this is
        //  up to Type implementation)
        static scalar propagationTol_;

        //- Used as default trackdata value to satisfy default template
        //  argument.
        static int dummyTrackData_;


    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Wall information for all points
        UList<Type>& allPointInfo_;

        //- Information on all mesh edges
        UList<Type>& allEdgeInfo_;

        //- Additional data to be passed into container
        TrackingData& td_;

        //- Has point changed
        boolList changedPoint_;

        //- List of changed points
        labelList changedPoints_;

        //- Number of changed points
        label nChangedPoints_;

        //- Edges that have changed
        boolList changedEdge_;
        labelList changedEdges_;
        label nChangedEdges_;

        //- Number of cyclic patches
        label nCyclicPatches_;

        //- Number of evaluations
        label nEvals_;

        //- Number of unvisited edges/points
        label nUnvisitedPoints_;
        label nUnvisitedEdges_;


    // Private Member Functions

        //- Adapt pointInfo for leaving domain
        void leaveDomain
        (
            const polyPatch&,
            const List<label>& patchPointLabels,
            List<Type>& pointInfo
        ) const;

        //- Adapt pointInfo for entering domain
        void enterDomain
        (
            const polyPatch&,
            const List<label>& patchPointLabels,
            List<Type>& pointInfo
        ) const;

        //- Transform. Implementation referred to Type
        void transform
        (
            const polyPatch& patch,
            const tensorField& rotTensor,
            List<Type>& pointInfo
        ) const;

        //- Updates pointInfo with information from neighbour. Updates all
        //  statistics.
        bool updatePoint
        (
            const label pointi,
            const label neighbourEdgeI,
            const Type& neighbourInfo,
            Type& pointInfo
        );

        //- Updates pointInfo with information from same point. Updates all
        //  statistics.
        bool updatePoint
        (
            const label pointi,
            const Type& neighbourInfo,
            Type& pointInfo
        );

        //- Updates edgeInfo with information from neighbour. Updates all
        //  statistics.
        bool updateEdge
        (
            const label edgeI,
            const label neighbourPointi,
            const Type& neighbourInfo,
            Type& edgeInfo
        );

        // Parallel, cyclic

            //- Has patches of certain type?
            template<class PatchType>
            label countPatchType() const;

            //- Merge data from across processor boundaries
            void handleProcPatches();

            //- Merge data from across cyclic boundaries
            void handleCyclicPatches();

            //- Explicitly sync all collocated points
            label handleCollocatedPoints();


        //- Disallow default bitwise copy construct
        PointEdgeWave(const PointEdgeWave&);

        //- Disallow default bitwise assignment
        void operator=(const PointEdgeWave&);


public:

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


    // Constructors

        //- Construct from mesh, list of changed points with the Type
        //  for these points. Gets work arrays to operate on, one of size
        //  number of mesh points, the other number of mesh edges.
        //  Iterates until nothing changes or maxIter reached.
        //  (maxIter can be 0)
        PointEdgeWave
        (
            const polyMesh& mesh,
            const labelList& initialPoints,
            const List<Type>& initialPointsInfo,
            UList<Type>& allPointInfo,
            UList<Type>& allEdgeInfo,
            const label maxIter,
            TrackingData& td = dummyTrackData_
        );

        //- Construct from mesh. Use setPointInfo and iterate() to do
        //  actual calculation
        PointEdgeWave
        (
            const polyMesh& mesh,
            UList<Type>& allPointInfo,
            UList<Type>& allEdgeInfo,
            TrackingData& td = dummyTrackData_
        );


    //- Destructor
    ~PointEdgeWave();


    // Member Functions

        //- Access allPointInfo
        UList<Type>& allPointInfo() const
        {
            return allPointInfo_;
        }

        //- Access allEdgeInfo
        UList<Type>& allEdgeInfo() const
        {
            return allEdgeInfo_;
        }

        //- Additional data to be passed into container
        const TrackingData& data() const
        {
            return td_;
        }

        //- Get number of unvisited edges, i.e. edges that were not (yet)
        //  reached from walking across mesh. This can happen from
        //  - not enough iterations done
        //  - a disconnected mesh
        //  - a mesh without walls in it
        label getUnsetEdges() const;

        label getUnsetPoints() const;

        //- Copy initial data into allPointInfo_
        void setPointInfo
        (
            const labelList& changedPoints,
            const List<Type>& changedPointsInfo
        );

        //- Propagate from point to edge. Returns total number of edges
        //  (over all processors) changed.
        label pointToEdge();

        //- Propagate from edge to point. Returns total number of points
        //  (over all processors) changed.
        label edgeToPoint();

        //- Iterate until no changes or maxIter reached. Returns actual
        //  number of iterations.
        label iterate(const label maxIter);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class listUpdateOp Declaration
\*---------------------------------------------------------------------------*/

//- List update operation
template<class Type, class TrackingData = int>
class listUpdateOp
{
    //- Additional data to be passed into container

    const scalar tol_;

    TrackingData& td_;

public:
    listUpdateOp(const scalar tol, TrackingData& td)
    :
        tol_(tol),
        td_(td)
    {}

    void operator()(List<Type>& x, const List<Type>& y) const
    {
        forAll(x, i)
        {
            if (y[i].valid(td_))
            {
                x[i].updatePoint(y[i], tol_, td_);
            }
        }
    }
};

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PointEdgeWave.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
